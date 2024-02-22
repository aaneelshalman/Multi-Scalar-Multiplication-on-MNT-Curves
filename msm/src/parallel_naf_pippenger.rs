use crate::operations::{add_points, scalar_multiply};
use ark_ec::Group;
use ark_ff::Zero;
use ark_mnt4_298::{G1Projective, Fr};
use std::collections::HashMap;
use std::ops::Neg;
use std::thread;

// Main function for Pippenger with parallelism and 2-NAF Decomposition
pub fn parallel_naf_pippenger(points: &[G1Projective], scalars: &[u32], window_size: usize) -> G1Projective {
    assert_eq!(points.len(), scalars.len(), "Points and scalars must have the same length");
    
    let partitions = parallel_naf_partition_msm(scalars, window_size);
    let decomposed_partitions = parallel_naf_decompose_partitions(&partitions, window_size);
    parallel_naf_combine_partitioned_msm(&decomposed_partitions, points)
}
pub struct ParallelNafMsmPartition {
    pub bit_index: usize,
    pub window_values: Vec<u32>,
}

impl Clone for ParallelNafMsmPartition {
    fn clone(&self) -> ParallelNafMsmPartition {
        ParallelNafMsmPartition {
            bit_index: self.bit_index,
            window_values: self.window_values.clone(),
        }
    }
}
pub struct ParallelNafMsmPartitionDecomposed {
    pub bit_index: usize,
    pub window_values: Vec<i64>, // New struct adjusted for signed integer decomposition
}

impl Clone for ParallelNafMsmPartitionDecomposed {
    fn clone(&self) -> ParallelNafMsmPartitionDecomposed {
        ParallelNafMsmPartitionDecomposed {
            bit_index: self.bit_index,
            window_values: self.window_values.clone(),
        }
    }
}


pub fn parallel_naf_partition_msm(scalars: &[u32], window_size: usize) -> Vec<ParallelNafMsmPartition> {
    let num_partitions = (32 + window_size - 1) / window_size;
    let mut partitions = Vec::new();

    for partition_index in 0..num_partitions {
        let bit_index = partition_index * window_size;
        let window_values: Vec<u32> = scalars.iter().map(|&scalar| {
            (scalar >> bit_index) & ((1 << window_size) - 1)
        }).collect();
        partitions.push(ParallelNafMsmPartition { bit_index, window_values });
    }

    partitions
}

pub fn parallel_naf_decompose_partitions(partitions: &[ParallelNafMsmPartition], window_size: usize) -> Vec<ParallelNafMsmPartitionDecomposed> {
    let base = 2u32.pow(window_size as u32); // 2^(window_size)
    let threshold = base / 2; // 2^(window_size-1)

    // Initialize decomposed partitions with the same structure but empty window values
    let mut decomposed_partitions: Vec<ParallelNafMsmPartitionDecomposed> = partitions.iter()
        .map(|p| ParallelNafMsmPartitionDecomposed {
            bit_index: p.bit_index,
            window_values: vec![0i64; p.window_values.len()], // Initialize with zeros
        })
        .collect();

    // Append an extra partition for overflow handling
    let last_bit_index = partitions.last().unwrap().bit_index + window_size;
    decomposed_partitions.push(ParallelNafMsmPartitionDecomposed {
        bit_index: last_bit_index,
        window_values: vec![0i64; partitions[0].window_values.len()], // Initialize with zeros for overflow handling
    });

    // Iterate through each position of window values across all partitions, not window_values in the same bit_index
    for i in 0..partitions[0].window_values.len() {
        let mut carry = 0i64;

        for j in 0..partitions.len() {
            let window_value = partitions[j].window_values[i] as i64;
            let adjusted_value = window_value + carry;
            carry = 0; // Reset carry for the next window value

            if adjusted_value >= threshold as i64 {
                decomposed_partitions[j].window_values[i] = adjusted_value - base as i64;
                // Forward carry to the next partition's same position
                if j < partitions.len() - 1 {
                    carry = 1;
                }
            } else {
                decomposed_partitions[j].window_values[i] = adjusted_value;
            }
        }

        if carry > 0 {
            decomposed_partitions[partitions.len()].window_values[i] += carry;
        }
    }

    decomposed_partitions
}

pub fn parallel_naf_compute_msm_for_partition(partition: &ParallelNafMsmPartitionDecomposed, points: &[G1Projective]) -> G1Projective {
    let mut buckets: HashMap<u32, Vec<(usize, i64)>> = HashMap::new(); // Use magnitude for keys and keep sign with index
    
    // Assign points to buckets based on the absolute value of their window value, but keep track of the original value's sign
    for (index, &value) in partition.window_values.iter().enumerate() {
        if value != 0 {
            let abs_value = value.abs() as u32; // Use absolute value for bucket key
            buckets.entry(abs_value).or_insert_with(Vec::new).push((index, value));
        }
    }

    let mut msm_result = G1Projective::zero();
    for (&abs_value, index_sign_pairs) in &buckets {
        let sum_of_points: G1Projective = index_sign_pairs.iter()
            .map(|&(i, sign)| {
                if i >= points.len() {
                    return G1Projective::zero(); // Return a zero point to avoid panic
                }

                let mut point = points[i];
                if sign < 0 {
                    point = point.neg(); // Negate the point if the original value was negative
                }
                point
            })
            .fold(G1Projective::zero(), |acc, p| add_points(acc, p));

        msm_result = add_points(msm_result, scalar_multiply(sum_of_points, Fr::from(abs_value as u32)));
    }

    msm_result
}

pub fn parallel_naf_combine_partitioned_msm(partitions: &[ParallelNafMsmPartitionDecomposed], points: &[G1Projective]) -> G1Projective {
    let mut handles = Vec::new();

    // Spawn a thread for each partition
    for partition in partitions {
        let partition_clone = partition.clone();
        let points_clone = points.to_vec();

        let handle = thread::spawn(move || {
            let msm_result = parallel_naf_compute_msm_for_partition(&partition_clone, &points_clone);
            (msm_result, partition_clone.bit_index)
        });

        handles.push(handle);
    }

    // Collect results from each thread and combine
    let mut final_result = G1Projective::zero();
    for handle in handles {
        let (mut partition_result, bit_index) = handle.join().unwrap();
        
        // Iteratively double the partition result bit_index times
        for _ in 0..bit_index {
            partition_result = partition_result.double();
        }

        // Add the iteratively doubled result to the final result
        final_result = add_points(final_result, partition_result);
    }

    final_result
}