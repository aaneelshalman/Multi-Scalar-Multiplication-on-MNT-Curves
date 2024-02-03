use crate::operations::{add_points, scalar_multiply};
use ::ark_ff::Field;
use ark_ff::Zero;
use ark_mnt4_298::{G1Projective, Fr};
use std::collections::HashMap;
use std::ops::Neg;
use std::thread;

// Main function for Pippenger with parallelism
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
    let base = 2u32.pow(window_size as u32);
    let threshold = base / 2; // Equivalent to 2^(window_size-1)

    // Extend partitions to include an extra bit index for overflow handling
    let mut extended_partitions = partitions.to_vec();
    let last_bit_index = extended_partitions.last().map_or(0, |p| p.bit_index + window_size);
    extended_partitions.push(ParallelNafMsmPartition {
        bit_index: last_bit_index,
        window_values: vec![0; partitions[0].window_values.len()], // Initialize with zeros
    });

    extended_partitions.iter().map(|partition| {
        let mut decomposed_window_values = vec![0i64; partition.window_values.len()]; // Pre-allocate space for decomposed values

        // Iterate through window values correctly for signed integer decomposition
        for (i, &window_value) in partition.window_values.iter().enumerate() {
            let signed_value = if window_value >= threshold {
                // Adjust values above threshold by subtracting 'base' and adding 1 to the next value if not the last
                if i < partition.window_values.len() - 1 {
                    decomposed_window_values[i + 1] += 1; // Add carry to the next value
                }
                window_value as i64 - base as i64
            } else {
                window_value as i64
            };
            decomposed_window_values[i] = signed_value;
        }

        ParallelNafMsmPartitionDecomposed {
            bit_index: partition.bit_index,
            window_values: decomposed_window_values,
        }
    }).collect()
}


pub fn parallel_naf_compute_msm_for_partition(partition: &ParallelNafMsmPartitionDecomposed, points: &[G1Projective]) -> G1Projective {
    let mut buckets: HashMap<u32, Vec<(usize, i64)>> = HashMap::new(); // Use magnitude for keys and keep sign with index
    
    // Assign points to buckets based on the absolute value of their window value, but keep track of the original value's sign
    for (index, &value) in partition.window_values.iter().enumerate() {
        if value !=0 {
            let abs_value = value.abs() as u32; // Use absolute value for bucket key
            buckets.entry(abs_value).or_insert_with(Vec::new).push((index, value));
        }    
    }

    let mut msm_result = G1Projective::zero();
    for (&abs_value, index_sign_pairs) in &buckets {
        let sum_of_points: G1Projective = index_sign_pairs.iter()
            .map(|&(i, sign)| {
                let mut point = points[i];
                if sign < 0 { point = point.neg(); } // Negate the point if the original value was negative
                point
            })
            .fold(G1Projective::zero(), |acc, p| add_points(acc, p));
        // Multiply the aggregated point by its scalar and add to the result
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
        let bit_index = partition.bit_index;

        let handle = thread::spawn(move || {
            let msm_result = parallel_naf_compute_msm_for_partition(&partition_clone, &points_clone);
            (msm_result, bit_index)
        });

        handles.push(handle);
    }

    // Collect results from each thread and combine
    let mut final_result = G1Projective::zero();
    for handle in handles {
    let (partition_result, bit_index) = handle.join().unwrap();
    // Convert bit_index to 2^bit_index within the scalar field directly
    let exponentiation_result = Fr::from(2).pow(&[bit_index as u64, 0, 0, 0]); // Assuming Fr supports pow
    final_result = add_points(final_result, scalar_multiply(partition_result, exponentiation_result));
}

    final_result
}

impl Clone for ParallelNafMsmPartitionDecomposed {
    fn clone(&self) -> ParallelNafMsmPartitionDecomposed {
        ParallelNafMsmPartitionDecomposed {
            bit_index: self.bit_index,
            window_values: self.window_values.clone(),
        }
    }
}
