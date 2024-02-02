use crate::operations::{add_points, scalar_multiply};
use ark_ff::Zero;
use ark_mnt4_298::{G1Projective, Fr};
use std::collections::HashMap;
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
    let threshold = 2u32.pow(window_size as u32 - 1);
    
    let decomposed_partitions: Vec<ParallelNafMsmPartitionDecomposed> = partitions.iter().map(|partition| {
        let mut decomposed_window_values: Vec<i64> = Vec::with_capacity(partition.window_values.len() + 1); // +1 in case we need to append an extra value for the carry

        let mut carry: i64 = 0;

        for &value in partition.window_values.iter() {
            let adjusted_value = value as i64 + carry;
            carry = 0; // Reset carry for each window_value

            if adjusted_value >= threshold as i64 {
                carry = 1; // Indicate a carry for the next value
                decomposed_window_values.push(adjusted_value - base as i64);
            } else {
                decomposed_window_values.push(adjusted_value);
            }
        }

        // Handle outstanding carry by appending an extra value if necessary
        if carry == 1 {
            decomposed_window_values.push(1); // Append '1' to represent the final carry
        }

        ParallelNafMsmPartitionDecomposed {
            bit_index: partition.bit_index,
            window_values: decomposed_window_values,
        }
    }).collect();

    decomposed_partitions
}

pub fn parallel_naf_compute_msm_for_partition(partition: &ParallelNafMsmPartitionDecomposed, points: &[G1Projective]) -> G1Projective {
    let mut buckets: HashMap<u32, Vec<(usize, i64)>> = HashMap::new(); // Use magnitude for keys and keep sign with index
    
    // Assign points to buckets based on the absolute value of their window value, but keep track of the original value's sign
    for (index, &value) in partition.window_values.iter().enumerate() {
        let abs_value = value.abs() as u32; // Use absolute value for bucket key
        buckets.entry(abs_value).or_insert_with(Vec::new).push((index, value));
    }

    let mut msm_result = G1Projective::zero();
    // Iterate over buckets to compute the MSM contribution of each
    for (&abs_value, index_sign_pairs) in &buckets {
        let sum_of_points: G1Projective = index_sign_pairs.iter()
            .map(|&(i, sign)| {
                let mut point = points[i];
                if sign < 0 { point = -point; } // Negate the point if the original value was negative
                point
            })
            .fold(G1Projective::zero(), |acc, p| add_points(acc, p));
        // Multiply the aggregated point by its scalar and add to the result
        msm_result = add_points(msm_result, scalar_multiply(sum_of_points, Fr::from(abs_value as u64)));
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
        final_result = add_points(final_result, scalar_multiply(partition_result, Fr::from(1 << bit_index)));
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
