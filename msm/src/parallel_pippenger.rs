use crate::operations::{add_points, scalar_multiply};
use ark_ff::Zero;
use ark_mnt4_298::{G1Projective, Fr};
use std::collections::HashMap;
use std::thread;

// Main function for Pippenger with parallelism
pub fn parallel_pippenger(points: &[G1Projective], scalars: &[u32], window_size: usize) -> G1Projective {
    assert_eq!(points.len(), scalars.len(), "Points and scalars must have the same length");
    
    let partitions = parallel_partition_msm(scalars, window_size);
    parallel_combine_partitioned_msm(&partitions, points)
}

pub struct ParallelMsmPartition {
    pub bit_index: usize,
    pub window_values: Vec<u32>,
}

impl Clone for ParallelMsmPartition {
    fn clone(&self) -> ParallelMsmPartition {
        ParallelMsmPartition {
            bit_index: self.bit_index,
            window_values: self.window_values.clone(),
        }
    }
}

pub fn parallel_partition_msm(scalars: &[u32], window_size: usize) -> Vec<ParallelMsmPartition> {
    let num_partitions = (32 + window_size - 1) / window_size;
    let mut partitions = Vec::new();

    for partition_index in 0..num_partitions {
        let bit_index = partition_index * window_size;
        let window_values: Vec<u32> = scalars.iter().map(|&scalar| {
            (scalar >> bit_index) & ((1 << window_size) - 1)
        }).collect();
        partitions.push(ParallelMsmPartition { bit_index, window_values });
    }

    partitions
}

pub fn parallel_compute_msm_for_partition(partition: &ParallelMsmPartition, points: &[G1Projective]) -> G1Projective {
    let mut buckets: HashMap<u32, Vec<usize>> = HashMap::new();
    for (index, &value) in partition.window_values.iter().enumerate() {
        if value != 0 {
        buckets.entry(value).or_insert_with(Vec::new).push(index);
        }
    }

    let mut msm_result = G1Projective::zero();
    for (&value, indexes) in &buckets {
        let sum_of_points: G1Projective = indexes.iter()
            .map(|&i| points[i])
            .fold(G1Projective::zero(), |acc, p| add_points(acc, p));
        msm_result = add_points(msm_result, scalar_multiply(sum_of_points, Fr::from(value)));
    }

    msm_result
}


pub fn parallel_combine_partitioned_msm(partitions: &[ParallelMsmPartition], points: &[G1Projective]) -> G1Projective {
    let mut handles = Vec::new();

    // Spawn a thread for each partition
    for partition in partitions {
        let partition_clone = partition.clone();
        let points_clone = points.to_vec();
        let bit_index = partition.bit_index;

        let handle = thread::spawn(move || {
            let msm_result = parallel_compute_msm_for_partition(&partition_clone, &points_clone);
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
