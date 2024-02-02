use crate::operations::{add_points, scalar_multiply};
use ark_ff::Zero;
use ark_mnt4_298::{G1Projective, Fr};
use std::collections::HashMap;
use std::thread;

// Adjusted to take &[u32] for scalars
pub fn parallel_naf_pippenger(points: &[G1Projective], scalars: &[u32], window_size: usize) -> G1Projective {
    assert_eq!(points.len(), scalars.len(), "Points and scalars must have the same length");
    
    let partitions = parallel_naf_partition_msm(scalars, window_size);
    parallel_naf_combine_partitioned_msm(&partitions, points)
}

pub struct ParallelNafMsmPartition {
    pub bit_index: usize,
    pub window_values: Vec<i64>,
}

// Adjusted to take u32 scalar, perform NAF without needing to convert to i64
pub fn to_2_naf(window_size: usize, scalar: u32) -> Vec<i64> {
    let mut naf = Vec::new();
    let mut temp_scalar = scalar as i64; // Convert to i64 if necessary for signed operations
    let base = 1 << window_size;

    while temp_scalar != 0 {
        if temp_scalar & (base / 2) != 0 {
            let adjustment = base * if temp_scalar < 0 { 1 } else { -1 };
            naf.push(temp_scalar - adjustment);
            temp_scalar += adjustment;
        } else {
            naf.push(0);
        }
        temp_scalar >>= window_size;
    }

    naf
}

// Adjusted to take &[u32] for scalars
pub fn parallel_naf_partition_msm(scalars: &[u32], window_size: usize) -> Vec<ParallelNafMsmPartition> {
    let num_partitions = (32 + window_size - 1) / window_size;
    let mut partitions = vec![ParallelNafMsmPartition { bit_index: 0, window_values: vec![0; scalars.len()] }; num_partitions];

    for (scalar_idx, &scalar) in scalars.iter().enumerate() {
        let naf = to_2_naf(window_size, scalar);
        for (partition_index, &naf_value) in naf.iter().enumerate() {
            if partition_index < num_partitions {
                partitions[partition_index].window_values[scalar_idx] = naf_value;
            }
        }
    }

    partitions
}

pub fn parallel_naf_compute_msm_for_partition(partition: &ParallelNafMsmPartition, points: &[G1Projective]) -> G1Projective {
    let mut buckets: HashMap<i64, Vec<usize>> = HashMap::new();
    for (index, &value) in partition.window_values.iter().enumerate() {
        if value != 0 {
            let bucket = if value > 0 { value as u32 } else { (-value) as u32 };
            buckets.entry(value).or_insert_with(Vec::new).push(index);
        }
    }

    let mut msm_result = G1Projective::zero();
    for (&value, indexes) in &buckets {
        let sum_of_points: G1Projective = indexes.iter()
            .map(|&i| if value > 0 { points[i] } else { -points[i] })
            .fold(G1Projective::zero(), |acc, p| add_points(acc, p));
        msm_result = add_points(msm_result, scalar_multiply(sum_of_points, Fr::from(value.abs() as u64)));
    }

    msm_result
}

pub fn parallel_naf_combine_partitioned_msm(partitions: &[ParallelNafMsmPartition], points: &[G1Projective]) -> G1Projective {
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

impl Clone for ParallelNafMsmPartition {
    fn clone(&self) -> ParallelNafMsmPartition {
        ParallelNafMsmPartition {
            bit_index: self.bit_index,
            window_values: self.window_values.clone(),
        }
    }
}
