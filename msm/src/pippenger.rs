use crate::operations::{add_points, scalar_multiply};
use ark_ff::Zero;
use ark_mnt4_298::{G1Projective, Fr};
use std::collections::HashMap;

// Main pippenger function
pub fn pippenger(points: &[G1Projective], scalars: &[u32], window_size: usize) -> G1Projective {

    let partitions = partition_msm(scalars, window_size);
    combine_partitioned_msm(&partitions, points)
}

struct MsmPartition {
    bit_index: usize,
    window_values: Vec<u32>,
}

/// Step 1: Split MSM with 32-bit scalars into 32/c MSMs with c-bit scalars. c == window_size
fn partition_msm(scalars: &[u32], window_size: usize) -> Vec<MsmPartition> {
    // Calculate the total number of partitions based on window size
    let num_partitions = 32 / window_size;

    // Vector to hold information on partitions
    let mut partitions = Vec::new();

    for partition_index in 0..num_partitions {
        // Calculate the starting bit index for the current partition
        let bit_index = partition_index * window_size;

        // Collect the corresponding window values for each scalar
        let window_values: Vec<u32> = scalars.iter().map(|&scalar| {
            // Shift the scalar to right align the desired bits and mask off the rest
            // Eg: 10110110 (182 in decimal)
            // Bit_index = 0, LHS: shift scalar by 0 so 10110110, RHS: 00000011, so you get 10 by bitwise AND operator
            // Bit index = 2, LHS: shift scalar by 2 so 00101101, RHS: 00000011, so you get 01 by bitwise AND operator
            // Bit index = 4, LHS: shift scalar by 4 so 00001011, RHS: 00000011 so you get 11 by bitwise AND operator
            // Bit index = 6, LHS: shift scalar by 6 so 00000010, RHS: 00000011, so you get 10 by bitwise AND operator
            (scalar >> bit_index) & ((1 << window_size) - 1)
        }).collect();

        // Push the partition information to the list of partitions
        partitions.push(MsmPartition { bit_index, window_values });
    }

    // Return the list of partitions, each containing its bit index and window values
    partitions
}

// Step 2: Compute MSM for each partition
fn compute_msm_for_partition(partition: &MsmPartition, points: &[G1Projective]) -> G1Projective {
    // HashMap to bucket indexes of similar scalar values
    let mut buckets: HashMap<u32, Vec<usize>> = HashMap::new();

    // Iterating over each window value and its index
    for (index, &value) in partition.window_values.iter().enumerate() {
        // Grouping indexes with the same scalar value
        buckets.entry(value).or_insert_with(Vec::new).push(index);
    }

    // Variable to store the computed MSM for this partition
    let mut msm_result = G1Projective::zero();

    // Iterating over each bucket
    for (&value, indexes) in &buckets {
        // Summing up the points corresponding to the indexes in the bucket
        let sum_of_points: G1Projective = indexes.iter()
            .map(|&i| points[i]) // Mapping index to point
            .fold(G1Projective::zero(), |acc, p| add_points(acc, p)); // Summing points
        
        // Adding the summed points to the MSM result after multiplying by the scalar value
        msm_result = add_points(msm_result, scalar_multiply(sum_of_points, Fr::from(value)));
    }

    // Returning the MSM result for this partition
    msm_result
}

// Step 3: Compute the final MSM result by combining all partitions
fn combine_partitioned_msm(partitions: &[MsmPartition], points: &[G1Projective]) -> G1Projective {
    // Variable to store the final MSM result
    let mut final_result = G1Projective::zero();

    // Iterating over each partition
    for partition in partitions {
        // Computing MSM for the current partition
        let partition_msm = compute_msm_for_partition(partition, points);

        // Scaling the partition MSM by 2^bit_index and adding it to the final result
        final_result = add_points(final_result, scalar_multiply(partition_msm, Fr::from(1 << partition.bit_index)));
    }

    // Returning the final combined MSM result
    final_result
}