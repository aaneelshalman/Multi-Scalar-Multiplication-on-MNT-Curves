use crate::operations::add_points;
use ark_ec::Group;
use ark_ff::Zero;
use ark_mnt4_298::G1Projective;
use std::collections::BTreeMap;

// Main pippenger function
pub fn subsum_pippenger(points: &[G1Projective], scalars: &[u32], window_size: usize) -> G1Projective {

    // Ensure points and scalars have the same length
    assert_eq!(points.len(), scalars.len(), "Points and scalars must have the same length");
    
    let partitions = subsum_partition_msm(scalars, window_size);
    subsum_combine_partitioned_msm(&partitions, points, window_size)
}

pub struct SubsumMsmPartition {
    pub bit_index: usize,
    pub window_values: Vec<u32>,
}

// Step 1: Split MSM with 32-bit scalars into 32/c MSMs with c-bit scalars. c == window_size
pub fn subsum_partition_msm(scalars: &[u32], window_size: usize) -> Vec<SubsumMsmPartition> {
    
    // Calculate the total number of partitions based on window size
    // (32 + window_size - 1) / window_size is used so that if 32 divides window_size, num_partitions will return the divisor
    // But if 32 does not divide window_size, then num_partitions will round up the float instead of round down, the default in Rust
    let num_partitions = (32 + window_size - 1) / window_size;

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
        partitions.push(SubsumMsmPartition { bit_index, window_values });
    }

    // Return the list of partitions, each containing its bit index and window values
    partitions
}

// Step 2: Compute MSM for each partition using subsum accumulation
pub fn subsum_compute_msm_for_partition(partition: &SubsumMsmPartition, points: &[G1Projective], window_size: usize) -> G1Projective {
    let mut buckets: BTreeMap<u32, Vec<usize>> = BTreeMap::new();

    // Populate buckets with indexes grouped by their scalar value
    for (index, &value) in partition.window_values.iter().enumerate() {
        if value != 0 {
            buckets.entry(value).or_insert_with(Vec::new).push(index);
        }
    }

    // Determine the maximum scalar value (m) for array initialization
    let m = (1 << window_size) - 1;
    let mut tmp = vec![G1Projective::zero(); m + 1];

    // Next scalar placeholder for comparison
    let mut next_scalar: Option<u32> = None;

    // Reverse iterate through the sorted buckets
    for (&scalar, indexes) in buckets.iter().rev() {
        let sum_of_points = indexes.iter()
            .map(|&index| points[index])
            .fold(G1Projective::zero(), |acc, point| add_points(acc, point));

        // Place the sum of points in tmp[0] initially
        tmp[0] = add_points(tmp[0], sum_of_points);

        if let Some(next_scalar_val) = next_scalar {
            // Calculate k as the difference between the current scalar and the next scalar
            let k = (scalar - next_scalar_val) as usize;
            if k < tmp.len() {
                tmp[k] = add_points(tmp[k], tmp[0]);
            }
        }
        // Update next_scalar for the next iteration
        next_scalar = Some(scalar);
    }

    // Initialize temp_point and msm_result to zero points
    let mut temp_point = G1Projective::zero();
    let mut msm_result = G1Projective::zero();

    // Accumulate subsums from tmp array
    for i in (0..=m).rev() {
        temp_point = add_points(temp_point, tmp[i]);
        msm_result = add_points(msm_result, temp_point);
    }

    msm_result
}


// Step 3: Compute the final MSM result by combining all partitions
pub fn subsum_combine_partitioned_msm(partitions: &[SubsumMsmPartition], points: &[G1Projective], window_size: usize) -> G1Projective {
    // Variable to store the final MSM result
    let mut final_result = G1Projective::zero();

    // Iterating over each partition in reverse to ensure doubling mimics scaling accurately
    for partition in partitions.iter().rev() {
        // Computing MSM for the current partition
        let partition_msm = subsum_compute_msm_for_partition(partition, points, window_size);

        // Double the final result window_size times to mimic scaling by 2^bit_index
        for _ in 0..window_size {
            final_result = final_result.double();
        }

        // Adding the partition MSM to the final result
        final_result = add_points(final_result, partition_msm);
    }

    // Returning the final combined MSM result
    final_result
}