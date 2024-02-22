use crate::operations::{add_points, scalar_multiply};
use ark_ec::Group;
use ark_ff::Zero;
use ark_mnt4_298::{G1Projective, Fr};
use std::collections::HashMap;
use std::ops::Neg;

// Main function for Pippenger with 2-NAF Decomposition
pub fn naf_pippenger(points: &[G1Projective], scalars: &[u32], window_size: usize) -> G1Projective {
    assert_eq!(points.len(), scalars.len(), "Points and scalars must have the same length");
    
    let partitions = naf_partition_msm(scalars, window_size);
    let decomposed_partitions = naf_decompose_partitions(&partitions, window_size);
    naf_combine_partitioned_msm(&decomposed_partitions, points)
}

pub struct NafMsmPartition {
    pub bit_index: usize,
    pub window_values: Vec<u32>,
}

impl Clone for NafMsmPartition {
    fn clone(&self) -> NafMsmPartition {
        NafMsmPartition {
            bit_index: self.bit_index,
            window_values: self.window_values.clone(),
        }
    }
}

pub struct NafMsmPartitionDecomposed {
    pub bit_index: usize,
    pub window_values: Vec<i64>, // New struct adjusted for signed integer decomposition
}


pub fn naf_partition_msm(scalars: &[u32], window_size: usize) -> Vec<NafMsmPartition> {
    let num_partitions = (32 + window_size - 1) / window_size;
    let mut partitions = Vec::new();

    for partition_index in 0..num_partitions {
        let bit_index = partition_index * window_size;
        let window_values: Vec<u32> = scalars.iter().map(|&scalar| {
            (scalar >> bit_index) & ((1 << window_size) - 1)
        }).collect();
        partitions.push(NafMsmPartition { bit_index, window_values });
    }

    partitions
}

pub fn naf_decompose_partitions(partitions: &[NafMsmPartition], window_size: usize) -> Vec<NafMsmPartitionDecomposed> {
    let base = 2u32.pow(window_size as u32); // 2^(window_size)
    let threshold = base / 2; // 2^(window_size-1)

    // Initialize decomposed partitions with the same structure but empty window values
    let mut decomposed_partitions: Vec<NafMsmPartitionDecomposed> = partitions.iter()
        .map(|p| NafMsmPartitionDecomposed {
            bit_index: p.bit_index,
            window_values: vec![0i64; p.window_values.len()], // Initialize with zeros
        })
        .collect();

    // Append an extra partition for overflow handling
    let last_bit_index = partitions.last().unwrap().bit_index + window_size;
    decomposed_partitions.push(NafMsmPartitionDecomposed {
        bit_index: last_bit_index,
        window_values: vec![0i64; partitions[0].window_values.len()], // Initialize with zeros for overflow handling
    });

    // Iterate through each position of window values across all partitions
    for i in 0..partitions[0].window_values.len() {
        let mut carry = 0i64;

        for j in 0..partitions.len() {
            let window_value = partitions[j].window_values[i] as i64;
            let adjusted_value = window_value + carry;
            carry = 0; // Reset carry for the next window value

            if adjusted_value >= threshold as i64 {
                decomposed_partitions[j].window_values[i] = adjusted_value - base as i64;
                // Ensure carry is forwarded to the next partition's same position
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

pub fn naf_compute_msm_for_partition(partition: &NafMsmPartitionDecomposed, points: &[G1Projective]) -> G1Projective {
    let mut buckets: HashMap<u32, Vec<(usize, i64)>> = HashMap::new(); // Use magnitude for keys and keep sign with index for values
    
    // Assign points to buckets based on the absolute value while keeping track of the original value's sign
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


pub fn naf_combine_partitioned_msm(partitions: &[NafMsmPartitionDecomposed], points: &[G1Projective]) -> G1Projective {
    partitions.iter().fold(G1Projective::zero(), |acc, partition| {
        let mut partition_msm = naf_compute_msm_for_partition(partition, points);

        // Iteratively double the partition result bit_index times
        for _ in 0..partition.bit_index {
            partition_msm = partition_msm.double();
        }

        // Add the iteratively doubled result to the accumulated result
        add_points(acc, partition_msm)
    })
}