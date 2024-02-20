extern crate ark_ff;
extern crate ark_ec;
extern crate ark_mnt4_298;
extern crate ark_std;
extern crate std;

mod operations;
mod pippenger;
mod parallel_pippenger;
mod naive;
mod parallel_naf_pippenger;
mod naf_pippenger;

use ark_mnt4_298::G1Projective;
use ark_std::{UniformRand, test_rng};
use rand::{Rng, thread_rng};
use std::time::Instant;
use parallel_pippenger::parallel_pippenger;
use pippenger::pippenger;
use naive::naive_msm;
use naf_pippenger::naf_pippenger;
use parallel_naf_pippenger::parallel_naf_pippenger;
// use operations::{add_points, scalar_multiply};

fn main() {

    fn generate_points(num_points: usize) -> Vec<G1Projective> {
        let mut rng = test_rng();
        (0..num_points).map(|_| G1Projective::rand(&mut rng)).collect()
    }
    
    fn generate_scalars(num_scalars: usize) -> Vec<u32> {
        let mut rng = thread_rng();
        (0..num_scalars).map(|_| rng.gen_range(0..65536)).collect()
    }

    let points = generate_points(1000);
    let scalars = generate_scalars(1000);
    let window_size = 2;

    let start = Instant::now();
    let result_naive = naive_msm(&points, &scalars);
    let duration_naive = start.elapsed();

    let start = Instant::now();
    let result_pippenger = pippenger(&points, &scalars, window_size);
    let duration_pippenger = start.elapsed();

    let start = Instant::now();
    let result_parallel = parallel_pippenger(&points, &scalars, window_size);
    let duration_parallel = start.elapsed();

    let start = Instant::now();
    let result_naf = naf_pippenger(&points, &scalars, window_size);
    let duration_naf = start.elapsed();

    let start = Instant::now();
    let result_parallel_naf = parallel_naf_pippenger(&points, &scalars, window_size);
    let duration_parallel_naf = start.elapsed();

    assert_eq!(result_naive, result_pippenger, "Results of pippenger and naive MSM should match");
    assert_eq!(result_naive, result_parallel, "Results of pippenger with parallelism and naive MSM should match");
    assert_eq!(result_naive, result_naf, "Results of pippenger with 2-NAF Decomposition should match with naive MSM");
    assert_eq!(result_naive, result_parallel_naf, "Results of pippenger with 2-NAF Decomposition and parallelism should match with naive MSM");

    println!("Naive: {:?}", duration_naive);
    println!("Pippenger: {:?}", duration_pippenger);
    println!("Parallel: {:?}", duration_parallel);
    println!("NAF: {:?}", duration_naf);
    println!("Parallel NAF: {:?}", duration_parallel_naf);
}