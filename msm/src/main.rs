extern crate ark_ff;
extern crate ark_ec;
extern crate ark_mnt4_298;
extern crate ark_std;

mod operations;
mod pippenger;
mod parallel_pippenger;
mod naive;

use ark_mnt4_298::G1Projective;
use ark_std::{UniformRand, test_rng};
// use operations::{add_points, scalar_multiply};

fn main() {
    let mut rng = test_rng();

    let points = vec![
        G1Projective::rand(&mut rng),
        G1Projective::rand(&mut rng),
        G1Projective::rand(&mut rng),
        G1Projective::rand(&mut rng),
        G1Projective::rand(&mut rng),
    ];

    let scalars = vec![2, 4, 6, 8, 10];

    let window_size = 2;

    let msm_result = naive::naive_msm(&points, &scalars);
    
    let msm_result_2 = pippenger::pippenger(&points, &scalars, window_size);

    let msm_result_3 = parallel_pippenger::parallel_pippenger(&points, &scalars, window_size);

    assert_eq!(msm_result, msm_result_2);
    assert_eq!(msm_result, msm_result_3);
}