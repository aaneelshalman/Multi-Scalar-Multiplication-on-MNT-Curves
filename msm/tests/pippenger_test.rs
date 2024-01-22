// design naive method, check if output the same
// check if partition msm returns correctly, give random integer, print bit value, check if partitions match expected value
// check if compute msm works correctly, generate random partitions struct, give random points list, check if its correct
// check if combine msm works correctly, generate random partitions struct, give random points list, check if its correct
// check if these 3 methods output same as naive

use msm::naive::naive_msm;
use ark_mnt4_298::G1Projective;
use ark_std::{test_rng, UniformRand};

#[test]
fn test_naive_msm() {
    let mut rng = test_rng();

    // Generate a list of 10 random points
    let num_points = 10;
    let points: Vec<G1Projective> = (0..num_points)
        .map(|_| G1Projective::rand(&mut rng))
        .collect();

    // Create a list of scalars with first element being 1 and the rest 0
    let mut scalars = vec![0; num_points];
    scalars[0] = 1;

    let msm_result = naive_msm(&points, &scalars);

    // The result of MSM should be equal to the first point
    assert_eq!(msm_result, points[0], "MSM result does not match the first point");
}
