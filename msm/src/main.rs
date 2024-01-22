extern crate ark_ff;
extern crate ark_ec;
extern crate ark_mnt4_298;
extern crate ark_std;

mod operations;
mod pippenger;

use ark_mnt4_298::{Fr, G1Projective};
use ark_std::{UniformRand, test_rng};
use operations::{add_points, scalar_multiply};

fn main() {
    let mut rng = test_rng();
    
    // Initialise two points on the MNT4-298 curve
    let point1 = G1Projective::rand(&mut rng);
    let point2 = G1Projective::rand(&mut rng);

    println!("Point 1: {:?}", point1);
    println!("Point 2: {:?}", point2);

    let point_sum = add_points(point1, point2);

    let scalar = Fr::from(123456789u128);
    let point_scalar_mul = scalar_multiply(point1, scalar);

    println!("Point Addition Result: {:?}", point_sum);
    println!("Scalar Multiplication Result: {:?}", point_scalar_mul);

    let points = vec![
        point1,
        point2,
        G1Projective::rand(&mut rng),
        G1Projective::rand(&mut rng),
        G1Projective::rand(&mut rng),
    ];

    println!("Point 0: {:?}", points[0]);
    println!("Point 1: {:?}", points[1]);
    println!("Point 2: {:?}", points[2 ]);
    println!("Point 3: {:?}", points[3]);
    println!("Point 4: {:?}", points[4]);

    let scalars = vec![2, 4, 6, 8, 10];

    let window_size = 2;

    let msm_result = pippenger::pippenger(&points, &scalars, window_size);

    println!("MSM Result: {:?}", msm_result);
}