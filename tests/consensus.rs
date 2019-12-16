#![cfg(feature = "consensus")]
use approx::assert_relative_eq;
use arrsac::{Arrsac, Config};
use nalgebra::{Isometry3, Point3, Translation, UnitQuaternion, Vector3};
use p3p::nordberg::{NordbergEstimator, Sample};
use rand::{rngs::SmallRng, SeedableRng};
use sample_consensus::Consensus;

const EPSILON_APPROX: f32 = 1e-2;

fn sample_conv(point: [f32; 3], bearing: [f32; 3]) -> Sample {
    Sample { point, bearing }
}

#[test]
fn arrsac() {
    let mut arrsac = Arrsac::new(Config::new(0.01), SmallRng::from_seed([0; 16]));

    // Define some points in camera coordinates (with z > 0).
    let p1_cam = [-0.228_125, -0.061_458_334, 1.0];
    let p2_cam = [0.418_75, -0.581_25, 2.0];
    let p3_cam = [1.128_125, 0.878_125, 3.0];
    let p4_cam = [-0.128_125, 0.578_125, 1.5];

    // Define the camera pose.
    let rot = UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3);
    let trans = Translation::from(Vector3::new(0.1, 0.2, 0.3));
    let pose = Isometry3::from_parts(trans, rot);

    // Compute world coordinates.
    let p1_world = (pose.inverse() * Point3::from(p1_cam)).coords.into();
    let p2_world = (pose.inverse() * Point3::from(p2_cam)).coords.into();
    let p3_world = (pose.inverse() * Point3::from(p3_cam)).coords.into();
    let p4_world = (pose.inverse() * Point3::from(p4_cam)).coords.into();

    let samples = [
        sample_conv(p1_world, p1_cam),
        sample_conv(p2_world, p2_cam),
        sample_conv(p3_world, p3_cam),
        sample_conv(p4_world, p4_cam),
    ];

    // Estimate potential poses with P3P.
    // Arrsac should use the fourth point to filter and find only one model from the 4 generated.
    let pose = arrsac
        .model(&NordbergEstimator, &samples)
        .unwrap()
        .to_iso3();

    // Compare the pose to ground truth.
    assert_relative_eq!(rot, pose.rotation, epsilon = EPSILON_APPROX);
    assert_relative_eq!(trans, pose.translation, epsilon = EPSILON_APPROX);
}
