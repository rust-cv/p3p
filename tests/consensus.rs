#![cfg(feature = "consensus")]
use approx::assert_relative_eq;
use arraymap::ArrayMap;
use arrsac::{Arrsac, Config};
use nalgebra::{Isometry3, Point3, Translation, UnitQuaternion, Vector3};
use p3p::nordberg::{NordbergEstimator, Sample};
use rand::{rngs::SmallRng, SeedableRng};
use sample_consensus::Consensus;

const EPSILON_APPROX: f32 = 1e-2;

fn sample_conv(world: [f32; 3], camera: [f32; 2]) -> Sample {
    Sample { world, camera }
}

#[test]
fn arrsac_manual() {
    let mut arrsac = Arrsac::new(Config::new(0.01), SmallRng::from_seed([0; 16]));

    // Define some points in camera coordinates (with z > 0).
    let camera_depth_points = [
        [-0.228_125, -0.061_458_334, 1.0],
        [0.418_75, -0.581_25, 2.0],
        [1.128_125, 0.878_125, 3.0],
        [-0.528_125, 0.178_125, 2.5],
    ]
    .map(|&p| Point3::from(p));

    // Define the camera pose.
    let rot = UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3);
    let trans = Translation::from(Vector3::new(0.1, 0.2, 0.3));
    let pose = Isometry3::from_parts(trans, rot);

    // Compute world coordinates.
    let world_points = camera_depth_points.map(|p| pose.inverse() * p);

    // Compute normalized image coordinates.
    let normalized_image_coordinates = camera_depth_points.map(|p| (p / p.z).xy());

    let samples: Vec<Sample> = world_points
        .iter()
        .zip(&normalized_image_coordinates)
        .map(|(world, image)| sample_conv(world.coords.into(), image.coords.into()))
        .collect();

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
