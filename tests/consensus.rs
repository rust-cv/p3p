#![cfg(feature = "consensus")]
use approx::{assert_relative_eq, relative_eq};
use arrsac::{Arrsac, Config};
use itertools::Itertools;
use nalgebra::{Isometry3, Point3, Translation, UnitQuaternion, Vector2, Vector3};
use p3p::nordberg::{NordbergEstimator, Sample};
use quickcheck_macros::quickcheck;
use rand::{rngs::SmallRng, SeedableRng};
use sample_consensus::Consensus;

const EPSILON_APPROX: f32 = 1e-2;

fn sample_conv(world: [f32; 3], camera: [f32; 3]) -> Sample {
    Sample { world, camera }
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

type V3 = (f32, f32, f32);

/// This test is ignored because it is random and may fail in CI.
/// Run `cargo test -- --ignored` to test it.
#[quickcheck]
fn non_degenerate_case(rot: V3, trans: V3, p1: V3, p2: V3, p3: V3, p4: V3) -> bool {
    let mut arrsac = Arrsac::new(Config::new(0.01), SmallRng::from_seed([0; 16]));

    // Use EPSILON_APPROX to force minimum distance.
    let p1_cam = [p1.0, p1.1, EPSILON_APPROX + p1.2.abs()];
    let p2_cam = [p2.0, p2.1, EPSILON_APPROX + p2.2.abs()];
    let p3_cam = [p3.0, p3.1, EPSILON_APPROX + p3.2.abs()];
    let p4_cam = [p4.0, p4.1, EPSILON_APPROX + p4.2.abs()];

    // 2d keypoints
    let p1_2d = Vector2::new(p1_cam[0], p1_cam[1]);
    let p2_2d = Vector2::new(p2_cam[0], p2_cam[1]);
    let p3_2d = Vector2::new(p3_cam[0], p3_cam[1]);
    let p4_2d = Vector2::new(p4_cam[0], p4_cam[1]);
    let p2ds = [p1_2d, p2_2d, p3_2d, p4_2d];

    // Stop if the keypoint's location on the frame is too close.
    for (a, b) in p2ds.iter().cartesian_product(&p2ds) {
        if (a - b).norm() < EPSILON_APPROX {
            return true;
        }
    }

    // Define the camera pose.
    let rotation = UnitQuaternion::from_euler_angles(rot.0, rot.1, rot.2);
    let translation = Translation::from(Vector3::new(trans.0, trans.1, trans.2));
    let pose = Isometry3::from_parts(translation, rotation);

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
    relative_eq!(rotation, pose.rotation, epsilon = EPSILON_APPROX)
        && relative_eq!(translation, pose.translation, epsilon = EPSILON_APPROX)
}
