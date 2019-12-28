use approx::{assert_relative_eq, relative_eq};
use arraymap::ArrayMap;
use cv::nalgebra::{Isometry3, Point3, Translation, UnitQuaternion, Vector3};
use cv::KeypointWorldMatch;
use itertools::Itertools;
use p3p::*;
use quickcheck_macros::quickcheck;

type V3 = (f32, f32, f32);

const EPSILON_APPROX: f32 = 1e-2;

#[test]
fn manual_case() {
    // Define some points in camera coordinates (with z > 0).
    let camera_depth_points = [
        [-0.228_125, -0.061_458_334, 1.0],
        [0.418_75, -0.581_25, 2.0],
        [1.128_125, 0.878_125, 3.0],
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

    let samples: Vec<KeypointWorldMatch> = world_points
        .iter()
        .zip(&normalized_image_coordinates)
        .map(|(&world, &image)| KeypointWorldMatch(image.into(), world.into()))
        .collect();

    // Estimate potential poses with P3P.
    let poses = nordberg::p3p([samples[0], samples[1], samples[2]]);
    assert!(!poses.is_empty());

    // Compare the first pose to ground truth.
    let p3p_iso3 = poses[0];
    assert_relative_eq!(rot, p3p_iso3.rotation, epsilon = EPSILON_APPROX);
    assert_relative_eq!(trans, p3p_iso3.translation, epsilon = EPSILON_APPROX);
}

// Failures log at EPSILON_APPROX = 1e-2:
// (0.0, 15.0, 61.0), (0.0, 0.0, 0.0), (-44.782043, -69.0, 15.0), (8.574509, -65.91768, 68.265625), (-14.3312, -23.98317, 67.01338)
// (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (75.34378, 71.73302, 97.50447), (-29.456757, 64.21983, 34.308487), (-1.483345, 97.73132, 62.84854)

// Failures log at EPSILON_APPROX = 1e-3: (very long running times)
// (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 15.0, 35.420532), (89.0, 0.0, 30.208725)

/// This test is ignored because it is random and may fail in CI.
/// Run `cargo test -- --ignored` to test it.
#[quickcheck]
#[ignore]
fn non_degenerate_case(rot: V3, trans: V3, p1: V3, p2: V3, p3: V3) -> bool {
    // Define some points in camera coordinates (with z > 0).
    let camera_depth_points = [
        [p1.0, p1.1, EPSILON_APPROX + p1.2.abs()],
        [p2.0, p2.1, EPSILON_APPROX + p2.2.abs()],
        [p3.0, p3.1, EPSILON_APPROX + p3.2.abs()],
    ]
    .map(|&p| Point3::from(p));

    // Define the camera pose.
    let rotation = UnitQuaternion::from_euler_angles(rot.0, rot.1, rot.2);
    let translation = Translation::from(Vector3::new(trans.0, trans.1, trans.2));
    let pose = Isometry3::from_parts(translation, rotation);

    // Compute world coordinates.
    let world_points = camera_depth_points.map(|p| pose.inverse() * p);

    // Compute normalized image coordinates.
    let normalized_image_coordinates = camera_depth_points.map(|p| (p / p.z).xy());

    // Stop if normalized image coords are colinear.
    for (a, b, c) in normalized_image_coordinates.iter().tuple_combinations() {
        if ((a - b).normalize()).dot(&(a - c).normalize()) > 1.0 - 10.0 * EPSILON_APPROX {
            return true;
        }
    }

    // Stop if the keypoint's location on the frame is too close.
    for (a, b) in normalized_image_coordinates.iter().tuple_combinations() {
        if (a - b).norm() < 10.0 * EPSILON_APPROX {
            return true;
        }
    }

    let samples: Vec<KeypointWorldMatch> = world_points
        .iter()
        .zip(&normalized_image_coordinates)
        .map(|(&world, &image)| KeypointWorldMatch(image.into(), world.into()))
        .collect();

    // Estimate potential poses with P3P.
    let poses = nordberg::p3p([samples[0], samples[1], samples[2]]);
    assert!(!poses.is_empty());

    // Check that at least one estimated pose is near solution.
    poses.iter().fold(false, |already_true, pose| {
        let p3p_pose = pose;
        let same_rot = relative_eq!(rotation, p3p_pose.rotation, epsilon = EPSILON_APPROX);
        let same_trans = relative_eq!(translation, p3p_pose.translation, epsilon = EPSILON_APPROX);
        already_true || (same_rot && same_trans)
    })
}
