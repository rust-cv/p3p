use crate::nordberg::{Pose, Sample};
use sample_consensus::{Estimator, Model};

impl Model for Pose {
    type Data = Sample;

    fn residual(&self, data: &Self::Data) -> f32 {
        self.error(data)
    }
}

/// This implements the [`sample_consensus::Estimator`] trait.
///
/// This uses the algorithm from the paper
/// "Lambda Twist: An Accurate Fast Robust Perspective Three Point (P3P) Solver"
/// to estimate four potential poses for the p3p problem.
pub struct NordbergEstimator;

impl Estimator for NordbergEstimator {
    type Model = Pose;
    type ModelIter = Vec<Pose>;
    const MIN_SAMPLES: usize = 3;
    fn estimate<'a, I>(&self, mut data: I) -> Self::ModelIter
    where
        I: Iterator<Item = &'a Sample> + Clone,
    {
        crate::nordberg::solve([
            data.next().unwrap().clone(),
            data.next().unwrap().clone(),
            data.next().unwrap().clone(),
        ])
    }
}
