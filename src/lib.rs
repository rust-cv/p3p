#![warn(missing_docs)]

//! This package provides functions to solve camera pose estimation
//! given a set of three 3D points and their corresponding pixel coordinates.
//!
//! This problem is generally known as ["Perspective-n-Point" (PnP)][pnp].
//! We focus on the P3P case (n = 3) which is the
//! minimal amount of points required to solve the problem with a finite number of solutions.
//! There are many approaches to solve P3P. In this package, we have implemented:
//!
//!  - Lambda Twist: An Accurate Fast Robust Perspective Three Point (P3P) Solver.
//!    Mikael Persson, Klas Nordberg. ECCV 2018. ([paper][lambda-twist])
//!
//! [pnp]: https://en.wikipedia.org/wiki/Perspective-n-Point
//! [lambda-twist]: http://openaccess.thecvf.com/content_ECCV_2018/html/Mikael_Persson_Lambda_Twist_An_ECCV_2018_paper.html

pub mod nordberg;
