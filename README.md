# rust-photogrammetry/p3p

Camera pose estimation given 3D points and corresponding pixel coordinates.

Implementation based on
"Lambda Twist: An Accurate Fast Robust Perspective Three Point (P3P) Solver"
Persson, M. and Nordberg, K. ECCV 2018.

TODO:
 - [ ] Improve readme
 - [ ] Improve documentation


Run clippy with:
```
touch src/lib.rs; cargo clippy -- -W clippy::all -W clippy::nursery -W clippy::pedantic
```
