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

## License

This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.

This rewrite is based on the adaptation of the original code (GPL-3.0)
into [OpenMVG, published under MPL-2.0 with the original author agreement][p3p-openmvg].

[p3p-openmvg]: https://github.com/openMVG/openMVG/pull/1500
