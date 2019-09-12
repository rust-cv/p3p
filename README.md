# rust-photogrammetry/p3p

[![][badge-crate]][crate]
[![][badge-doc]][doc]
[![][badge-license]][license]

[badge-crate]: https://img.shields.io/badge/crate-latest-yellow.svg?style=flat-square
[badge-doc]: https://img.shields.io/badge/documentation-latest-yellow.svg?style=flat-square
[badge-license]: https://img.shields.io/badge/license-MPL--2.0-blue.svg?style=flat-square

[crate]: https://crates.io/crates/p3p
[doc]: https://docs.rs/p3p
[license]: https://www.mozilla.org/en-US/MPL/2.0/

Camera pose estimation given 3D points and corresponding pixel coordinates.

Implementation based on
"Lambda Twist: An Accurate Fast Robust Perspective Three Point (P3P) Solver"
Persson, M. and Nordberg, K. ECCV 2018.

TODO:
 - [ ] Improve readme
 - [ ] Improve documentation

## License

This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.

This rewrite is based on the adaptation of the original code (GPL-3.0)
into [OpenMVG, published under MPL-2.0 with the original author agreement][p3p-openmvg].

[p3p-openmvg]: https://github.com/openMVG/openMVG/pull/1500

## Contributions

All forms of contributions are welcomed, **preferably first as github issues**.

- Questions
- Documentation
- Tests
- Benchmarks
- Features

In case of contribution to source code,
it needs to use [rustfmt][rustfmt] and [clippy][clippy].
To run clippy:

```
touch src/lib.rs; cargo clippy -- -W clippy::all -W clippy::nursery -W clippy::pedantic
```

[rustfmt]: https://github.com/rust-lang/rustfmt
[clippy]: https://github.com/rust-lang/rust-clippy
