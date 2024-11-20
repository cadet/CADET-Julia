[![DOI](https://zenodo.org/badge/690988707.svg)](https://doi.org/10.5281/zenodo.14068012)
![CI](https://github.com/cadet/cadet-julia/actions/workflows/CI.yml/badge.svg)

CADET-Julia
===========

CADET-Julia is a numerical solver for chromatography and is a streamlined and efficient adaption of the more comprehensive [CADET-Core](https://github.com/modsim/CADET) project, and is based on a spatial Discontinuous Galerkin Spectral Element Method (DGSEM), similar to the one implemented in CADET-Core.
If you find it useful for your own work, we would appreciate acknowledgements of this software and citations of our papers:

- Jesper Frandsen, Jan M. Breuer, Johannes Schmölder, Jakob K. Huusom, Krist V. Gernaey, Jens Abildskov, Eric von Lieres: [CADET-Julia: Efficient and versatile, open-source simulator for batch chromatography in Julia](https://doi.org/10.1016/j.compchemeng.2024.108913), Computers & Chemical Engineering, Volume 192 (2025), 108913

Overview
--------

[CADET-Core](https://github.com/modsim/CADET) is the foundational C++ project, with a comprehensive codebase and extensive functionality.
Recognizing the need for a more accessible and agile solution, this Julia implementation was created.
It retains the computationally powerful DGSEM but focuses on providing a simplified and readable alternative with reduced functionality.
Our main targets with this project include
- **Rapid prototyping:** Whether you are familiar with the original C++ project or exploring CADET-Core for the first time, this Julia implementation is designed for swift experimentation and prototyping, enabling a more agile development workflow.
- **Benchmarking:** We leverage the reduced codebase to benchmark the performance of methods implemented in both C++ and Julia. The streamlined Julia implementation not only enhances efficiency but also allows for the exploration of improved performance, e.g., through alternative time integration methods. Additionally,  the two code bases can be used for mutual verification, reinforcing our confidence in the correct solution of the problem.

For more details about the CADET-project, including the original CADET-Core, we refer to the [webpage](https://cadet.github.io/master/index.html#).
 

Getting started
---------------

Clone this repository: `git clone https://github.com/jespfra/CADET-Julia.git`

Install Julia, preferably as an extension to your IDE such as visual studio code.

Open a Julia REPL at the project directory and run the following commands:

```julia
using Pkg
Pkg.instantiate()
```

You can compile and execute any of the simulations defined in the test folder.

To set up and run your own simulations, please refer to the test folder.
 
License
-------

Released under GPL v3. License (see [here](https://github.com/jespfra/CADET-Julia/blob/main/LICENSE))
