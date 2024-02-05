CADET-Julia
======

CADET-Julia is a numerical solver for chromatography and is a streamlined and efficient adaption of the more comprehensive [CADET](https://github.com/modsim/CADET) project, and is based on an implementation of the Discontinuous Galerkin Spectral Element Method (DGSEM), similar to the one in CADET.

Overview
-------

[CADET](https://github.com/modsim/CADET) is the foundational C++ project, with a comprehensive codebase and extensive functionality.
Recognizing the need for a more accessible and agile solution, this Julia implementation was created.
It retains the powerful DGSEM but focuses on providing a simplified and readable alternative with reduced functionality.
Our main targets with this project include
- **Rapid prototyping:** Whether you are familiar with the original C++ project or exploring CADET for the first time, this Julia implementation is designed for swift experimentation and prototyping, enabling a more agile development workflow.
- **Benchmarking:** We leverage the reduced codebase to benchmark the performance of methods implemented in both C++ and Julia. The streamlined Julia implementation not only enhances efficiency but also allows for the exploration of improved performance, e.g., through alternative time integration methods. Additionally,  the two code bases can be used for mutual verification, reinforcing our confidence in the correct solution of the problem.

For more details about the original CADET, we refer to the [webpage](https://cadet.github.io/master/index.html#).
 

Getting started
-------

Clone this repository: `git clone https://github.com/jespfra/CADET-Julia.git`

To set up and run simulations, please refer to the example folder.
 
License
-------

Released under GPL v3. License (see `LICENSE.txt <https://github.com/jespfra/CADET-Julia/blob/main/LICENSE.txt>`_)::

Acknowledgments
---------------

### Contributors
* [Jesper Frandsen](https://github.com/jespfra) (Technical University of Denmark, Copenhagen, Denmark)
* [Jan Breuer](https://github.com/jbreue16) (Forschungszentrum Juelich GmbH, IBG-1: Biotechnology, Juelich, Germany)
* [Johannes Schmölder](https://github.com/schmoelder) (Forschungszentrum Juelich GmbH, IBG-1: Biotechnology, Juelich, Germany)

