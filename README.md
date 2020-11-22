# article-powmde

This repository includes source for experiments in the preprint [Computing the matrix fractional power based on the double exponential formula, arXiv: hoge.fuga](https://example.com/).

## How to run the code?
Our code is tested on the following environment:

- OS: Linux (Ubuntu 20.04)
- CPU: Intel(R) Core(TM) i7-9700K
- Memory: 16GB
- Julia version: 1.5.1

The source code is written in [Julia Language](https://julialang.org/).
Hence, Julia Language has to be installed in your computer.

In addition, some additional packages of Julia are required.

- [Arpack.jl](https://github.com/JuliaLinearAlgebra/Arpack.jl)
- [Elliptic.jl](https://github.com/nolta/Elliptic.jl)
- [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl)
- [MatrixMarket.jl](https://github.com/JuliaSparse/MatrixMarket.jl)
- [SuiteSparse.jl](https://github.com/JuliaLinearAlgebra/SuiteSparse.jl)

To install these packages, run `julia install_packages.jl` in your terminal.

The test matrices have to be in `matrix` directory in `.mtx` format, for example, `matrix/ex5.mtx`.
If you want to download, extract and put these matrices automatticaly, run `julia download_matrices.jl`.

In order to generate test matrices, `SPD_(well|ill)`, `NS_(well|ill)`, `poisson200`, run `generate_matrices.jl`.

In order to run all tests, run `julia run_alltests.jl` in your terminal.
