# article-powmde

This repository includes programs for experiments in the preprint [Computing the matrix fractional power based on the double exponential formula, arXiv: hoge.fuga](https://example.com/).

## How to run the code?

1. Install Julia and set path.
1. `git clone https://github.com/f-ttok/article-powmde.git` (or download zip file from https://github.com/f-ttok/article-powmde/archive/master.zip, and extract the zip file.)
1. `cd article-powmde`
1. `julia install_packages.jl`
1. `julia download_matrices.jl`
1. `julia generate_matrices.jl`
1. `julia run_alltests.jl`

### Details
Our environment is as follows:

- OS: Linux (Ubuntu 20.04)
- CPU: Intel(R) Core(TM) i7-9700K
- Memory: 16GB
- Julia version: 1.5.1

The programs are written in [Julia Language](https://julialang.org/).
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

In order to generate test matrices, `SPD_(well|ill)`, `NS_(well|ill)`, `poisson200`, run `julia generate_matrices.jl`.

In order to run all tests, run `julia run_alltests.jl` in your terminal.


## Some Comments
### Test 5
Julia compile the code before running the code for the first time.
Hence, it may affect the computational time.
In order to avoid the effect, we measure the computational time twice for `poisson200` and `cell1`.

### Other comments
Figures are plotted via Jupyter notebook.
Their results are in `notebooks` directory.