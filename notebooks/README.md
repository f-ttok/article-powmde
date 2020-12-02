In order to plot figures, we use both of Python (with Matplotlib) and Julia (with PyPlot).

1. Install Python
1. `pip install numpy pandas matplotlib jupyter notebook`
1. Install Julia packages
    - `julia`
    - `ENV["PYTHON"] = python; ENV["JUPYTER"] = jupyter`
    - `(@1.5) pkg> add IJulia PyPlot` in the Pkg mode (You can enter the Pkg mode by pressing `]` in the REPL. See [the official document](https://docs.julialang.org/en/v1/stdlib/Pkg/).).
    - After installing IJulia and PyPlot, press backspace to get back to the REPL.
    - Ctrl + C to stop Julia Repl
1. Run `jupyter notebook` in `notebooks` directory

