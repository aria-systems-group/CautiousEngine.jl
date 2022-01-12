# CautiousEngine.jl
This package contains end-to-end tools for the verification and synthesis of data driven systems using Gaussian process regression and IMDP abstraction. 

## Installation
Installation has only been tested for Julia 1.5.x. Instantiation fails for Julia 1.6+, likely due to old dependencies. 

Since this package is not public yet, simply clone and add the repository to the Julia load path (or add it to `startup.jl`) via:

`push!(LOAD_PATH, "path/to/CautiousEngine.jl")`. 

If `~/.julia/config/startup.jl` does not exist, you can create it. A good location for packages under development is `~/.julia/dev`.

Next, use the following commands to activate the project environment and then instantiate the project and install dependencies: 

`pkg"activate path/to/CautiousEngine.jl"` then `pkg"instantiate"`.

This will install all the required packages.

### BMDP Tool
This package depends on the `bmdp-tool` here: https://github.com/aria-systems-group/bmdp-tool

The tool should be compiled using Make and the `synthesis` executable moved to a location on the user executable path e.g. `/usr/local/bin`. CautiousEngine does not use the included MATLAB tools.

## Usage
Full documentation is in progress, but here are some notes to get started.

Full examples are provided in `examples/` for synthesis and verification applications. The synthesis example initially uses a small number of points to test that the package is working. Copy this file and modify it as needed. 

**Note how measurement and process noise are handled in the scripts.** The synthesis procedure can only use process noise, and the verification procedure can only use measurement noise.

**The synthesis example uses an optional known component of the dynamics.** This is currently limited to knowing `x_k` in `x_{k+1} = x_k + g(x_k)` but not `g(x_k)`. 

## Contributing
Contributions are welcome, especially those related to testing, code optimization, documentation, etc. 
