# PostSPH.jl
Easing the process of post-processing simulations done in DualSPHysics

# Installation
NOTE: Nothing is available yet, readme is just being developed
The code has been programmed for the stable v.1.1 of Julia and currently installation is done through:

```julia
using Pkg
Pkg.add("PostSPH")
```

# Usage
After installation the new functions which are provided are:

```julia
readVtkArray(filenames::Array{String,1},cat::Enum)
```

* Example:
Go in and test MovingSquare in DualSPHysics main and then from Julia navigate to the directory in which the .vtk files are stored and use the command:

```julia
Just have to do it, takes some time
```

# Current Implementation (Version 0.0)

Currently it is possible to:
1. Read single or all .vtk files existing in a directory, matching a string ie. "PartFloating".
2. Automatically calculate mean force in all directions and magnitude of a fluid / solid

# To do

Realistic ideas for the future:

1. Develop algorithm to give a number describing particle spacing ie. has the simulation coarsened over time or kept uniform
1. Implement a version of DualSPHysics own tools for post-processing in Julia ie. velocity measurement, surface elevation etc.

