# NOTICE
I am actively developing it, some features especially related to vtk have been
removed for now. If you were using this package, please use an older version
from late 2019, to keep the same functions as before.

# PostSPH.jl
Easing the process of post-processing simulations done in DualSPHysics. Currently supporting version 5.0 of DualSPHysics. Please make a issue if you are not able read bi4's generated using version 5.

# Installation
The code has been tested on version 1.5 of Julia and currently installation is done through:

```julia
] add https://github.com/AhmedSalih3d/PostSPH.jl
```
Bracket, "]", is to enter "pkg" mode in Julia shell. The reason for install being this way currently is that it is an unofficial package right now. To use the package write, ```using PostSPH```.

If you want to test the "dev" version or any other branch, then instead write:

```julia
] add https://github.com/AhmedSalih3d/PostSPH.jl#dev
```

If you encounter an issue using a development branch, please post an issue.

# Usage
After installation some different functions are available, which are specifically related to the bi4 file type. Previously working with vtk's was also supported (2019 repo), but now focus is on bi4.

```julia
readBi4Array(typ::Cat)
```

* Example:
Please go to the "example" folder in this repo. Feel free to add your own examples or use cases if you believe it is relevant for others too.

# Subfunctions in this package

To make this module easier to maintain functionality has been split into subfunctions
which the main functions utilize. Every function preprended with a underscore, ' _ ' is not meant
to be used by the user, but is still usable. 

# Performance Tip

Due to the way Julia functions, it will recompile a function if input/output changes. A basic example would be running a function with a "Float32" and then a "Float64" input. If the input type is hardcoded Julia will under the hood recompile a type specific version for each input type. Therefore, as an user, if you experience slow initial performance, when trying to load a new variable, try in the first run only to read one file, ie. ```vel_array = readBi4Array("Part_0000.bi4",Cat(3))``` and then afterwards read all files, by using only "Part".

# Current Implementation (Version 0.3.0)

Currently it is possible to:
1. Read single or all .bi4 files in a directory, simply mentioning the array in question

# To do

Realistic ideas for the future:

1. Develop algorithm to give a number describing particle spacing ie. has the simulation coarsened over time or kept uniform
2. Implement a version of DualSPHysics own tools for post-processing in Julia ie. velocity measurement, surface elevation etc.
3. Further automatic post-processing calculations

# Credits

Version 0.3.0 has been developed by;

* Ahmed Salih (AhmedSalih3d) - Idea instigator, user of DualSPHysics for about a year, primary maintainer  of code now
* Saif Salih (sayfsal) - Developer of initial custom .vtk reader
