# PostSPH.jl
Easing the process of post-processing simulations done in DualSPHysics. Currently supporting version 4.4 of DualSPHysics.

# Installation
# NOTE: Nothing is available yet, readme is just being developed
The code has been programmed for the stable v.1.1 of Julia and currently installation is done through:

```julia
using Pkg
Pkg.add("PostSPH")
```

# Usage
After installation the new function which is provided is:

```julia
PostSPH.readVtkArray(filename::Array{String,1},Cat::Enum)
```

* Example:
Go run the MovingSquare example in DualSPHysics (..\DualSPHysics_v4.4\examples\main\03_MovingSquare) and then from Julia navigate to the directory in which the .vtk files are stored. An example would be:

```julia
using PostSPH #Write this once in terminal or script file, to use function in PostSPH.jl
#Move to correct folder
cd(raw"..\DualSPHysics_v4.4\examples\main\03_MovingSquare\CaseMovingSquare_out\particles")
#Extracts all velocity data from all .vtk files and stores in array
#Second argument is velocity, clear a few lines further down
velSquare = readVtkArray("PartSquare",Cat(2)) 
```

To access data do, ```velSquare[1]``` to get vel data from "PartSquare_0000.vtk" and etc. Do ```velSquare[1][:,1]``` to get all x-components, 2 for y-components and 3 for z-components of velocity. 

The first argument, "filename", takes a string and searches for all similar strings. If only one file is wanted to be read in do;


```julia
velSquare = readVtkArray("PartSquare_0000.vtk",Cat(2))
```

Ie. provide full file name. All hard-coded properties are able to be extracted with the use of this function. The second argument "Cat", can therefore be:

Points,Idp,Vel,Rhop,Mass,Press,Vol,Ace,Vor,Typ and Mk, which each have assigned an integer, which can be found using the command ```Cat```:

```julia
Cat
Enum Cat:
Points = 0
Idp = 1
Vel = 2
Rhop = 3
Mass = 4
Press = 5
Vol = 6
Ace = 7
Vor = 8
Typ = 9
Mk = 10
```

Note it is possible to define ``` Idp = Cat(0) ``` etc. if it is needed to increase readability.

# Performance Tip

Due to the way Julia functions, it will recompile a function if input/output changes. A basic example would be running a function with a "Float32" and then a "Float64" input. If the input type is hardcoded Julia will under the hood recompile a type specific version for each input type. Therefore, as an user, if you experience slow initial performance, when trying to load a new variable, try in the first run only to read one file, ie. ```velSquare = readVtkArray("PartSquare_0000.vtk",Cat(2))``` and then afterwards read all files, by using only "PartSquare". 

# Current Implementation (Version 0.1.0)

Currently it is possible to:
1. Read single or all .vtk files existing in a directory, matching a string ie. "PartSquare".

# To do

Realistic ideas for the future:

1. Develop algorithm to give a number describing particle spacing ie. has the simulation coarsened over time or kept uniform
1. Implement a version of DualSPHysics own tools for post-processing in Julia ie. velocity measurement, surface elevation etc.
1. Further automatic post-processing calculations

# Credits

Version 0.0 has been developed by;

* Ahmed Salih (AhmedSalih3d) - Idea instigator, user of DualSPHysics for about a year, primary maintainer  of code now
* Saif Salih (Sayfsal) - Developer of initial custom .vtk reader
