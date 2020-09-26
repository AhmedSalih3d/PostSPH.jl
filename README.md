# NOTICE
I am actively developing it, some features especially related to vtk have been
removed for now. If you were using this package, please use an older version
from late 2019, to keep the same functions as before.

# PostSPH.jl
Easing the process of post-processing simulations done in DualSPHysics. Currently supporting version 4.4 of DualSPHysics.

# Installation
The code has been programmed for the stable v.1.1.1 of Julia and currently installation is done through:

```julia
] add https://github.com/AhmedSalih3d/PostSPH.jl
```
Bracket, "]", is to enter "pkg" mode in Julia shell. The reason for install being this way currently is that it is an unofficial package right now. To use the package write, ```using PostSPH```.

If you want to test the "dev" version, then instead write:

```julia
] add https://github.com/AhmedSalih3d/PostSPH.jl#dev
```

This is highly unrecommended but possible.

# Usage
After installation some different functions are available:

```julia
readVtkArray(filename::Array{String,1},Cat::Enum)
readVtkParticles(filename::Array{String,1})
readBi4Array(typ::Cat)
```

readVtkArray is a function which loads an array in from all simulation files,
while readVtkParticles allows an user to get the number of particles in each .vtk
file, which can be quite useful. readBi4Array is a function which lets you extract
all particles from a specific array type.

* Example:
Go run the MovingSquare example in DualSPHysics (..\DualSPHysics_v4.4\examples\main\03_MovingSquare) and then from Julia navigate to the directory in which the .vtk files are stored. An example would be:

```julia
using PostSPH #Write this once in terminal or script file, to use function in PostSPH.jl
#Move to correct folder
cd(raw"..\DualSPHysics_v4.4\examples\main\03_MovingSquare\CaseMovingSquare_out\particles")
#Extracts all velocity data from all .vtk files and stores in array
#Second argument is velocity, clear a few lines further down
velSquare = readVtkArray("PartSquare",Cat(2))
#To see number of particles in each .vtk file do then:
nSquare   = readVtkParticles("PartSquare")
```

Note for non-fluid .vtk files the number of particles should be constant.

To access data do, ```velSquare[1]``` to get vel data from "PartSquare_0000.vtk" and etc. Do ```velSquare[1][:,1]``` to get all x-components, 2 for y-components and 3 for z-components of velocity.

The first argument, "filename", takes a string and searches for all similar strings. If only one file is wanted to be read in do;


```julia
#Different ways to write the same thing
velSquare = readVtkArray("PartSquare_0000.vtk",Cat(2))
velSquare = readVtkArray("PartSquare_0000.vtk",PostSPH.Vel)
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

To use the readBi4Array functionality simply do:

```julia
rhop_array = readBi4Array(PostSPH.Rhop)
```
So now the density of all particles in each time step has been saved into the variable "rhop_array".

# Subfunctions in this package

To make this module easier to maintain functionality has been split into subfunctions
which the main functions utilize. Every function preprended with a underscore, ' _ ' is not meant
to be used by the user.

# Performance Tip

Due to the way Julia functions, it will recompile a function if input/output changes. A basic example would be running a function with a "Float32" and then a "Float64" input. If the input type is hardcoded Julia will under the hood recompile a type specific version for each input type. Therefore, as an user, if you experience slow initial performance, when trying to load a new variable, try in the first run only to read one file, ie. ```velSquare = readVtkArray("PartSquare_0000.vtk",Cat(2))``` and then afterwards read all files, by using only "PartSquare". This also applies to readVtkParticles.

# Current Implementation (Version 0.2.0)

Currently it is possible to:
1. Read single or all .vtk files existing in a directory, matching a string ie. "PartSquare", using "readVtkArray".
2. Extract number of particles in each simulation step, matching a string ie. "PartSquare", using "readVtkParticles"
3. Extract information about total mass in each vtk file at each time step, matching a string ie. "PartFluid", using "MassVtk"
4. Calculate forces (in Newton) exerted on the all particles included in a vtk file, matching a string, ie. "PartSquare", using "ForceVtk"
5. Read single or all .bi4 files in a directory, simply mentioning the array in question

Some further tools have been developed to ease useability:

1. "readVtkVariables" has been developed which allows an user to check the available parameters in
    a single vtk file.
2. "readVtkNames" has been developed so the user can get an array holdning the existing vtk files.

# To do

Realistic ideas for the future:

1. Develop algorithm to give a number describing particle spacing ie. has the simulation coarsened over time or kept uniform
2. Implement a version of DualSPHysics own tools for post-processing in Julia ie. velocity measurement, surface elevation etc.
3. Further automatic post-processing calculations

# Credits

Version 0.2.0 has been developed by;

* Ahmed Salih (AhmedSalih3d) - Idea instigator, user of DualSPHysics for about a year, primary maintainer  of code now
* Saif Salih (sayfsal) - Developer of initial custom .vtk reader
