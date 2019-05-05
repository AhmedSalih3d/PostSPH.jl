#File holding some older tests

function test2()
filenames   = readdir()
filterFiles = split.(filenames,"_")

k::Array{String,1} = []
for i = 1:length(filterFiles)
   k = push!(k,filterFiles[i][1])
end
deleteat!(k,length(k))

return unique(k)
end


function getfilestrings(namelist, pattern::Regex=".*")  # default is to accept all strings
    filenames = Vector{String}()  # this makes sure that the list is unique
    for filename in namelist
        m = match(pattern, filename)
        isnothing(m) && continue
        str = first(m.captures) #? captures for help, it just gets the "PartFluid", component
        !isempty(str) && push!(filenames, str)
    end
    return unique(filenames)
end

#list = getfilestrings(readdir(), r"([a-zA-Z]+).*\.vtk$")

function test3()
function glob_thing(pattern, to_search::AbstractVector)::AbstractVector
    filter(x -> occursin(pattern, x), to_search)
end

function name_beginning(name::AbstractString)
    names_no_extension = splitext(name)[1]
    first_part = split(names_no_extension, "_")[1]
    return first_part
end

files = readdir()
vtk_files = glob_thing(".vtk", files)
name_beginnings = name_beginning.(vtk_files)
return unique(name_beginnings)
end
