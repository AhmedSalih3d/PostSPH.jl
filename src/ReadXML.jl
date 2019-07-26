using LightXML

struct Fixed
    mkbound
    mk
    beg
    count
    property
end

struct Moving
    mkbound
    mk
    beg
    count
    property
    refmotion
end

struct Floating
    mkbound
    mk
    beg
    count
    property
end

struct Fluid
    mkfluid
    mk
    beg
    count
end

const BodyType  =  Dict("fixed" => Fixed,"moving" => Moving, "floating" => Floating, "fluid" => Fluid)

## Accesses the relevant part of the XML file array depending on BodyType ie. "k"
# is the array which "RelevantXmlSnip" produced. "BodySelect" is then used to
# grab the right structure.
function Access(k::Array{XMLElement,1},Body::String)
    n = length(k)
    j = Array{BodyType[Body],1}(undef,n)
    for i = 1:n
        j[i] = BodySelect(k[i],Body)
    end
    return j
end

##To ease conversion... works for now.
function p(m)
    parse(Int64,m)
end

## Depending on choice of body seperate values are stored
function BodySelect(k,Body)
    mkbound   = attribute(k,"mkbound")
    mkfluid   = attribute(k,"mkfluid")
    mk        = attribute(k,"mk")
    beg       = attribute(k,"begin")
    count     = attribute(k,"count")
    property  = attribute(k,"property")
    refmotion = attribute(k,"refmotion")
    if Body == "fixed"
        Fixed(p(mkbound),p(mk),p(beg),p(count),property)
    elseif Body == "moving"
        Moving(p(mkbound),p(mk),p(beg),p(count),property,p(refmotion))
    elseif Body == "floating"
        Floating(p(mkbound),p(mk),p(beg),p(count),property)
    elseif Body == "fluid"
        Fluid(p(mkfluid),p(mk),p(beg),p(count))
    end
end

## This function returns the relevant xml snippet for each BodyType ie. "fluid",
## "floating" etc. Also returns the parsed xml file, so that it can be closed later.
function RelevantXmlSnip(filename::String,Bodies::Array{String,1}=["fixed","moving","floating","fluid"])
    nBody = length(Bodies)
    j = Array{Array{XMLElement,1}}(undef,nBody)
    xdoc = parse_file(filename)
    xroot = root(xdoc)
    ces = xroot["execution"]
    e1  = ces[1]
    t   = find_element(e1, "particles")

    for i = 1:nBody
        j[i] = get_elements_by_tagname(t, Bodies[i])
    end

    return j,xdoc
end


## Generates final array based on former functions
function MkArray(filename::String,Bodies::Array{String,1}=["fixed","moving","floating","fluid"])
    y,xdoc = RelevantXmlSnip(filename)

    j = Array{Any,1}(undef,length(Bodies))

    for i = 1:length(j)
        j[i] = Access(y[i],Bodies[i])
    end

    free(xdoc)

    return j
end
