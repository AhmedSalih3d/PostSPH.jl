using EzXML
using XMLDict
using OrderedCollections

file = raw"D:\DualSPHysics_v5.0\examples\main\03_MovingSquare\CaseMovingSquare_Def.xml"
xmlfile = []
open(file,"r") do io
    xmlfile = readxml(io)
    xmlroot = root(xmlfile)
    a = []
    
    # Casedef and Execution
    for i in eachelement(xmlroot)
        # Constantsdef, mkconfig, parameters etc.
        for j in eachelement(i)

                

                dct = OrderedDict{String,OrderedDict{String,String}}()
                
                if haselement(j) == false
                    atr_name = nodename(j)
                    dct[atr_name] = OrderedDict{String,String}();
                    for l in eachattribute(j)
                        dct[atr_name][l.name] = l.content
                    end
                else
                    for k in eachelement(j)
                        for kk in eachelement(k)
                            loop!(dct,kk)
                        end
                    end
                    #loop!(dct,j)
                end

            push!(a,dct)
        end
    end

end

function loop!(dct,j)
    for k in eachelement(j)
        atr_name = nodename(k)
        dct[atr_name] = OrderedDict{String,String}();
        for l in eachattribute(k)
            dct[atr_name][l.name] = l.content
        end
    end
end

case = firstnode(xmlroot)

case_elements     = elements(xmlroot)

case_sub_elements = elements.(case_elements) 

has_elements_bool1  = map(x->haselement.(x),case_sub_elements)
has_elements_bool2  =map(x->hasnextelement.(x),case_sub_elements)


### ###


file = raw"D:\DualSPHysics_v5.0\examples\main\03_MovingSquare\CaseMovingSquare_Def.xml"
xmlfile = readxml(file)
d = xml_dict(xmlfile)

mainlist_fix = d["case"]["casedef"]["geometry"]["commands"]["mainlist"][""]
deleteat!(d["case"]["casedef"]["geometry"]["commands"]["mainlist"][""],1:2:length(mainlist_fix))


### ### ### 

doc = parsexml("""
<?xml version="1.0" encoding="UTF-8" ?>
<case>
</case>
""")

rd  = root(doc)

constantsdef = addelement!(rd,"constantsdef")

gravity      = addelement!(constantsdef,"gravity")

gx = AttributeNode("x","0")
gy = AttributeNode("y","0")
gz = AttributeNode("z","-9.81")

link!(gravity,gx)
link!(gravity,gy)
link!(gravity,gz)

comment=AttributeNode("comment","Gravitational acceleration") 
units_comment=AttributeNode("units_comment","m/s^2")

link!(gravity,comment)
link!(gravity,units_comment)

rhop0      = addelement!(constantsdef,"rhop0")

rhop0_value = AttributeNode("value","1000")


link!(rhop0,rhop0_value)

comment=AttributeNode("comment","Reference density of the fluid") 
units_comment=AttributeNode("units_comment","kg/m^3")

link!(rhop0,comment)
link!(rhop0,units_comment)

write("Test.xml",doc)


function adddictelement!(parentelement,dict)
    k = collect(keys(dict))
    for i in k
        p = addelement!(parentelement,i)
        for j in k
            tmp = AttributeNode(i,j)
            print(tmp)
            link!(p,tmp)
        end
    end
end