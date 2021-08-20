using EzXML

function constructElement(name,attr_dct::OrderedDict=OrderedDict{String,String}(); content::String="")
    p = ElementNode(name)
    p.content = content

    for (key,value) in attr_dct
        a = AttributeNode(key,value)
        link!(p,a)
    end



    return p
end

function getConstantsDef()
    doc = parsexml("""
    <constantsdef>            
        <gravity x="0" y="0" z="0" comment="Gravitational acceleration" units_comment="m/s^2" />
        <rhop0 value="1000" comment="Reference density of the fluid" units_comment="kg/m^3" />
        <hswl value="0" auto="true" comment="Maximum still water level to calculate speedofsound using coefsound" units_comment="metres (m)" />
        <gamma value="7" comment="Polytropic constant for water used in the state equation" />
        <speedsystem value="0" auto="true" comment="Maximum system speed (by default the dam-break propagation is used)" />
        <coefsound value="20" comment="Coefficient to multiply speedsystem" />
        <speedsound value="0" auto="true" comment="Speed of sound to use in the simulation (by default speedofsound=coefsound*speedsystem)" />
        <b value="1.1200e+05" auto="false" />
        <coefh value="1.0" comment="Coefficient to calculate the smoothing length (h=coefh*sqrt(3*dp^2) in 3D)" />
        <cflnumber value="0.2" comment="Coefficient to multiply dt" />
    </constantsdef>
    """)


end