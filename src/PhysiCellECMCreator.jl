module PhysiCellECMCreator

using Optimization, OptimizationOptimJL, Polynomials, LightXML, CSV, DataFrames

export generateICECM, createICECMXMLTemplate

"""
    generateICECM(path_to_ic_ecm_xml::String, path_to_ic_ecm_csv::String, config_dict::Dict{String,Float64})

Generate an initial condition ECM file from an XML file.

# Arguments
- `path_to_ic_ecm_xml::String`: the path to the XML file that describes the ECM
- `path_to_ic_ecm_csv::String`: the path to the CSV file that will be created
- `config_dict::Dict{String,Float64}`: a dictionary with the following keys:
    - `x_min::Float64`: the minimum x-coordinate of the PhysiCell domain
    - `x_max::Float64`: the maximum x-coordinate of the PhysiCell domain
    - `dx::Float64`: the spacing in the x-direction of voxels in the PhysiCell domain
    - `y_min::Float64`: the minimum y-coordinate of the PhysiCell domain
    - `y_max::Float64`: the maximum y-coordinate of the PhysiCell domain
    - `dy::Float64`: the spacing in the y-direction of voxels in the PhysiCell domain
    - `z0::Float64`: optional, the z-coordinate of the PhysiCell domain (default: 0.0)
"""
function generateICECM(path_to_ic_ecm_xml::String, path_to_ic_ecm_csv::String, config_dict::Dict{String,Float64})
    xml_doc = parse_file(path_to_ic_ecm_xml)
    ic_ecm = root(xml_doc)
    df = parseECM(ic_ecm, config_dict)
    free(xml_doc)
    CSV.write(path_to_ic_ecm_csv, df)
end

function parseECM(ic_ecm::XMLElement, config_dict::Dict{String,Float64})
    df = initializeDataFrame(config_dict)
    for layer in child_elements(ic_ecm)
        parseLayer!(df, layer, config_dict)
    end
    @assert all(.!df.is_missing) "All voxels must have their ECM properties defined."
    select!(df, Not(:is_missing))
    return df
end

function initializeDataFrame(config_dict::Dict{String,Float64})
    coords = []
    for d in ["x", "y"]
        c0 = config_dict["$(d)_min"]
        c1 = config_dict["$(d)_max"]
        dc = config_dict["d$(d)"]
        nc = 1e-16 + (c1 - c0) / dc |> ceil |> Int
        cc = c0 .+ (0.5:1:nc) .* dc
        push!(coords, cc)
    end
    x = repeat(coords[1]; outer=length(coords[2]))
    y = repeat(coords[2]; inner=length(coords[1]))
    z = "z" in keys(config_dict) ? config_dict["z0"] : 0.0
    n = length(x)
    return DataFrame(x=x, y=y, z=z, ecm_density=zeros(Float64, n), ecm_orientation_x=zeros(Float64, n), ecm_orientation_y=zeros(Float64, n), is_missing=trues(n))
end

function parseLayer!(df::DataFrame, layer::XMLElement, config_dict::Dict{String,Float64})
    layer_df = initializeDataFrame(config_dict)
    for patch in child_elements(layer)
        parsePatch!(layer_df, patch, config_dict)
    end
    updateDataFrame!(df, layer_df; allow_overwrite=true)
end

function parsePatch!(layer_df::DataFrame, patch::XMLElement, config_dict::Dict{String,Float64})
    patch_df = patchECM(patch, config_dict)
    updateDataFrame!(layer_df, patch_df; allow_overwrite=false)
end

function patchECM(patch::XMLElement, config_dict::Dict{String,Float64})
    patch_type = attribute(patch, "type")
    if patch_type == "everywhere"
        return parseEverywherePatch(patch, config_dict)
    elseif patch_type == "ellipse"
        return parseEllipsePatch(patch, config_dict)
    elseif patch_type == "elliptical_disc"
        return parseEllipticalDiscPatch(patch, config_dict)
    else
        throw(ArgumentError("Patch type $patch_type not recognized."))
    end
end

function parseEverywherePatch(patch::XMLElement, config_dict::Dict{String,Float64})
    df = initializeDataFrame(config_dict)
    density, orientation, anisotropy = parseECMFeatures(patch)
    df.ecm_density .= density
    if orientation == "random"
        theta = 2 * pi * rand(length(df.x))
        df.ecm_orientation_x .= anisotropy * cos.(theta)
        df.ecm_orientation_y .= anisotropy * sin.(theta)
    else
        throw(ArgumentError("Orientation $orientation not recognized for an everywhere ECM patch.\nRecognized orientations are: `random`"))
    end
    df.is_missing .= false
    return df
end

function parseEllipsePatch(patch, config_dict)
    df = initializeDataFrame(config_dict)
    e = parseEllipseParameters(patch)
        
    in_ellipse = [positionRelativeToEllipse(e, x, y)==:thickness for (x, y) in zip(df.x, df.y)]
    n = sum(in_ellipse)
    density, orientation, anisotropy = parseECMFeatures(patch)
    df.ecm_density[in_ellipse] .= density
    if orientation == "random"
        theta = 2 * pi * rand(n)
        df.ecm_orientation_x[in_ellipse] .= anisotropy * cos.(theta)
        df.ecm_orientation_y[in_ellipse] .= anisotropy * sin.(theta)
    elseif orientation == "parallel"
        all_orientations = [anisotropy * parallelFiberOrientation(e, x, y) for (x, y) in zip(df.x[in_ellipse], df.y[in_ellipse])]
        df.ecm_orientation_x[in_ellipse] .= [v[1] for v in all_orientations]
        df.ecm_orientation_y[in_ellipse] .= [v[2] for v in all_orientations]
    elseif orientation == "perpendicular"
        all_orientations = [anisotropy * perpendicularFiberOrientation(e, x, y) for (x, y) in zip(df.x[in_ellipse], df.y[in_ellipse])]
        df.ecm_orientation_x[in_ellipse] .= [v[1] for v in all_orientations]
        df.ecm_orientation_y[in_ellipse] .= [v[2] for v in all_orientations]
    else
        throw(ArgumentError("Orientation $orientation not recognized for an ellipse ECM patch.\nRecognized orientations are: `random`, `parallel`, `perpendicular`"))
    end
    df.is_missing[in_ellipse] .= false
    return df
end

function parseEllipticalDiscPatch(patch, config_dict)
    df = initializeDataFrame(config_dict)
    e = parseEllipseParameters(patch)

    in_disc = [insideEllipse(e, x, y) for (x, y) in zip(df.x, df.y)]
    n = sum(in_disc)
    density, orientation, anisotropy = parseECMFeatures(patch)
    df.ecm_density[in_disc] .= density
    if orientation == "random"
        theta = 2 * pi * rand(n)
        df.ecm_orientation_x[in_disc] .= anisotropy * cos.(theta)
        df.ecm_orientation_y[in_disc] .= anisotropy * sin.(theta)
    else
        throw(ArgumentError("Orientation $orientation not recognized for an elliptical disc ECM patch.\nRecognized orientations are: `random`"))
    end
    df.is_missing[in_disc] .= false
    return df
end

function parseEllipseParameters(patch::XMLElement)
    x0 = parse(Float64, find_element(patch, "x0") |> content)
    y0 = parse(Float64, find_element(patch, "y0") |> content)
    a = parse(Float64, find_element(patch, "a") |> content)
    b = parse(Float64, find_element(patch, "b") |> content)
    θ = parseRotation(patch)
    thickness = parseThickness(patch)
    return EllipticalECM((x0, y0), (a, b), θ, thickness)
end

function parseRotation(patch::XMLElement)
    rotation_element = find_element(patch, "rotation")
    if isnothing(rotation_element)
        return 0.0
    end
    θ = content(rotation_element)
    units = attribute(rotation_element, "units")
    uses_pi = contains(θ, "pi") || contains(θ, "π")
    if isnothing(units)
        is_radians = uses_pi
    else
        is_radians = startswith(units, "rad")
    end
    @assert !uses_pi || is_radians """
    If you use pi or π, you cannot specify that the units are degrees. Solutions:
    - Omit the units attribute: <rotation>π/6</rotation> # ECMCreator will interpret this as radians
    - Specify the units as radians: <rotation units="rad">π/6</rotation> OR <rotation units="radians">π/6</rotation> # anything that starts with "rad"
    - Convert to degrees: <rotation units="deg">30</rotation> OR <rotation>30</rotation> # ECMCreator will interpret this as degrees since it does not contain pi or π
    """
    if !is_radians
        return parse(Float64, θ) * pi / 180
    end
    if uses_pi
        try
            return eval(Meta.parse(θ))
        catch e
            println("""
            Error: the rotation ($θ) had pi or π but could not be parsed as a number.
            Try, for example,
                <rotation>π/6</rotation>
            Also, you can specify the units as degrees:
                <rotation units="deg">30</rotation>
            """
            )
            rethrow(e)
        end
    end
    return parse(Float64, θ)
end

function parseThickness(patch::XMLElement)
    thickness_element = find_element(patch, "thickness")
    if isnothing(thickness_element)
        return 0.0
    end
    thickness = content(thickness_element)
    return parse(Float64, thickness)
end

function parseECMFeatures(patch::XMLElement)
    density = parse(Float64, find_element(patch, "density") |> content)
    orientation = find_element(patch, "orientation") |> content
    anisotropy = parse(Float64, find_element(patch, "anisotropy") |> content)
    return density, orientation, anisotropy
end

function updateDataFrame!(df::DataFrame, new_df::DataFrame; allow_overwrite::Bool=false)
    if !allow_overwrite
        @assert all(df.is_missing .| new_df.is_missing) "Overlapping patches are not allowed within the same layer.\nRemove the overlap or change the patch layers to indicate which one is on top (a higher layer)."
    end
    new_inds = .!new_df.is_missing
    df[new_inds, :] .= new_df[new_inds, :]
end

struct EllipticalECM
    center::Tuple{Real,Real}
    axes::Tuple{Real,Real}
    rotation::Real # radians
    thickness::Real

    M_from_elliptical_coords::Matrix{Float64}
    M_to_elliptical_coords::Matrix{Float64}
end

function EllipticalECM(center::Tuple{Real,Real}, axes::Tuple{Real,Real}, rotation::Real, thickness::Real)
    M_from_elliptical_coords = [cos(rotation) -sin(rotation); sin(rotation) cos(rotation)]
    M_to_elliptical_coords = [cos(rotation) sin(rotation); -sin(rotation) cos(rotation)]
    return EllipticalECM(center, axes, rotation, thickness, M_from_elliptical_coords, M_to_elliptical_coords)
end

function toEllipticalCoordinates(e::EllipticalECM, x::Real, y::Real)
    return e.M_to_elliptical_coords * ([x, y] .- e.center)
end

function fromEllipticalCoordinates(e::EllipticalECM, x::Real, y::Real)
    return e.M_from_elliptical_coords * [x, y] .+ e.center
end

function _ellipseFunctional(x::Real, y::Real, a::Real, b::Real)
    return (x / a)^2 + (y / b)^2
end

function insideEllipse(e::EllipticalECM, x::Real, y::Real)
    x, y = toEllipticalCoordinates(e, x, y)
    return _ellipseFunctional(x, y, e.axes[1], e.axes[2]) <= 1
end

function insideEllipseThickness(e::EllipticalECM, x::Real, y::Real)
    x, y = toEllipticalCoordinates(e, x, y)
    in_thickness = _ellipseFunctional(x, y, e.axes[1] + e.thickness, e.axes[2] + e.thickness) <= 1 &&
                   _ellipseFunctional(x, y, e.axes[1], e.axes[2]) >= 1
    return in_thickness
end

function outsideEllipseThickness(e::EllipticalECM, x::Real, y::Real)
    x, y = toEllipticalCoordinates(e, x, y)
    return _ellipseFunctional(x, y, e.axes[1] + e.thickness, e.axes[2] + e.thickness) > 1
end

function positionRelativeToEllipse(e::EllipticalECM, x::Real, y::Real)
    x, y = toEllipticalCoordinates(e, x, y)
    if _ellipseFunctional(x, y, e.axes[1], e.axes[2]) <= 1
        return :inside
    elseif _ellipseFunctional(x, y, e.axes[1] + e.thickness, e.axes[2] + e.thickness) <= 1
        return :thickness
    else
        return :outside
    end
end

function perpendicularFiberOrientation(e::EllipticalECM, x::Real, y::Real)
    x, y = toEllipticalCoordinates(e, x, y)
    a, b = e.axes
    objfn(theta, x, y, a, b) = (x - a * cos(theta))^2 + (y - b * sin(theta))^2
    objfn(theta, p) = objfn(theta[1], p...)
    optf = OptimizationFunction(objfn, Optimization.AutoForwardDiff())
    prob = OptimizationProblem(optf, [atan(y, x)], [x; y; a; b], lb=[-pi], ub=[pi])
    sol = solve(prob, BFGS())
    optimal_theta = sol.u[1]
    vec = [x, y] - [a * cos(optimal_theta), b * sin(optimal_theta)]
    vec ./= sqrt(vec[1]^2 + vec[2]^2)
    return e.M_from_elliptical_coords * vec
end

function parallelFiberOrientation(e::EllipticalECM, x::Real, y::Real; method::Symbol=:rotate_perpendicular)
    if method==:solve_polynomial
        @warn("ECMCreator: the :solve_polynomial method seems to have a bug. Better to use :rotate_perpendicular at this time.")
        x, y = toEllipticalCoordinates(e, x, y)
        a, b = e.axes

        P1 = Polynomial([a, 1])^2
        P2 = Polynomial([b, 1])^2
        P = x^2 * P2 + y^2 * P1 - P1 * P2
        r = roots(P)
        r = r[imag.(r).==0]
        r = Real.(r)
        r = r[r.>=0]

        ratio = (a + r[1]) / (b + r[1])
        return e.M_from_elliptical_coords * [-ratio * y, ratio * x]
    elseif method==:rotate_perpendicular
        return [0;1;;-1;0] * perpendicularFiberOrientation(e, x, y) # rotate 90 degrees counterclockwise
    else
        throw(ArgumentError("Method $method not recognized."))
    end
end

"""
    createICECMXMLTemplate

Create folder and create a template XML file for IC ECM.

Create the whole path to `path_to_folder` and then create `ecm.xml` in that folder.
"""
function createICECMXMLTemplate(path_to_folder::String)
    path_to_ic_ecm_xml = joinpath(path_to_folder, "ecm.xml")
    mkpath(dirname(path_to_ic_ecm_xml))
    xml_doc = XMLDocument()
    xml_root = create_root(xml_doc, "ic_ecm")

    e_layer = new_child(xml_root, "layer")
    set_attribute(e_layer, "ID", "1")

    e_patch = new_child(e_layer, "patch")
    set_attributes(e_patch, Dict("ID"=>"1", "type"=>"everywhere"))
    for (name, value) in [("density", "0.4"), ("orientation", "random"), ("anisotropy", "0.3")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    e_layer = new_child(xml_root, "layer")
    set_attribute(e_layer, "ID", "2")

    e_patch = new_child(e_layer, "patch")
    set_attributes(e_patch, Dict("ID"=>"1", "type"=>"ellipse"))
    for (name, value) in [("x0", "200.0"), ("y0", "-200.0"), ("a", "120.0"), ("b", "60.0"), ("rotation", "pi/8"), ("thickness", "50.0"), ("density", "0.5"), ("orientation", "parallel"), ("anisotropy", "0.8")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    e_patch = new_child(e_layer, "patch")
    set_attributes(e_patch, Dict("ID"=>"2", "type"=>"elliptical_disc"))
    for (name, value) in [("x0", "200.0"), ("y0", "-200.0"), ("a", "120.0"), ("b", "60.0"), ("rotation", "pi/8"), ("density", "0.2"), ("orientation", "random"), ("anisotropy", "0.5")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    e_patch = new_child(e_layer, "patch")
    set_attributes(e_patch, Dict("ID"=>"3", "type"=>"ellipse"))
    for (name, value) in [("x0", "-200.0"), ("y0", "200.0"), ("a", "200.0"), ("b", "40.0"), ("rotation", "-pi/2"), ("thickness", "100.0"), ("density", "0.8"), ("orientation", "perpendicular"), ("anisotropy", "1.0")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    e_patch = new_child(e_layer, "patch")
    set_attributes(e_patch, Dict("ID"=>"4", "type"=>"elliptical_disc"))
    for (name, value) in [("x0", "-200.0"), ("y0", "200.0"), ("a", "200.0"), ("b", "40.0"), ("rotation", "-pi/2"), ("density", "0.6"), ("orientation", "random"), ("anisotropy", "0.7")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    save_file(xml_doc, path_to_ic_ecm_xml)
    free(xml_doc)
end

end