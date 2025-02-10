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

function parseLayer!(df::AbstractDataFrame, layer::XMLElement, config_dict::Dict{String,Float64})
    layer_df = initializeDataFrame(config_dict)
    for patch_collection in child_elements(layer)
        parsePatchCollection!(layer_df, patch_collection, config_dict)
    end
    updateDataFrame!(df, layer_df; allow_overwrite=true)
end

function parsePatchCollection!(layer_df::AbstractDataFrame, patch_collection::XMLElement, config_dict::Dict{String,Float64})
    patch_type = attribute(patch_collection, "type")
    if patch_type == "everywhere"
        patch_parser = parseEverywherePatch
    elseif patch_type == "ellipse"
        patch_parser = parseEllipsePatch
    elseif patch_type == "elliptical_disc"
        patch_parser = parseEllipticalDiscPatch
    elseif patch_type == "ellipse_with_shell"
        patch_parser = parseEllipseWithShellPatch
    else
        throw(ArgumentError("Patch type $patch_type not recognized."))
    end
    for patch in child_elements(patch_collection)
        patch_df = patch_parser(patch, config_dict)
        updateDataFrame!(layer_df, patch_df; allow_overwrite=false)
    end
end

function parseEverywherePatch(patch::XMLElement, config_dict::Dict{String,Float64})
    df = initializeDataFrame(config_dict)
    parseEverywherePatch!(df, patch)
    return df
end

function parseEverywherePatch!(df::AbstractDataFrame, patch::XMLElement)
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
end

function parseEllipsePatch(patch::XMLElement, config_dict::Dict{String,Float64})
    df = initializeDataFrame(config_dict)
    parseEllipsePatch!(df, patch)
    return df
end

function parseEllipsePatch!(df::AbstractDataFrame, patch::XMLElement; check_coords::Bool=true)
    e = parseEllipseParameters(patch)

    if check_coords
        in_ellipse = [insideEllipseThickness(e, x, y) for (x, y) in zip(df.x, df.y)]
    else
        in_ellipse = trues(length(df.x))
    end
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
end

function parseEllipticalDiscPatch(patch::XMLElement, config_dict::Dict{String,Float64})
    df = initializeDataFrame(config_dict)
    parseEllipticalDiscPatch!(df, patch)
    return df
end

function parseEllipticalDiscPatch!(df::AbstractDataFrame, patch::XMLElement; check_coords::Bool=true)
    e = parseEllipseParameters(patch)

    if check_coords
        in_disc = [insideEllipse(e, x, y) for (x, y) in zip(df.x, df.y)]
    else
        in_disc = trues(length(df.x))
    end
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
end

function parseEllipseWithShellPatch(patch::XMLElement, config_dict::Dict{String,Float64})
    df = initializeDataFrame(config_dict)

    e = parseEllipseParameters(patch)
    positions = [positionRelativeToEllipse(e, x, y) for (x, y) in zip(df.x, df.y)]

    interior_element = find_element(patch, "interior")
    if isnothing(interior_element)
        throw(ErrorException("An interior element must be specified for an ellipse with shell ECM patch.\nIf only setting the shell of the ellipse, use `ellipse` instead of `ellipse_with_shell`."))
    end
    for (tag, content_) in [("x0", e.center[1]), ("y0", e.center[2]), ("a", e.axes[1]), ("b", e.axes[2]), ("rotation", e.rotation)]
        new_element = new_child(interior_element, tag)
        set_content(new_element, string(content_))
    end
    interior_df = @view df[positions.==:inside, :]
    parseEllipticalDiscPatch!(interior_df, interior_element; check_coords=false)

    shell_element = find_element(patch, "shell")
    if isnothing(shell_element)
        throw(ErrorException("A shell element must be specified for an ellipse with shell ECM patch.\nIf only setting the interior of the ellipse, use `elliptical_disc` instead of `ellipse_with_shell`."))
    end
    for (tag, content_) in [("x0", e.center[1]), ("y0", e.center[2]), ("a", e.axes[1] + e.thickness), ("b", e.axes[2] + e.thickness), ("rotation", e.rotation), ("thickness", e.thickness)]
        new_element = new_child(shell_element, tag)
        set_content(new_element, string(content_))
    end
    shell_df = @view df[positions.==:thickness, :]
    parseEllipsePatch!(shell_df, shell_element; check_coords=false)

    exterior_element = find_element(patch, "exterior")
    if isnothing(exterior_element)
        return df
    end

    # optional exterior element defining ECM outside the ellipse with shell
    in_exterior = positions.==:outside
    n = sum(in_exterior)
    density, orientation, anisotropy = parseECMFeatures(exterior_element)
    df.ecm_density[in_exterior] .= density
    if orientation == "random"
        theta = 2 * pi * rand(n)
        df.ecm_orientation_x[in_exterior] .= anisotropy * cos.(theta)
        df.ecm_orientation_y[in_exterior] .= anisotropy * sin.(theta)
    else
        throw(ArgumentError("Orientation $orientation not recognized for the exterior an ellipse with shell ECM patch.\nRecognized orientations are: `random`"))
    end
    df.is_missing[in_exterior] .= false
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

function updateDataFrame!(df::AbstractDataFrame, new_df::AbstractDataFrame; allow_overwrite::Bool=false)
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
The optional `monolayer` keyword argument will create a monolayer template instead of the default multilayer template.
The monolayer template creates the ECM based on relative position to an ellipse with a shell around it.
"""
function createICECMXMLTemplate(path_to_folder::String; monolayer::Bool=false)

    path_to_ic_ecm_xml = joinpath(path_to_folder, "ecm.xml")
    mkpath(dirname(path_to_ic_ecm_xml))
    xml_doc = XMLDocument()
    xml_root = create_root(xml_doc, "ic_ecm")

    xml_root |> (monolayer ? createMonolayerXMLTemplate! : createMultilayerXMLTemplate)

    save_file(xml_doc, path_to_ic_ecm_xml)
    free(xml_doc)
end

function createMonolayerXMLTemplate!(xml_root::XMLElement)
    e_layer = new_child(xml_root, "layer")
    set_attribute(e_layer, "ID", "1")

    e_ellipse_with_shell = new_child(e_layer, "patch_collection")
    set_attribute(e_ellipse_with_shell, "type", "ellipse_with_shell")

    e_patch = new_child(e_ellipse_with_shell, "patch")
    set_attribute(e_patch, "ID", "1")

    for (name, value) in [("x0", "0.0"), ("y0", "0.0"), ("a", "400.0"), ("b", "200.0"), ("rotation", "π/6"), ("thickness", "100.0")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    for tag in ["interior", "shell", "exterior"]
        e_location = new_child(e_patch, tag)
        for (name, value) in [("density", "0.4"), ("orientation", tag == "shell" ? "perpendicular" : "random"), ("anisotropy", "0.3")]
            e = new_child(e_location, name)
            set_content(e, value)
        end
    end
end

function createMultilayerXMLTemplate(xml_root::XMLElement)
    e_layer = new_child(xml_root, "layer")
    set_attribute(e_layer, "ID", "1")

    ## layer 1 everyhweres
    e_everyhwere = new_child(e_layer, "patch_collection")
    set_attribute(e_everyhwere, "type", "everywhere")

    e_patch = new_child(e_everyhwere, "patch")
    set_attributes(e_patch, Dict("ID" => "1"))
    for (name, value) in [("density", "0.4"), ("orientation", "random"), ("anisotropy", "0.3")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    # layer 2
    e_layer = new_child(xml_root, "layer")
    set_attribute(e_layer, "ID", "2")

    ## layer 2 ellipses
    e_ellipse = new_child(e_layer, "patch_collection")
    set_attribute(e_ellipse, "type", "ellipse")

    e_patch = new_child(e_ellipse, "patch")
    set_attributes(e_patch, Dict("ID" => "1"))
    for (name, value) in [("x0", "200.0"), ("y0", "-200.0"), ("a", "120.0"), ("b", "60.0"), ("rotation", "pi/8"), ("thickness", "50.0"), ("density", "0.5"), ("orientation", "parallel"), ("anisotropy", "0.8")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    e_patch = new_child(e_ellipse, "patch")
    set_attributes(e_patch, Dict("ID" => "2"))
    for (name, value) in [("x0", "-200.0"), ("y0", "200.0"), ("a", "200.0"), ("b", "40.0"), ("rotation", "-pi/2"), ("thickness", "100.0"), ("density", "0.8"), ("orientation", "perpendicular"), ("anisotropy", "1.0")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    ## layer 2 elliptical discs
    e_elliptical_disc = new_child(e_layer, "patch_collection")
    set_attribute(e_elliptical_disc, "type", "elliptical_disc")

    e_patch = new_child(e_elliptical_disc, "patch")
    set_attributes(e_patch, Dict("ID" => "1"))
    for (name, value) in [("x0", "200.0"), ("y0", "-200.0"), ("a", "120.0"), ("b", "60.0"), ("rotation", "pi/8"), ("density", "0.2"), ("orientation", "random"), ("anisotropy", "0.5")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    e_patch = new_child(e_elliptical_disc, "patch")
    set_attributes(e_patch, Dict("ID" => "2"))
    for (name, value) in [("x0", "-200.0"), ("y0", "200.0"), ("a", "200.0"), ("b", "40.0"), ("rotation", "-pi/2"), ("density", "0.6"), ("orientation", "random"), ("anisotropy", "0.7")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    ## layer 2 ellipses with shell
    e_ellipse_with_shell = new_child(e_layer, "patch_collection")
    set_attribute(e_ellipse_with_shell, "type", "ellipse_with_shell")

    e_patch = new_child(e_ellipse_with_shell, "patch")
    set_attribute(e_patch, "ID", "1")
    for (name, value) in [("x0", "350.0"), ("y0", "350.0"), ("a", "100.0"), ("b", "100.0"), ("rotation", "0.0"), ("thickness", "80.0")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end
    e_interior = new_child(e_patch, "interior")
    for (name, value) in [("density", "0.4"), ("orientation", "random"), ("anisotropy", "0.3")]
        e = new_child(e_interior, name)
        set_content(e, value)
    end
    e_shell = new_child(e_patch, "shell")
    for (name, value) in [("density", "0.8"), ("orientation", "random"), ("anisotropy", "0.3")]
        e = new_child(e_shell, name)
        set_content(e, value)
    end
end
end