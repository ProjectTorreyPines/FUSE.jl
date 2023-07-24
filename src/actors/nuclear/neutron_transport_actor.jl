using HDF5
using NeutronTransport
using GridapGmsh
using GridapGmsh: gmsh, GmshDiscreteModel
using Gridap
import NeutronTransport: MoCSolution
import Gridap: writevtk
import Gridap.Geometry: get_triangulation

const XSs = CrossSections

function get_xss_from_hdf5(
    filename::String,
    material::String,
    material_tag::String=material;
    )
    file = h5open(filename, "r")
    νΣf = read(file["material"][material]["nu-fission"]["average"])
    Σt = read(file["material"][material]["nu-transport"]["average"])
    Σs0 = permutedims(read(file["material"][material]["consistent nu-scatter matrix"]["average"]))
    return CrossSections(material_tag, length(νΣf); νΣf, Σt, Σs0)
end

function get_xss_from_hdf5(filename::String, material::String, xs_type::Vector{String})
    xs_dict = Dict{String,Vector{Float64}}()
    file = h5open(filename, "r")
    for xs_name in xs_type
        xs = read(file["material"][material][xs_name]["average"])
        xs_dict[xs_name] = xs
    end
    close(file)
    return xs_dict
end

function get_xss_from_hdf5(filename::String, materials::Vector{String}, xs_type::Vector{String})
    xs_dict=Dict{String, Dict{String,Vector{Float64}}}()
    for mat in materials
        mat_xs_dict = get_xss_from_hdf5(filename, mat, xs_type)
        xs_dict[mat]=mat_xs_dict
    end
    return xs_dict
end

function concentric_circles(layer_thicknesses::Vector{T}, materials::Vector{String}, layers::Vector{String}, data_path::String, data_filename::String, save_path::String=pwd(), save_filename::String="concentric_circles") where T
    gmsh.initialize()
    model_bound = sum(layer_thicknesses) + 1

    mats = deepcopy(materials)
    thicknesses = deepcopy(layer_thicknesses)
    names = deepcopy(layers)
    push!(mats, "Air-(dry,-near-sea-level)")
    push!(names, "Exterior")

    append!(thicknesses, 1.0)
    factory = gmsh.model.geo
    lc = model_bound/30

    radius=0.0

    data_filename=joinpath(data_path, data_filename)

    for (idx, thickness) in enumerate(thicknesses)
        radius+=thickness
        mats[idx] = replace(mats[idx], " " => "-")
        if mats[idx] == "Vacuum"
            mats[idx] = "Air-(dry,-near-sea-level)"
        end 
        material = factory.addPhysicalGroup(2, [idx], idx)
        GridapGmsh.gmsh.model.setPhysicalName(2, idx, names[idx])

        if idx == length(mats)
            # square cell
            top_right = factory.addPoint(radius, radius, 0, lc, 5*(idx-1)+1)
            bottom_right = factory.addPoint(radius, -radius, 0, lc, 5*(idx-1)+2)
            bottom_left = factory.addPoint(-radius, -radius, 0, lc, 5*(idx-1)+3)
            top_left = factory.addPoint(-radius, radius, 0, lc, 5*(idx-1)+4)
            right = factory.addLine(top_right, bottom_right,  4*(idx-1)+1)
            bottom = factory.addLine(bottom_right, bottom_left, 4*(idx-1)+2)
            left = factory.addLine(bottom_left, top_left, 4*(idx-1)+3)
            top = factory.addLine(top_left, top_right, 4*(idx-1)+4)
            boundary = factory.addCurveLoop([right, bottom, left, top], idx)

            # boundaries
            right_bound = factory.addPhysicalGroup(1,  [right], material+1)
            bottom_bound = factory.addPhysicalGroup(1, [bottom], material+2)
            left_bound = factory.addPhysicalGroup(1, [left], material+3)
            top_bound = factory.addPhysicalGroup(1, [top], material+4)
            GridapGmsh.gmsh.model.setPhysicalName(1, right_bound, "right")
            GridapGmsh.gmsh.model.setPhysicalName(1, bottom_bound, "bottom")
            GridapGmsh.gmsh.model.setPhysicalName(1, left_bound, "left")
            GridapGmsh.gmsh.model.setPhysicalName(1, top_bound, "top")
        else
            # make inner circle points
            center = factory.addPoint(0, 0, 0, lc, 5*(idx-1)+1)
            right = factory.addPoint(radius, 0, 0, lc, 5*(idx-1)+2)
            top = factory.addPoint(0, radius, 0, lc, 5*(idx-1)+3)
            left = factory.addPoint(-radius, 0, 0, lc, 5*(idx-1)+4)
            bottom = factory.addPoint(0, -radius, 0, lc, 5*(idx-1)+5)

            # make arcs and circle
            right_top = factory.addCircleArc(right, center, top, 4*(idx-1)+1)
            top_left = factory.addCircleArc(top, center, left, 4*(idx-1)+2)
            left_bottom = factory.addCircleArc(left, center, bottom, 4*(idx-1)+3)
            bottom_right = factory.addCircleArc(bottom, center, right, 4*(idx-1)+4)
            circle = factory.addCurveLoop([right_top, top_left, left_bottom, bottom_right], idx)
        end

        # make surfaces and materials
        if idx==1
            surface = factory.addPlaneSurface([idx], idx)
        else
            surface = factory.addPlaneSurface([idx, idx-1], idx)
        end
    end

    xss = [get_xss_from_hdf5(data_filename, mats[idx], names[idx]) for idx in eachindex(thicknesses)]

    factory.synchronize()

    gmsh.model.mesh.generate(2)

    mshfile = joinpath(save_path, save_filename * ".msh")

    gmsh.write(mshfile)

    # gmsh.fltk.run()

    gmsh.finalize()

    model = GmshDiscreteModel(mshfile; renumber=true)
    jsonfile = joinpath(save_path, save_filename * ".json")
    Gridap.Io.to_json_file(model, jsonfile)

    jsonfile = joinpath(save_path, save_filename * ".json")
    geometry = DiscreteModelFromFile(jsonfile)

    # number of azimuthal angles
    nφ = 4

    # azimuthal spacing
    δ = 3.0 

    # boundary conditions
    bcs = BoundaryConditions(top=Vaccum, left=Vaccum, bottom=Vaccum, right=Vaccum)

    # initialize track generator
    tg = TrackGenerator(geometry, nφ, δ, bcs=bcs)
    
    # perform ray tracing
    trace!(tg)

    # proceed to segmentation
    segmentize!(tg)

    # polar quadrature
    pq = NeutronTransport.TabuchiYamamoto(2)

    # define the problem
    prob = MoCProblem(tg, pq, xss)

    # define fixed source material
    fixed_sources = set_fixed_source_material(prob, "plasma", 8 , 1)

    # solve
    sol = NeutronTransport.solve(prob, fixed_sources, debug=true, max_residual=0.05, max_iterations=500)

    return sol
end

"""
    get_cell_φ(sol, cell, g)

Returns the average φ in a given cell and energy group.
"""
function get_cell_φ(sol::MoCSolution{T}, cell::Int, g::Int) where T
    φs = sol(g)[findall(c -> c == cell, sol.prob.fsr_tag)]
    # vols = sol.prob.trackgenerator.volumes[findall(c -> c == cell, sol.prob.fsr_tag)]
    
    # vol_integrated_φ = sum(φs .* vols / 31415)
    # return vol_integrated_φ

    mean_φ = sum(φs)/length(φs)
    return mean_φ
end

function get_cell_φ(sol::MoCSolution{T}, cell::Int, gs::UnitRange{Int64}) where T
    φs = Float64[]
    for g in gs
        append!(φs, get_cell_φ(sol, cell, g))
    end
    return φs
end

function get_tbr(sol::MoCSolution{T}, data_path::String, data_filename::String) where T
    NeutronTransport.@unpack φ, prob = sol
    NeutronTransport.@unpack trackgenerator, fsr_tag, xss = prob
    NeutronTransport.@unpack volumes = trackgenerator
    NGroups = NeutronTransport.ngroups(prob)
    NRegions = NeutronTransport.nregions(prob)

    vol_integrated_φ = Vector{Float64}()
    tbr=0.0
    plasma_vol = sum(volumes[findall(t -> t == 1, fsr_tag)])

    for i in 1:NRegions
        tbr_xs_dict = get_xss_from_hdf5(joinpath(data_path,data_filename), xss[fsr_tag[i]].name, ["(n,Xt)"])
        for g in 1:NGroups
            ig = NeutronTransport.@region_index(i, g)
            push!(vol_integrated_φ, sol.φ[ig] * volumes[i] / plasma_vol)
            tbr+=vol_integrated_φ[ig]*tbr_xs_dict["(n,Xt)"][g]
        end        
    end
    return tbr
end

function write_neutron_flux_vtk(sol::MoCSolution{T}, groups::Union{Vector{Int},UnitRange{Int}}) where T
    NeutronTransport.@unpack prob = sol
    NeutronTransport.@unpack trackgenerator = prob
    triang = get_triangulation(trackgenerator.mesh.model)
    writevtk(triang, "fluxes", cellfields=[string(g) => sol(g) for g in groups])
end

function neutron_transport_1d(dd::IMAS.dd, data_path::String, data_filename::String, save_path::String=pwd(), save_filename::String="fuse_neutron_transport")
    thicknesses = [layer.thickness*100 for layer in dd.build.layer if -1 <= layer.fs <= 2]
    materials = [layer.material for layer in dd.build.layer if -1 <= layer.fs <= 2]
    layers = [layer.name for layer in dd.build.layer if -1 <= layer.fs <= 2]
    return concentric_circles(thicknesses, materials, layers, data_path, data_filename, save_path, save_filename)
end