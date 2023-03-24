#= ======== =#
#  Plotting  #
#= ======== =#
abstract type component end

mutable struct coords
    name::String
    x::Vector{Float64}
    y::Vector{Float64}
end

function endpoints(c::coords)
    start_pt = [c.x[1], c.y[1]]
    end_pt = [c.x[end], c.y[end]]
    return start_pt, end_pt
end

Base.:+(c1::coords, c2::coords) = coords(c1.name, vcat(c1.x, c2.x), vcat(c1.y, c2.y))

# port
mutable struct sys_coords
    outer_wall::coords
    blanket_wall::coords
    inner_wall::coords
    upper_divertor::coords
    lower_divertor::coords
    inner_blanket::coords
    outer_blanket::coords
end

function sys_coords(layer_coords::Vector{coords}, structure_coords::Vector{coords})
    outer_wall = [c for c in layer_coords if c.name == "outer_wall"]
    blanket_wall = [c for c in layer_coords if c.name == "blanket_wall"]
    inner_wall = [c for c in layer_coords if c.name == "inner_wall"]
    upper_divertor = [c for c in structure_coords if c.name == "upper_divertor"]
    lower_divertor = [c for c in structure_coords if c.name == "lower_divertor"]
    inner_blanket = [c for c in structure_coords if c.name == "inner_blanket"]
    outer_blanket = [c for c in structure_coords if c.name == "outer_blanket"]
    return sys_coords(outer_wall[1], blanket_wall[1], inner_wall[1], upper_divertor[1], lower_divertor[1], inner_blanket[1], outer_blanket[1])
end

function sys_coords(dd::IMAS.dd)
    # to_find = ["hfs blanket" ,"lfs gap vacuum vessel low temp shield","lfs first wall"]
    to_find = ["hfs blanket", "hfs TF", "lfs first wall"]
    re_name = Dict(
        "hfs blanket" => "blanket_wall",
        "hfs TF" => "outer_wall",
        "lfs first wall" => "inner_wall",
        "Upper divertor" => "upper_divertor",
        "Lower divertor" => "lower_divertor",
        "LFS blanket" => "outer_blanket",
        "HFS blanket" => "inner_blanket")

    lay_coords = [coords(re_name[lay.name], lay.outline.r, lay.outline.z) for lay in dd.build.layer if (lay.name in to_find)]
    struct_coords = [coords(re_name[lay.name], lay.outline.r, lay.outline.z) for lay in dd.build.structure]
    sys_coords(lay_coords, struct_coords)
end

function sys2vec(core::sys_coords)
    return [core.outer_wall, core.blanket_wall, core.inner_wall, core.upper_divertor, core.lower_divertor, core.inner_blanket, core.outer_blanket]
end

@recipe function plot_sys(core::sys_coords)
    cvec = sys2vec(core)
    fill_col = [:gray :gray :white :mediumpurple1 :mediumpurple1 :orange :orange]
    fill_alpha = [0.5 1 1 1 1 1 1]
    count = 1
    for c in cvec
        # @show c.name
        @series begin
            seriestype --> :shape
            fillcolor --> fill_col[count]
            fillalpha --> fill_alpha[count]
            label --> nothing
            linewidth --> 1.0
            linecolor --> :black
            aspect_ratio --> :equal
            Shape(c.x, c.y)
        end
        count = count + 1
    end

    @series begin
        seriestype --> :path
        fill --> false
        label --> nothing
        linewidth --> 3
        linecolor --> :black
        aspect_ratio --> :equal
        [0.0, 0.0, 10.0, 10.0, 0], [7.5, -7.5, -7.5, 7.5, 7.5]
    end

end

function offset_entities(c::coords; dir=1.0, scale=0.1)
    ref_vec = [0.0, 0.0, dir]

    xpts = c.x
    ypts = c.y

    npts = length(xpts)

    off_x = 0
    off_y = 0

    first_point_found = false

    for i = 1:npts
        idx_a = i - 1
        idx_b = i + 1

        if idx_a < 1
            idx_a = npts
        end
        if idx_b > npts
            idx_b = 1
        end
        x_dir = xpts[idx_b] - xpts[idx_a]
        y_dir = ypts[idx_b] - ypts[idx_a]


        path_vec = [x_dir, y_dir, 0.0] ./ (sqrt(x_dir^2 + y_dir^2))

        N = cross(path_vec, ref_vec) .* scale
        newx = xpts[i] + N[1]
        newy = ypts[i] + N[2]

        min_dist = IMAS.minimum_distance_two_shapes([newx], [newy], xpts, ypts)
        if min_dist >= norm(N)
            if !first_point_found
                off_x = [newx]
                off_y = [newy]
                first_point_found = true
            else
                off_x = vcat(off_x, newx)
                off_y = vcat(off_y, newy)
            end
        end

    end
    return off_x, off_y
end

function reorder_c(c::coords; fromPt="top")
    xpts = c.x
    ypts = c.y
    idx = 1
    if fromPt == "top"
        idx = findfirst(y -> y == maximum(ypts), ypts)
    elseif fromPt == "bot"
        idx = findfirst(y -> y == minimum(ypts), ypts)
    elseif fromPt == "right"
        idx = findfirst(x -> x == maximum(xpts), xpts)
    end

    if idx != 1
        c.x = vcat(c.x[idx:end], c.x[1:idx])
        c.y = vcat(c.y[idx:end], c.y[1:idx])
    end
end

function breeder_cooling_route(sc::sys_coords)
    reorder_c(sc.outer_wall; fromPt="bot")
    ow = sc.outer_wall
    xm, ym = offset_entities(ow; dir=1.0, scale=0.2)
    c_supply = coords("off", xm, ym)

    xmr, ymr = offset_entities(c_supply; dir=1.0, scale=0.5)
    c_return = coords("blanket_return", xmr, ymr)
    reorder_c(c_return; fromPt="bot")
    miny = minimum(sc.lower_divertor.y)

    sup_idx = findall(y -> y > miny, c_supply.y)
    ret_idx = findall(y -> y > miny, c_return.y)

    newx = vcat(c_supply.x[sup_idx], reverse(c_return.x[ret_idx]))
    newy = vcat(c_supply.y[sup_idx], reverse(c_return.y[ret_idx]))
    cp = coords("breeder_cooling", newx, newy)
    reorder_c(cp)

    cp.x = cp.x[10:(end-5)]
    cp.y = cp.y[10:(end-5)]

    st, en = endpoints(cp)
    if st[1] < en[1]
        cp.x = reverse(cp.x)
        cp.y = reverse(cp.y)
    end

    st, en = endpoints(cp)
    x1, y1 = st
    x2, y2 = en

    y_st = [y1 + 0.75, y1 + 0.75, y1]
    y_en = [y2, y2 + 1.25, y1 + 1.25]
    x_st = [11, x1, x1]
    x_en = [x2, x2, 11]

    c_before = coords(cp.name, x_st, y_st)
    c_after = coords(cp.name, x_en, y_en)
    cout = c_before + cp
    cout = cout + c_after

    return cout
end

function wall_cooling_route(core::sys_coords)
    reorder_c(core.blanket_wall; fromPt="bot")
    x_start = core.blanket_wall.x[1]

    x_right = maximum(core.lower_divertor.x)
    y_tt = minimum(core.upper_divertor.y)
    y_bb = maximum(core.lower_divertor.y)

    # OUTER BLANKET WALL
    x_in, y_in = offset_entities(core.blanket_wall; dir=1.0, scale=0.35) #inner
    idx_right = findall(x -> x > x_right, x_in)
    x_ob_supply = x_in[idx_right]
    y_ob_supply = y_in[idx_right]
    ob_topIdx = findfirst(y -> y == maximum(y_ob_supply), y_ob_supply)
    ob_botIdx = findfirst(y -> y == minimum(y_ob_supply), y_ob_supply)

    x_in, y_in = offset_entities(core.blanket_wall; dir=1.0, scale=0.2) #inner
    idx_left = findall(x -> x < x_start, x_in)
    x_ib_supply = x_in[idx_left]
    y_ib_supply = y_in[idx_left]

    xloc = minimum(x_ib_supply)
    x_ib_supply = [xloc, xloc]
    y_ib_supply = [y_tt - 0.25, y_bb + 0.25]

    x_in, y_in = offset_entities(core.blanket_wall; dir=-1.0, scale=0.25) #outer
    if x_in[2] < x_in[1]
        x_in = reverse(x_in)
        y_in = reverse(y_in)
    end
    it1, it2 = IMAS.minimum_distance_two_shapes([x_ob_supply[ob_topIdx]], [y_ob_supply[ob_topIdx]], x_in, y_in; return_index=true)
    ib1, ib2 = IMAS.minimum_distance_two_shapes([x_ob_supply[ob_botIdx]], [y_ob_supply[ob_botIdx]], x_in, y_in; return_index=true)

    itt1, itt2 = IMAS.minimum_distance_two_shapes([x_ib_supply[1]], [y_ib_supply[1]], x_in, y_in; return_index=true)
    ibb1, ibb2 = IMAS.minimum_distance_two_shapes([x_ib_supply[2]], [y_ib_supply[2]], x_in, y_in; return_index=true)

    x = vcat(x_in[1:ib2], reverse(x_ob_supply), x_in[it2:itt2], x_in[itt2], x_ib_supply, x_in[ibb2], x_in[ibb2:end])
    y = vcat(y_in[1:ib2], reverse(y_ob_supply), y_in[it2:itt2], y_ib_supply[1], y_ib_supply, y_ib_supply[2], y_in[ibb2:end])
    cout = coords("wall_cooling", x, y)
    reorder_c(cout; fromPt="right")
    cout.x = cout.x[4:(end-4)]
    cout.y = cout.y[4:(end-4)]

    st, en = endpoints(cout)
    if st[2] > en[2]
        cout.x = reverse(cout.x)
        cout.y = reverse(cout.y)
    end
    st, en = endpoints(cout)
    x1, y1 = st
    x2, y2 = en

    y_st = [y1, y1]
    x_st = [11, x1]

    y_en = [y2, y2]
    x_en = [x2, 11]

    c_before = coords(cout.name, x_st, y_st)
    c_after = coords(cout.name, x_en, y_en)
    cnew = c_before + cout
    cnew = cnew + c_after
    return cnew
end

function divertor_flow_path(core::sys_coords)
    xmin = minimum(core.lower_divertor.x)
    xmax = maximum(core.lower_divertor.x)

    idxmin = findfirst(x -> x == xmax, core.lower_divertor.x)    #right most point
    yb = core.lower_divertor.y[idxmin] - 0.3

    xcriteria = findall(x -> x > xmin + 0.3, core.lower_divertor.x)
    ycriteria = findall(y -> y < yb, core.lower_divertor.y)
    goodIdx = intersect(xcriteria, ycriteria)

    xd = core.lower_divertor.x[goodIdx]
    yd = core.lower_divertor.y[goodIdx]
    xact = vcat(xd, reverse(xd))
    yact = vcat(yd .+ 0.35, reverse(yd))
    cd = coords("divflow", xact, yact)
    reorder_c(cd; fromPt="bot")

    xuse = cd.x[12:end-6]
    yuse = cd.y[12:end-6]

    cd = coords("divertor_cooling", xuse, yuse)
    st, en = endpoints(cd)
    if st[1] > en[1]
        cd.x = reverse(cd.x)
        cd.y = reverse(cd.y)
    end

    st, en = endpoints(cd)
    x1, y1 = st
    x2, y2 = en

    y_st = [y1 - 1, y1 - 1, y1]
    y_en = [y2, y2 - 0.5, y2 - 0.5]
    x_st = [11, x1, x1]
    x_en = [x2, x2, 11]

    c_before = coords(cd.name, x_st, y_st)
    c_after = coords(cd.name, x_en, y_en)
    cout = c_before + cd
    cout = cout + c_after
    return cout
end

mutable struct turb <: component
    name::String
    coords::coords
    h::Float64
    w::Float64
end

mutable struct comp <: component
    name::String
    coords::coords
    h::Float64
    w::Float64
end

mutable struct heat_exchanger <: component
    name::String
    coords::coords
    h::Float64
    w::Float64
    outline_path::coords
    hot_path::coords
    cold_path::coords
end

mutable struct intercooler <: component
    name::String
    coords::coords
    h::Float64
    w::Float64
end

function circle_coords(xcenter::Float64, ycenter::Float64, r::Float64)
    ang = LinRange(0, 2 * π, 500)
    x = xcenter .+ r .* cos.(ang)
    y = ycenter .+ r .* sin.(ang)
    return x, y
end

function turbPts(t::turb)
    half_height = t.h / 2
    half_width = t.w / 2
    y = t.coords.y[1]
    x = t.coords.x[1]

    y1 = y + half_height / 2
    y2 = y - half_height / 2
    y3 = y - half_height
    y4 = y + half_height

    x1 = x - half_width
    x2 = x + half_width

    xx = [x1, x1, x2, x2, x1]
    yy = [y1, y2, y3, y4, y1]
    return xx, yy
end

function init_hx_from_port(w::Float64, h::Float64, xport::Float64, yport::Float64)
    half_height = h / 2
    half_width = w / 2
    w_offset = half_width * 0.9
    x_center = xport + w_offset
    y_center = yport - half_height * 1.2
    return heat_exchanger(x_center, y_center, h, w)
end

function compPts(t::comp)
    half_height = t.h / 2
    half_width = t.w / 2
    y = t.coords.y[1]
    x = t.coords.x[1]

    y1 = y + half_height
    y2 = y - half_height
    y3 = y - half_height / 2
    y4 = y + half_height / 2

    x1 = x - half_width
    x2 = x + half_width

    xx = [x1, x1, x2, x2, x1]
    yy = [y1, y2, y3, y4, y1]
    return xx, yy
end

function init_hx(name::String, x::Float64, y::Float64, h::Float64, w::Float64)
    half_height = h / 2
    half_width = w / 2

    cornerX = x .+ [-half_width, -half_width, half_width, half_width]
    cornerY = y .+ [-half_height, half_height, half_height, -half_height]
    outline_path = coords("outline", cornerX, cornerY)
    #offsets for the 2 hx paths
    h_offset = half_height / 2
    w_offset = half_width * 0.9
    m_offset = half_width * 0.8
    zig_offset = half_width * 0.7

    l_port_x = x - w_offset
    l_pre_zig = x - m_offset
    left_x = [l_port_x, l_port_x, l_pre_zig]

    r_port_x = x + w_offset
    r_pre_zig = x + m_offset
    right_x = [r_pre_zig, r_port_x, r_port_x]

    t_port_y = y + half_height * 1.2
    t_path_y = y + h_offset
    top_y = [t_port_y, t_path_y, t_path_y]

    b_port_y = y - half_height * 1.2
    b_path_y = y - h_offset
    bot_y = [b_port_y, b_path_y, b_path_y]

    zig_height = half_height / 4
    zig_points = 6

    zig_xpoints = LinRange(-zig_offset, zig_offset, 7) .+ x
    idx_range = LinRange(1, zig_points + 1.0, zig_points + 1)
    zig_ypoints = real((-1 + 0im) .^ idx_range) .* zig_height

    x_top = vcat(left_x, zig_xpoints, right_x)
    y_top = vcat(top_y, zig_ypoints .+ t_path_y, reverse(top_y))
    if name == "regen"
        x_top = reverse(x_top)
        y_top = reverse(y_top)
    end
    hot_path = coords("hot_path", reverse(x_top), reverse(y_top))
    x_bot = x_top
    y_bot = vcat(bot_y, zig_ypoints .+ b_path_y, reverse(bot_y))
    if name == "regen"
        x_bot = reverse(x_bot)
        y_bot = reverse(y_bot)
    end

    cold_path = coords("cold_path", x_bot, y_bot)
    return heat_exchanger(name, coords("center_pos", [x], [y]), h, w, outline_path, hot_path, cold_path)
end

@recipe function plot_turbine(t::turb)
    xpts, ypts = turbPts(t)
    @series begin
        seriestype --> :path
        fill --> true
        fillcolor --> :gray
        fillalpha --> 1.0
        linewidth --> 1.0
        color --> :black
        label --> nothing
        annotations --> (t.coords.x, t.coords.y, ("Turb", :center, 7))
        Shape(xpts, ypts)
    end

end

@recipe function plot_compressor(c::comp)
    xpts, ypts = compPts(c)

    @series begin
        seriestype --> :path
        fill --> true
        fillcolor --> :lightseagreen
        fillalpha --> 1.0
        linewidth --> 1.0
        color --> :black
        label --> nothing
        annotations --> (c.coords.x, c.coords.y, ("comp", :center, 7))
        Shape(xpts, ypts)
    end
end

@recipe function plot_hx(hxx::heat_exchanger)
    tp_col = :red
    bp_col = :blue
    if hxx.name == "regen"
        tp_col = :blue
        bp_col = :red
    end

    @series begin
        seriestype --> :shape
        fill --> false
        fillcolor --> :white
        fillalpha --> 0.5
        linewidth --> 1.0
        color --> :black
        label --> nothing
        Shape(hxx.outline_path.x, hxx.outline_path.y)
    end
    @series begin
        seriestype --> :path
        fillcolor --> :red
        fillalpha --> 0.5
        linewidth --> 2.5
        color --> tp_col
        label --> nothing
        hxx.hot_path.x, hxx.hot_path.y
    end

    @series begin
        seriestype --> :path
        fillcolor --> :red
        fillalpha --> 0.5
        linewidth --> 2.5
        color --> bp_col
        label --> nothing
        hxx.cold_path.x, hxx.cold_path.y
    end
end

@recipe function plot_intercooler(ic::intercooler)
    xcenter = ic.coords.x
    ycenter = ic.coords.y
    half_height = ic.h / 2
    half_width = ic.w / 2

    cornerX = xcenter .+ [-half_width, -half_width, half_width, half_width]
    cornerY = ycenter .+ [-half_height, half_height, half_height, -half_height]

    @series begin
        seriestype --> :shape
        fill --> true
        fillcolor --> :lightsteelblue1
        fillalpha --> 1
        linewidth --> 1.0
        color --> :black
        label --> nothing
        annotations --> (ic.coords.x, ic.coords.y, ("cool", :center, 7))
        Shape(cornerX, cornerY)
    end
end

function attach2hx(mainc::coords, hx::heat_exchanger)
    to_connect = hx.hot_path
    mainst, mainend = endpoints(mainc)
    hxst, hxend = endpoints(to_connect)

    mpointA = coords(mainc.name, [hxst[1]], [mainend[2]])
    mpointB = coords(mainc.name, [hxend[1], mainst[1]], [mainst[2], mainst[2]])
    mainc = mainc + mpointA + to_connect + mpointB
    return mainc
end

function cyclePath(part_vec::Vector{<:component})
    ca = part_vec[1]
    cp = coords("cycle_path", ca.coords.x, ca.coords.y .- 3.0)
    regen_found = false
    regen_device = 0
    for c in part_vec
        if c == ca
            cp = cp + c.coords
        elseif typeof(c) != heat_exchanger
            mpoint = coords(cp.name, [cp.x[end]], c.coords.y)
            if typeof(c) == turb
                mpoint = coords(cp.name, c.coords.x, [cp.y[end]])
            end
            cp = cp + mpoint + c.coords
        elseif typeof(c) == heat_exchanger
            if c.name == "regen"
                mpoint = coords(cp.name, [c.hot_path.x[1]], [cp.y[end]])
                cp = cp + mpoint + c.hot_path + coords(cp.name, [c.hot_path.x[end]], [c.hot_path.y[end] + 0.5])
                regen_found = true
                regen_device = c
            else
                mpoint = coords(cp.name, [c.cold_path.x[1]], [cp.y[end]])
                cp = cp + mpoint + c.cold_path
            end
            annotate!(c.coords.x, c.coords.y .- c.h, (c.name, :center, 6))
        end
    end
    if regen_found
        st, en = endpoints(regen_device.cold_path)
        if st[1] < en[1]
            regen_device.cold_path.x = reverse(regen_device.cold_path.x)
            regen_device.cold_path.y = reverse(regen_device.cold_path.y)
        end
        mpoint = coords(cp.name, [cp.x[end], regen_device.cold_path.x[1]], [cp.y[end] - 1.5, cp.y[end] - 1.5])
        cp = cp + mpoint + regen_device.cold_path
    end

    conpoint = coords(cp.name, [cp.x[end], cp.x[1]], [cp.y[1], cp.y[1]])
    cp = cp + conpoint

    return cp
end
