#= ========== =#
#  ActorZMQ    #
#= ========== =#
# Coupling actor for FUSE <-> GSLite/GSEvolve communication via ZeroMQ (REQ/REP pattern)
#
# In the dynamic plasma loop:
#   Phase 1 start: receive!(actor_zmq) — gets flux matrix + Ip from GSLite
#   Phase 2 end:   send!(actor_zmq)    — sends FRESCO equilibrium back to GSLite
#
# Message format: JSON over ZMQ REQ/REP
#
# Requires: ZMQ.jl (add to Project.toml: `import Pkg; Pkg.add("ZMQ")`)

import JSON
import Interpolations
import ZMQ

@actor_parameters_struct ActorZMQ{T} begin
    endpoint::Entry{String} = Entry{String}("-", "ZMQ endpoint (e.g., tcp://localhost:5555)"; default="tcp://localhost:5555")
    timeout_ms::Entry{Int} = Entry{Int}("ms", "ZMQ receive timeout"; default=30000)
    enabled::Entry{Bool} = Entry{Bool}("-", "Enable ZMQ coupling with external code"; default=false)
    nR::Entry{Int} = Entry{Int}("-", "Number of R grid points for psizr reshape"; default=33)
    nZ::Entry{Int} = Entry{Int}("-", "Number of Z grid points for psizr reshape"; default=33)
end

mutable struct ActorZMQ{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorZMQ{P}}
    context::Any   # ZMQ.Context when connected, nothing otherwise
    socket::Any    # ZMQ.Socket when connected, nothing otherwise
    is_connected::Bool
    first_receive::Bool
    # Previous-step values for time derivatives in send!
    prev_time::Float64
    prev_betap::Float64
    prev_li::Float64
    prev_psipla::Float64
end

function ActorZMQ(dd::IMAS.dd{D}, par::FUSEparameters__ActorZMQ{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorZMQ)
    par = OverrideParameters(par; kw...)
    return ActorZMQ{D,P}(dd, par, nothing, nothing, false, true, NaN, NaN, NaN, NaN)
end

"""
    ActorZMQ(dd::IMAS.dd, act::ParametersAllActors; kw...)

Coupling actor for exchanging data with external codes (e.g., GSLite/GSEvolve) via ZeroMQ.

This actor uses a REQ/REP pattern where FUSE initiates each exchange:
- `receive!`: FUSE sends "ready" → external code replies with DataForFUSE
- `send!`: FUSE sends DataFromFUSE → external code replies with "ack"

Intended integration in ActorDynamicPlasma:
- Phase 1 start: `receive!(actor_zmq)` — receive ψ(R,Z), Ip, NBI, gas, coil currents
- Phase 2 end: `send!(actor_zmq)` — send betap, li, p_res, derivatives, CO2 density
"""
function ActorZMQ(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorZMQ(dd, act.ActorZMQ; kw...)
    if actor.par.enabled
        connect!(actor)
    end
    return actor
end

#= ============ =#
#  Connection    #
#= ============ =#

"""
    connect!(actor::ActorZMQ)

Establish ZMQ REQ socket connection to external code.
"""
function connect!(actor::ActorZMQ)
    if actor.is_connected
        return actor
    end

    @info "ActorZMQ: connecting to $(actor.par.endpoint)"
    ctx = ZMQ.Context()
    sock = ZMQ.Socket(ctx, ZMQ.REQ)
    ZMQ.connect(sock, string(actor.par.endpoint))
    actor.context = ctx
    actor.socket = sock
    actor.is_connected = true
    return actor
end

"""
    disconnect!(actor::ActorZMQ)

Close ZMQ socket and context.
"""
function disconnect!(actor::ActorZMQ)
    if !actor.is_connected
        return actor
    end
    @info "ActorZMQ: disconnecting"

    ZMQ.close(actor.socket)
    ZMQ.close(actor.context)
    actor.socket = nothing
    actor.context = nothing
    actor.is_connected = false
    return actor
end

#= ===================== =#
#  Receive from GSLite    #
#= ===================== =#

"""
    receive!(actor::ActorZMQ)

Receive data from GSLite and update dd.
Called at the start of Phase 1 in the dynamic plasma loop.

DataForFUSE fields (matching C++ struct):
- `sim_time`:        double          — GSLite simulation time [s]
- `Ip_latest`:       double          — Latest Ip measurement [A] → pulse_schedule
- `Ip_avg`:          double          — Average Ip measurement [A] → dd._aux (for NN ne predictor)
- `Bt`:              double          — Vacuum toroidal field at R0 [T] (first message only)
- `pr15v`:           double          — PR15V Ip measurement [A] → dd._aux (for NN ne predictor)
- `I_coil`:          double[NCC_MEAS] — PF coil currents [A] → dd._aux (for NN ne predictor)
- `psizr`:           double[NGG]      — Flat flux matrix ψ(R,Z) [Wb/rad], reshaped to (nR, nZ)
- `pinj_per_beam`:   double[NNBI]     — NBI injected power per beam [W] → pulse_schedule.nbi
- `nbi_acc_voltage`: double[NNBI]     — NBI acceleration voltage per beam [eV] → pulse_schedule.nbi
- `gas_cal`:         double[NGAS]     — Gas calibration values → dd._aux (for NN ne predictor)
"""
function receive!(actor::ActorZMQ)
    if !actor.par.enabled || !actor.is_connected
        return actor
    end

    dd = actor.dd
    par = actor.par


    # REQ: send ready signal, then receive data
    ZMQ.send(actor.socket, JSON.json(Dict("status" => "ready", "time" => dd.global_time)))
    raw = String(ZMQ.recv(actor.socket))
    msg = JSON.parse(raw)

    @info "ActorZMQ: received data at sim_time=$(get(msg, "sim_time", "?")) s"

    # --- Update Ip in pulse_schedule so QED/FRESCO pick it up ---
    if haskey(msg, "Ip_latest")
        ip = Float64(msg["Ip_latest"])
        ps_fc = dd.pulse_schedule.flux_control
        IMAS.set_time_array(ps_fc.i_plasma, :reference, dd.global_time, ip)
    end

    # --- Store auxiliary Ip signals for NN ne predictor ---
    # Uses dd._aux (same pattern as FUSE workflow/logging metadata)
    aux = getfield(dd, :_aux)
    if haskey(msg, "Ip_avg")
        if :zmq_Ip_avg ∉ keys(aux)
            aux[:zmq_Ip_avg] = Dict{Float64,Float64}()
        end
        aux[:zmq_Ip_avg][dd.global_time] = Float64(msg["Ip_avg"])
    end
    if haskey(msg, "pr15v")
        if :zmq_pr15v ∉ keys(aux)
            aux[:zmq_pr15v] = Dict{Float64,Float64}()
        end
        aux[:zmq_pr15v][dd.global_time] = Float64(msg["pr15v"])
    end

    # --- Update Bt (first step only, constant per shot) ---
    if haskey(msg, "Bt")
        b0 = Float64(msg["Bt"])
        if length(dd.equilibrium.vacuum_toroidal_field.b0) == 0
            push!(dd.equilibrium.vacuum_toroidal_field.b0, b0)
        else
            dd.equilibrium.vacuum_toroidal_field.b0[end] = b0
        end
        @info "ActorZMQ: set Bt = $b0 T"
    end

    # --- Store PF coil currents from GSLite (isolated from dd.pf_active) ---
    # Stored in dd._aux only so FRESCO's own coil current solve is not disturbed.
    # NN ne predictor reads from aux[:zmq_I_coil].
    if haskey(msg, "I_coil")
        i_coil = Float64.(msg["I_coil"])
        if :zmq_I_coil ∉ keys(aux)
            aux[:zmq_I_coil] = Dict{Float64,Vector{Float64}}()
        end
        aux[:zmq_I_coil][dd.global_time] = i_coil
    end

    # --- Update NBI power per beam [W] ---
    # Uses IMAS.set_time_array to properly build up time history,
    # which is needed by smooth_beam_power (smooths over τ_thermalization)
    if haskey(msg, "pinj_per_beam")
        pinj = Float64.(msg["pinj_per_beam"])
        ps_nbi = dd.pulse_schedule.nbi
        n_beams = min(length(pinj), length(ps_nbi.unit))
        time0 = dd.global_time
        for k in 1:n_beams
            IMAS.set_time_array(ps_nbi.unit[k].power, :reference, time0, pinj[k])
        end
        @info "ActorZMQ: updated NBI power for $n_beams beams"
    end

    # --- Update NBI acceleration voltage per beam [eV] ---
    if haskey(msg, "nbi_acc_voltage")
        nbi_voltage = Float64.(msg["nbi_acc_voltage"])
        ps_nbi = dd.pulse_schedule.nbi
        n_beams = min(length(nbi_voltage), length(ps_nbi.unit))
        time0 = dd.global_time
        for k in 1:n_beams
            IMAS.set_time_array(ps_nbi.unit[k].energy, :reference, time0, nbi_voltage[k])
        end
        @info "ActorZMQ: updated NBI voltage for $n_beams beams"
    end

    # --- Store gas calibration values for NN ne predictor ---
    # gas_cal[NGAS]: index 1=gasA, 2=gasB, 3=gasC, 4=gasD, 5=gasE
    if haskey(msg, "gas_cal")
        gas_cal = Float64.(msg["gas_cal"])
        gas_names = (:zmq_gasa_cal, :zmq_gasb_cal, :zmq_gasc_cal, :zmq_gasd_cal, :zmq_gase_cal)
        for k in 1:min(length(gas_cal), length(gas_names))
            sym = gas_names[k]
            if sym ∉ keys(aux)
                aux[sym] = Dict{Float64,Float64}()
            end
            aux[sym][dd.global_time] = gas_cal[k]
        end
    end

    # --- Update equilibrium from flat psizr array ---
    if haskey(msg, "psizr")
        psizr_flat = Float64.(msg["psizr"])
        nR = par.nR
        nZ = par.nZ
        @assert length(psizr_flat) == nR * nZ "ActorZMQ: psizr length $(length(psizr_flat)) != nR*nZ = $(nR*nZ)"
        psi_rz = reshape(psizr_flat, nR, nZ)  # GSLite stores column-major (R varies fastest)

        eqt = dd.equilibrium.time_slice[]

        # Ensure profiles_2d exists
        if isempty(eqt.profiles_2d)
            resize!(eqt.profiles_2d, 1)
        end
        p2d = eqt.profiles_2d[1]

        # Use existing grid from dd if available, otherwise need r_grid/z_grid from GSLite
        if haskey(msg, "r_grid") && haskey(msg, "z_grid")
            p2d.grid.dim1 = Float64.(msg["r_grid"])
            p2d.grid.dim2 = Float64.(msg["z_grid"])
        elseif isempty(p2d.grid.dim1)
            @warn "ActorZMQ: psizr received but no R/Z grid available; skipping equilibrium update"
            return actor
        end
        p2d.psi = psi_rz
        p2d.grid_type.index = 1  # rectangular grid

        rgrid = range(p2d.grid.dim1[1], p2d.grid.dim1[end], length=length(p2d.grid.dim1))
        zgrid = range(p2d.grid.dim2[1], p2d.grid.dim2[end], length=length(p2d.grid.dim2))
        fw_r, fw_z = IMAS.first_wall(dd.wall)

        if actor.first_receive 
            # First step: full flux_surfaces (Method 2) to get all 1D profiles for QED
            # First find axis and boundary from psizr (Method 1)
            PSI_itp = Interpolations.cubic_spline_interpolation(
                (rgrid, zgrid), psi_rz;
                extrapolation_bc=Interpolations.Line())
            psi_sign = sign(PSI_itp(rgrid[1], zgrid[1]) - PSI_itp((rgrid[1]+rgrid[end])/2, (zgrid[1]+zgrid[end])/2))
            axis_result = IMAS.find_magnetic_axis(rgrid, zgrid, PSI_itp, psi_sign)
            Ψaxis = PSI_itp(axis_result.RA, axis_result.ZA)
            axis2bnd = psi_sign > 0 ? :increasing : :decreasing
            psi_bnd = IMAS.find_psi_boundary(
                rgrid, zgrid, psi_rz, Ψaxis, axis2bnd, axis_result.RA, axis_result.ZA, fw_r, fw_z;
                raise_error_on_not_open=false, raise_error_on_not_closed=false)
            Ψbnd = psi_bnd.last_closed

            # Set up 1D seed arrays for flux_surfaces
            eqt1d = eqt.profiles_1d
            n_psi = 101
            eqt1d.psi = collect(range(Ψaxis, Ψbnd, length=n_psi))
            eqt1d.f = fill(eqt.global_quantities.vacuum_toroidal_field.b0 * eqt.global_quantities.vacuum_toroidal_field.r0, n_psi)
            eqt1d.pressure = zeros(n_psi)
            eqt1d.f_df_dpsi = zeros(n_psi)

            # Set global quantities
            eqt.global_quantities.magnetic_axis.r = axis_result.RA
            eqt.global_quantities.magnetic_axis.z = axis_result.ZA
            eqt.global_quantities.psi_axis = Ψaxis
            eqt.global_quantities.psi_boundary = Ψbnd

            # Run full flux_surfaces to get all 1D profiles (gm1, gm9, q, volume, etc.)
            IMAS.flux_surfaces(eqt, fw_r, fw_z)
            @info "ActorZMQ: first step — full flux_surfaces from psizr ($(nR)×$(nZ))"
        else
            # Subsequent steps: extract boundary only (Method 1), FRESCO handles the rest
            PSI_itp = Interpolations.cubic_spline_interpolation(
                (rgrid, zgrid), psi_rz;
                extrapolation_bc=Interpolations.Line())
            psi_sign = sign(PSI_itp(rgrid[1], zgrid[1]) - PSI_itp((rgrid[1]+rgrid[end])/2, (zgrid[1]+zgrid[end])/2))
            axis_result = IMAS.find_magnetic_axis(rgrid, zgrid, PSI_itp, psi_sign)
            Ψaxis = PSI_itp(axis_result.RA, axis_result.ZA)
            axis2bnd = psi_sign > 0 ? :increasing : :decreasing
            psi_bnd = IMAS.find_psi_boundary(
                rgrid, zgrid, psi_rz, Ψaxis, axis2bnd, axis_result.RA, axis_result.ZA, fw_r, fw_z;
                raise_error_on_not_open=false, raise_error_on_not_closed=false)
            Ψbnd = psi_bnd.last_closed

            # Trace LCFS boundary
            psi_levels = Float64[Ψaxis, Ψbnd]
            surfaces = IMAS.trace_simple_surfaces(psi_levels, rgrid, zgrid, psi_rz, PSI_itp,
                axis_result.RA, axis_result.ZA, fw_r, fw_z)

            # Update boundary in equilibrium
            if !isempty(surfaces)
                eqt.boundary.outline.r = surfaces[end].r
                eqt.boundary.outline.z = surfaces[end].z
            end
            eqt.global_quantities.magnetic_axis.r = axis_result.RA
            eqt.global_quantities.magnetic_axis.z = axis_result.ZA
            eqt.global_quantities.psi_axis = Ψaxis
            eqt.global_quantities.psi_boundary = Ψbnd

            @info "ActorZMQ: updated boundary from psizr ($(nR)×$(nZ))"
        end
        actor.first_receive = false
    end

    return actor
end

#= ================== =#
#  Send to GSLite      #
#= ================== =#

"""
    send!(actor::ActorZMQ)

Send FUSE results back to GSLite after Phase 2.
Called at the end of Phase 2 in the dynamic plasma loop.

DataFromFUSE fields (matching C++ struct):
- `sim_time`:     double        — Current simulation time [s]
- `betap`:        double        — Poloidal beta (beta_pol)
- `betap_dot`:    double        — Time derivative of beta_pol [1/s]
- `li`:           double        — Internal inductance (li_1)
- `li_dot`:       double        — Time derivative of li_1 [1/s]
- `p_res`:        double        — Plasma resistance [Ohm] (circuit-model: dψ_plasma/dt / Ip)
- `dens_co2_sig`: double[NCO2]  — CO2 interferometer line-integrated density [m⁻²]
"""
function send!(actor::ActorZMQ)
    if !actor.par.enabled || !actor.is_connected
        return actor
    end

    dd = actor.dd
    eqt = dd.equilibrium.time_slice[]


    betap = eqt.global_quantities.beta_pol
    li = eqt.global_quantities.li_1
    time_now = dd.global_time

    # Time derivatives (0.0 on first step when prev values are NaN)
    dt = time_now - actor.prev_time
    if isnan(actor.prev_time) || dt <= 0.0
        betap_dot = 0.0
        li_dot = 0.0
        p_res = 1e-9
    else
        betap_dot = (betap - actor.prev_betap) / dt
        li_dot = (li - actor.prev_li) / dt
        # Plasma resistance: Rp = dψ_plasma/dt / Ip, clamped >= 1e-9
        psipla_now = _compute_psipla(eqt)
        Ip_val = eqt.global_quantities.ip
        p_res = max((psipla_now - actor.prev_psipla) / dt / Ip_val, 1e-9)
    end

    msg = Dict{String,Any}()
    msg["sim_time"] = time_now
    msg["betap"] = betap
    msg["betap_dot"] = betap_dot
    msg["li"] = li
    msg["li_dot"] = li_dot
    msg["p_res"] = p_res
    msg["dens_co2_sig"] = _compute_co2_density(dd)

    # REQ: send data, then receive acknowledgment
    ZMQ.send(actor.socket, JSON.json(msg))
    ack = String(ZMQ.recv(actor.socket))
    @info "ActorZMQ: sent betap=$betap, li=$li, p_res=$p_res at t=$(time_now) s"

    # Store current values for next step's derivatives
    actor.prev_time = time_now
    actor.prev_betap = betap
    actor.prev_li = li
    actor.prev_psipla = _compute_psipla(eqt)

    return actor
end

#= ========================= =#
#  Standard actor interface   #
#= ========================= =#


#= ========== =#
#  Utilities   #
#= ========== =#

"""
    _compute_psipla(eqt)

Current-density-weighted average poloidal flux:
    ψ_plasma = ∫(ψ × j_tor × dV) / ∫(j_tor × dV)
"""
function _compute_psipla(eqt)
    eqt1d = eqt.profiles_1d
    psi = eqt1d.psi
    j_tor = eqt1d.j_tor
    volume = eqt1d.volume
    nv = length(volume)
    if nv < 2
        return NaN
    end
    dV = diff(volume)
    psi_mid = 0.5 .* (psi[1:nv-1] .+ psi[2:nv])
    j_mid = 0.5 .* (j_tor[1:nv-1] .+ j_tor[2:nv])
    numerator = sum(psi_mid .* j_mid .* dV)
    denominator = sum(j_mid .* dV)
    if abs(denominator) < 1e-30
        return NaN
    end
    return numerator / denominator
end

"""
    _compute_co2_density(dd::IMAS.dd)

Compute CO2 interferometer line-averaged electron density for each channel.
Uses the same IMAS.line_average as ActorInterferometer.
Returns vector of line-averaged ne [m⁻³] for each interferometer channel.
"""
function _compute_co2_density(dd::IMAS.dd)
    intf = dd.interferometer
    if isempty(intf.channel)
        return Float64[]
    end
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    dens_co2 = Float64[]
    for ch in intf.channel
        try
            result = IMAS.line_average(eqt, cp1d.electrons.density_thermal, cp1d.grid.rho_tor_norm, ch.line_of_sight)
            push!(dens_co2, result.line_average)
        catch
            push!(dens_co2, NaN)
        end
    end
    return dens_co2
end
