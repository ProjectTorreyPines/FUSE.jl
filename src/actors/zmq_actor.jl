#= ========== =#
#  ActorZMQ    #
#= ========== =#
# Coupling actor for FUSE <-> GSLite/GSEvolve communication via ZeroMQ (REQ/REP pattern)
#
# In the dynamic plasma loop:
#   Phase 1 start: receive!(actor_zmq) — gets flux matrix + Ip from GSLite
#   Phase 2 end:   send!(actor_zmq)    — sends FRESCO equilibrium back to GSLite
#
# Message format: Protocol Buffers over ZMQ REQ/REP
# Schema defined in zmq_messages.proto

import Interpolations
import ZMQ
import ProtoBuf

include(joinpath(@__DIR__, "zmq_proto_generated", "zmq_messages_pb.jl"))
using .zmq_messages_pb: FUSERequest, WireDataForFUSE, WireDataFromFUSE, Ack

# Wire-contract version for FUSE↔GSLite ZMQ coupling.
# Bump in lockstep with the `schema_version` fields in zmq_messages.proto on any
# breaking schema change. FUSE sends this in every FUSERequest and refuses any
# WireDataForFUSE whose schema_version does not match.
const SCHEMA_VERSION = Int32(1)

# --- Protobuf over ZMQ helpers ---

function _pb_send(socket::ZMQ.Socket, msg)
    io = IOBuffer()
    e = ProtoBuf.ProtoEncoder(io)
    ProtoBuf.encode(e, msg)
    ZMQ.send(socket, take!(io))
end

function _pb_recv(socket::ZMQ.Socket, ::Type{T}) where {T}
    raw = ZMQ.recv(socket)
    io = IOBuffer(raw)
    d = ProtoBuf.ProtoDecoder(io)
    return ProtoBuf.decode(d, T)
end

@actor_parameters_struct ActorZMQ{T} begin
    endpoint::Entry{String} = Entry{String}("-", "ZMQ endpoint (e.g., tcp://localhost:5555)"; default="tcp://localhost:5555")
    timeout_ms::Entry{Int} = Entry{Int}("ms", "ZMQ receive timeout"; default=10000)
    enabled::Entry{Bool} = Entry{Bool}("-", "Enable ZMQ coupling with external code"; default=false)
    exit_on_timeout::Entry{Bool} = Entry{Bool}("-", "Call exit(1) after tearing the socket down if a ZMQ exchange fails (e.g. GSLite goes silent past timeout_ms)"; default=false)
end

mutable struct ActorZMQ{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorZMQ{P}}
    context::Union{Nothing,ZMQ.Context}
    socket::Union{Nothing,ZMQ.Socket}
    is_connected::Bool
    had_psizr::Bool
    # Previous-step values for time derivatives in send!
    prev_time::Float64
    prev_betap::Float64
    prev_li::Float64
    prev_psipla::Float64
    prev_co2::Vector{Float64}  # last valid CO2 density values (fallback if computation fails)
end

function ActorZMQ(dd::IMAS.dd{D}, par::FUSEparameters__ActorZMQ{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorZMQ)
    par = OverrideParameters(par; kw...)
    return ActorZMQ{D,P}(dd, par, nothing, nothing, false, false, NaN, NaN, NaN, NaN, Float64[])
end

"""
    ActorZMQ(dd::IMAS.dd, act::ParametersAllActors; kw...)

Coupling actor for exchanging data with external codes (e.g., GSLite/GSEvolve) via ZeroMQ.

This actor uses a REQ/REP pattern where FUSE initiates each exchange:
- `receive!`: FUSE sends "ready" → external code replies with WireDataForFUSE
- `send!`: FUSE sends WireDataFromFUSE → external code replies with "ack"

Intended integration in ActorDynamicPlasma:
- Phase 1 start: `receive!(actor_zmq)` — receive ψ(R,Z), Ip, NBI, gas, coil currents
- Phase 2 end: `send!(actor_zmq)` — send betap, li, p_res, derivatives, CO2 density

A timeout or protocol error during a ZMQ exchange terminates the coupled run; there is no automatic reconnect.

Each exchange begins with a wire-contract version check: FUSE sends `SCHEMA_VERSION`
in `FUSERequest`, and any mismatch with `WireDataForFUSE.schema_version` terminates the
run with a loud error (no negotiation). GSLite must reply with `Ack.ok = true` on
accepted steps; `ok = false` (e.g. NaN input, solver divergence) likewise terminates
the run, surfacing `Ack.error`.
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
    ZMQ.setsockopt(sock, ZMQ.RCVTIMEO, actor.par.timeout_ms)
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

WireDataForFUSE fields (matching C++ struct):
- `sim_time`:        double          — GSLite simulation time [s]
- `Ip_latest`:       double          — Latest Ip measurement [A] → pulse_schedule
- `Ip_avg`:          double          — Average Ip measurement [A] → dd._aux (for NN ne predictor)
- `Bt`:              double          — Vacuum toroidal field at R0 [T] (first message only)
- `pr15v`:           double          — PR15V Ip measurement [A] → dd._aux (for NN ne predictor)
- `I_coil`:          double[24]      — PCS coil currents [A] → dd._aux[:zmq_I_coil] (for NN pedestal predictor).
                                       Order matches DIII-D `PCSpcsRtnetCoilNames`:
                                         1:PCECOILA  2:PCE89DN   3:PCE567UP  4:PCECOILB
                                         5:PCE89UP   6:PCE567DN  7:PCF1A..15:PCF9A
                                        16:PCF1B..24:PCF9B
                                       PedestalPredictor FPE consumes ecoila (idx 1), ecoilb (idx 4),
                                       and f1a..f9b (idx 7..24); the C-coil entries (idx 2,3,5,6) are unused.
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

    # REQ: send ready signal, then receive protobuf data.
    # Fail-fast: any timeout / EFSM / decode error terminates the coupled run.
    msg = try
        _pb_send(actor.socket, FUSERequest("ready", dd.global_time, SCHEMA_VERSION))
        _pb_recv(actor.socket, WireDataForFUSE)
    catch e
        @error "ActorZMQ.receive!: exchange with GSLite failed at t=$(dd.global_time) s — terminating coupled run" exception=(e, catch_backtrace())
        disconnect!(actor)
        if actor.par.exit_on_timeout
            @error "ActorZMQ: exit_on_timeout=true — calling exit(1) to terminate the FUSE process"
            exit(1)
        end
        rethrow()
    end

    # --- Wire-contract version check: any mismatch is fatal, before any other handling ---
    if msg.schema_version != SCHEMA_VERSION
        error("ActorZMQ: ZMQ schema mismatch — FUSE sent version $SCHEMA_VERSION, GSLite replied with version $(msg.schema_version). Incompatible; rebuild the older side against the matching zmq_messages.proto.")
    end

    @info "ActorZMQ: received data at sim_time=$(msg.sim_time) s"

    # --- Sync FUSE clock to GSLite ---
    dd.global_time = msg.sim_time

    # --- Check for end-of-simulation signal from GSLite ---
    if msg.done
        @info "ActorZMQ: GSLite signaled end of simulation"
        disconnect!(actor)
        return actor
    end

    # --- Update Ip in pulse_schedule so QED/FRESCO pick it up ---
    if msg.has_Ip_latest
        ps_fc = dd.pulse_schedule.flux_control
        IMAS.set_time_array(ps_fc.i_plasma, :reference, dd.global_time, msg.Ip_latest)
    end

    # --- Store auxiliary signals for NN ne predictor ---
    # Uses dd._aux (same pattern as FUSE workflow/logging metadata)
    # Stored as (times=Float64[], values=...) parallel vectors to avoid Float64 dict keys
    aux = getfield(dd, :_aux)
    if msg.has_Ip_avg
        if :zmq_Ip_avg ∉ keys(aux)
            aux[:zmq_Ip_avg] = (times=Float64[], values=Float64[])
        end
        push!(aux[:zmq_Ip_avg].times, dd.global_time)
        push!(aux[:zmq_Ip_avg].values, msg.Ip_avg)
    end
    if msg.has_pr15v
        if :zmq_pr15v ∉ keys(aux)
            aux[:zmq_pr15v] = (times=Float64[], values=Float64[])
        end
        push!(aux[:zmq_pr15v].times, dd.global_time)
        push!(aux[:zmq_pr15v].values, msg.pr15v)
    end

    # --- Update Bt (first step only, constant per shot) ---
    if msg.has_Bt
        b0 = msg.Bt
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
    if !isempty(msg.I_coil)
        if :zmq_I_coil ∉ keys(aux)
            aux[:zmq_I_coil] = (times=Float64[], values=Vector{Float64}[])
        end
        push!(aux[:zmq_I_coil].times, dd.global_time)
        push!(aux[:zmq_I_coil].values, msg.I_coil)
    end

    # --- Update NBI power per beam [W] ---
    # Uses IMAS.set_time_array to properly build up time history,
    # which is needed by smooth_beam_power (smooths over τ_thermalization)
    if !isempty(msg.pinj_per_beam)
        ps_nbi = dd.pulse_schedule.nbi
        n_beams = min(length(msg.pinj_per_beam), length(ps_nbi.unit))
        time0 = dd.global_time
        for k in 1:n_beams
            IMAS.set_time_array(ps_nbi.unit[k].power, :reference, time0, msg.pinj_per_beam[k])
        end
        @info "ActorZMQ: updated NBI power for $n_beams beams"
    end

    # --- Update NBI acceleration voltage per beam [eV] ---
    if !isempty(msg.nbi_acc_voltage)
        ps_nbi = dd.pulse_schedule.nbi
        n_beams = min(length(msg.nbi_acc_voltage), length(ps_nbi.unit))
        time0 = dd.global_time
        for k in 1:n_beams
            IMAS.set_time_array(ps_nbi.unit[k].energy, :reference, time0, msg.nbi_acc_voltage[k])
        end
        @info "ActorZMQ: updated NBI voltage for $n_beams beams"
    end

    # --- Store gas calibration values for NN ne predictor ---
    # gas_cal[NGAS]: index 1=gasA, 2=gasB, 3=gasC, 4=gasD, 5=gasE
    if !isempty(msg.gas_cal)
        gas_names = (:zmq_gasa_cal, :zmq_gasb_cal, :zmq_gasc_cal, :zmq_gasd_cal, :zmq_gase_cal)
        for k in 1:min(length(msg.gas_cal), length(gas_names))
            sym = gas_names[k]
            if sym ∉ keys(aux)
                aux[sym] = (times=Float64[], values=Float64[])
            end
            push!(aux[sym].times, dd.global_time)
            push!(aux[sym].values, msg.gas_cal[k])
        end
    end

    # --- Update equilibrium from flat psizr array ---
    if !isempty(msg.psizr)
        psizr_flat = msg.psizr

        eqt = dd.equilibrium.time_slice[]

        # Ensure profiles_2d exists
        if isempty(eqt.profiles_2d)
            resize!(eqt.profiles_2d, 1)
        end
        p2d = eqt.profiles_2d[1]

        # Resolve grid: prefer the wire grid; fall back to the cached dd grid.
        if !isempty(msg.r_grid) && !isempty(msg.z_grid)
            p2d.grid.dim1 = msg.r_grid
            p2d.grid.dim2 = msg.z_grid
        elseif isempty(p2d.grid.dim1)
            @warn "ActorZMQ: psizr received but no R/Z grid available; skipping equilibrium update"
            return actor
        end
        nR = length(p2d.grid.dim1)
        nZ = length(p2d.grid.dim2)
        if length(psizr_flat) != nR * nZ
            error("ActorZMQ: psizr length $(length(psizr_flat)) != nR*nZ = $(nR*nZ) — check GSLite vs FUSE grid agreement")
        end
        psi_rz = reshape(psizr_flat, nR, nZ)  # GSLite stores column-major (R varies fastest)
        p2d.psi = psi_rz
        p2d.grid_type.index = 1  # rectangular grid

        rgrid = range(p2d.grid.dim1[1], p2d.grid.dim1[end], length=length(p2d.grid.dim1))
        zgrid = range(p2d.grid.dim2[1], p2d.grid.dim2[end], length=length(p2d.grid.dim2))
        fw_r, fw_z = IMAS.first_wall(dd.wall)

        if !actor.had_psizr
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
        actor.had_psizr = true
    end

    # --- Compute and store ohmic & NBI power in dd._aux for NN predictor ---
    try
        Pohm = IMAS.total_power_source(IMAS.ohmic_source!(dd).profiles_1d[])
        if :zmq_Pohm ∉ keys(aux)
            aux[:zmq_Pohm] = (times=Float64[], values=Float64[])
        end
        push!(aux[:zmq_Pohm].times, dd.global_time)
        push!(aux[:zmq_Pohm].values, Pohm)
    catch e
        @warn "ActorZMQ: failed to compute Pohm" exception=e
    end

    try
        Pnbi = @ddtime(dd.summary.heating_current_drive.power_launched_nbi.value)
        if :zmq_Pnbi ∉ keys(aux)
            aux[:zmq_Pnbi] = (times=Float64[], values=Float64[])
        end
        push!(aux[:zmq_Pnbi].times, dd.global_time)
        push!(aux[:zmq_Pnbi].values, Pnbi)
    catch e
        @warn "ActorZMQ: failed to compute Pnbi" exception=e
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

WireDataFromFUSE fields (matching C++ struct):
- `sim_time`:     double        — Current simulation time [s]
- `valid`:        bool          — false on first send (derivatives unreliable), true thereafter
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
    psipla_now = _compute_psipla(eqt)
    dt = time_now - actor.prev_time
    if isnan(actor.prev_time) || dt <= 0.0
        betap_dot = 0.0
        li_dot = 0.0
        p_res = 0.0
    else
        betap_dot = (betap - actor.prev_betap) / dt
        li_dot = (li - actor.prev_li) / dt
        # Plasma resistance: Rp = dψ_plasma/dt / Ip, clamped >= 1e-9
        Ip_val = eqt.global_quantities.ip
        p_res = max((psipla_now - actor.prev_psipla) / dt / Ip_val, 1e-9)
    end

    msg = WireDataFromFUSE(
        time_now,
        !isnan(actor.prev_time),  # valid: false on first send (no prior data for derivatives)
        betap,
        betap_dot,
        li,
        li_dot,
        p_res,
        _compute_co2_density(actor)
    )

    # REQ: send data, then receive acknowledgment.
    # Fail-fast: any timeout / EFSM / decode error terminates the coupled run.
    ack = try
        _pb_send(actor.socket, msg)
        _pb_recv(actor.socket, Ack)
    catch e
        @error "ActorZMQ.send!: exchange with GSLite failed at t=$(time_now) s — terminating coupled run" exception=(e, catch_backtrace())
        disconnect!(actor)
        if actor.par.exit_on_timeout
            @error "ActorZMQ: exit_on_timeout=true — calling exit(1) to terminate the FUSE process"
            exit(1)
        end
        rethrow()
    end

    # Application-level rejection: GSLite parsed the message but is refusing it
    # (NaN input, solver divergence, etc.). Surface its diagnostic and bail.
    if !ack.ok
        error("ActorZMQ.send!: GSLite rejected WireDataFromFUSE at t=$(time_now) s — status=$(ack.status), ok=$(ack.ok), error=$(ack.error)")
    end
    @info "ActorZMQ: sent betap=$betap, li=$li, p_res=$p_res at t=$(time_now) s"

    # Store current values for next step's derivatives
    actor.prev_time = time_now
    actor.prev_betap = betap
    actor.prev_li = li
    actor.prev_psipla = psipla_now

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
    _compute_co2_density(actor::ActorZMQ)

Compute CO2 interferometer line-integrated electron density for each channel.
Uses the same IMAS.line_average as ActorInterferometer.
On failure, logs a warning and falls back to last valid values.
Returns vector of line-integrated ne [m⁻²] for each interferometer channel.
"""
function _compute_co2_density(actor::ActorZMQ)
    dd = actor.dd
    intf = dd.interferometer
    if isempty(intf.channel)
        return Float64[]
    end
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    dens_co2 = Float64[]
    for (i, ch) in enumerate(intf.channel)
        try
            result = IMAS.line_average(eqt, cp1d.electrons.density_thermal, cp1d.grid.rho_tor_norm, ch.line_of_sight)
            push!(dens_co2, result.line_integral * 1e-4)
        catch e
            fallback = (i <= length(actor.prev_co2)) ? actor.prev_co2[i] : NaN
            @warn "ActorZMQ: CO2 channel $i failed, using fallback=$fallback" exception=e
            push!(dens_co2, fallback)
        end
    end
    actor.prev_co2 = dens_co2
    return dens_co2
end
