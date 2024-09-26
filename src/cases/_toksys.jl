function from_TokSys(dd::IMAS.dd, mat::Dict)    
    
    # ==========
    
    eqs = mat["eqs"]

    nr = Int(eqs["nr"][1,1])
    nz = Int(eqs["nz"][1,1])

    time = Float64.(eqs["time"][1,:])
    ip = Float64.(eqs["cpasma"][1,:])
    psimag = Float64.(eqs["psimag"][1,:])
    psibry = Float64.(eqs["psibry"][1,:])

    b0 = Float64.(eqs["bzero"][1,:])
    r0 = Float64.(eqs["rzero"][1,:][1])
    
    rmaxis = Float64.(eqs["rmaxis"][1,:])
    zmaxis = Float64.(eqs["zmaxis"][1,:])    

    ffprim = Float64.(vcat(eqs["ffprim"][1, :]...)')
    pprime = Float64.(vcat(eqs["pprime"][1, :]...)')
    pres = Float64.(vcat(eqs["pres"][1, :]...)')
    fpol = Float64.(vcat(eqs["fpol"][1, :]...)')
    q = Float64.(vcat(eqs["q"][1, :]...)')
    psi_norm = range(0.0, 1.1, length(ffprim[:,1]))

    rbbbs = Float64.(hcat(eqs["rbbbs"][1, :]...))
    zbbbs = Float64.(hcat(eqs["zbbbs"][1, :]...))

    rlim = Float64.(hcat(eqs["Rlim"][1, :]...)[:,1])
    zlim = Float64.(hcat(eqs["Zlim"][1, :]...)[:,1])
    rg = Float64.(eqs["rg"][1, :][1][1,:])
    zg = Float64.(eqs["zg"][1, :][1][:,1])

    psirz = [Float64.(x) for x in eqs["psizr"][1, :]]

    # ===========
    
    pf_time = Float64.(mat["Vct"])[:,1]
    pf_voltages = Float64.(mat["Vcd"])
    pf_names = mat["Tokamak"]["c"]["ccnames"]

    # ===========

    dd.equilibrium.time = time
    dd.equilibrium.vacuum_toroidal_field.r0 = r0
    dd.equilibrium.vacuum_toroidal_field.b0 = b0
    resize!(dd.equilibrium.time_slice, length(time))
    for (k,eqt) in enumerate(dd.equilibrium.time_slice)
        eq1d = eqt.profiles_1d
        eq2d = resize!(eqt.profiles_2d, 1)[1]
        eq2d.grid_type.index = 1

        eqt.global_quantities.magnetic_axis.r = rmaxis[k]
        eqt.global_quantities.magnetic_axis.z = zmaxis[k]

        #eqt.boundary.geometric_axis.r = g.rcentr
        #eqt.boundary.geometric_axis.z = g.zmid

        eqt.global_quantities.ip = ip[k]

        eq1d.psi = range(psimag[k], psibry[k], length(psi_norm))

        eq1d.q = q[:,k]
        eq1d.pressure = pres[:,k]
        eq1d.dpressure_dpsi = pprime[:,k]
        eq1d.f = fpol[:,k]
        eq1d.f_df_dpsi = ffprim[:,k]

        eq2d.grid.dim1 = rg
        eq2d.grid.dim2 = zg
        eq2d.psi = collect(psirz[k]')
    end
    
    resize!(dd.wall.description_2d, 1)
    resize!(dd.wall.description_2d[1].limiter.unit, 1)
    dd.wall.description_2d[1].limiter.unit[1].outline.r = rlim
    dd.wall.description_2d[1].limiter.unit[1].outline.z = zlim
    
    dd.pulse_schedule.profiles_control.time = time
    dd.pulse_schedule.profiles_control.psi_norm = psi_norm
    dd.pulse_schedule.profiles_control.dpressure_dpsi.reference = pprime
    dd.pulse_schedule.profiles_control.f_df_dpsi.reference = ffprim
    
    dd.pulse_schedule.flux_control.time = time
    dd.pulse_schedule.flux_control.i_plasma.reference = ip
    
    dd.pulse_schedule.position_control.time = time
    dd.pulse_schedule.position_control.magnetic_axis.r.reference = rmaxis
    dd.pulse_schedule.position_control.magnetic_axis.z.reference = zmaxis
    resize!(dd.pulse_schedule.position_control.boundary_outline, size(rbbbs)[1])
    for k in 1:size(rbbbs)[1]
        dd.pulse_schedule.position_control.boundary_outline[k].r.reference = rbbbs[k,:]
        dd.pulse_schedule.position_control.boundary_outline[k].z.reference = zbbbs[k,:]
    end
    
    dd.pulse_schedule.pf_active.time = pf_time
    resize!(dd.pulse_schedule.pf_active.coil, size(pf_voltages)[2])
    for (k,voltage) in enumerate(eachcol(pf_voltages))
        dd.pulse_schedule.pf_active.coil[k].identifier = pf_names[k]
        dd.pulse_schedule.pf_active.coil[k].voltage.reference = voltage
    end
    
    IMAS.flux_surfaces(dd.equilibrium, IMAS.first_wall(dd.wall)...)
end

"""
AR => Integral(dA/R) over cell
ARcell => Integral[dA/R] inside limiter
Acell => Area inside limiter
Ag => Area of a grid cell [m2]
Pcc => coil currents = Pcc * circuit currents
Pcci => circuit currents = Pcci * coil currents
Pvc => vessel currents = Pvc * vessel circuit currents
Pvci => vessel circuit currents = Pvci * vessel currents
Pvg => the output ivg = Pvg * vessel circuit currents
RA => Integral(R*dA) over cell
RAcell => Integral[R*dA] inside limiter
Rhat => Generalized resistance in evolution equation
Rlim => R of limiter for points at corners
Vhat => Generalized voltage in evolution equation
Zlim => Z of limiter for points at corners
amax => Radius for plasma with no closed flux surface
amin => Current exists outside last-closed-flux-surface if its aminor<amin
bpdata => position and orientation of magnetic probes
bpnames => names of magnetic probes
bzero => Bphi at reference point inside the machine
ccnames => names of E,F coils
concavel => Limiter corner points into machine
constraints => Reduces nsp to 1 and nsf to 2 if > 0
defaultf => True if constrained ffprim/mu0 is the default
defaultp => True if constrained pprime is the default = ones(nr,1)
dl => sqrt(drl.^2+dzl.^2)
dpsibarlimit => Maximum change before updating response
dr => Distance radially between grid points [m]
drl => diff(rl)
dtmin => smallest time step between updates
dz => Distance in Z-direction between grid points [m]
dzl => diff(zl)
evolve_option => Option for how to evolve equilibrium
faeq => Settings for fixed-area-equilibrium
fb => fpol at boundary = fb + vf*sf
fc => Settings for initial fixed circular plasma
fcdata => F coil geometry
fcnames => names of F coils
fldata => flux loop positions
flnames => names of flux loops
fnewmin => Old calculations can be used if fnew>0 (max is 1), 0 < fnewmin < 1
gammamax => Upper limit on gamma in gsevolve
gbc => Greens between magnetic probes and coils
gbv => Greens between magnetic probes and vessel
gpb => Greens between magnetic probes and grid
grdc => Br at diagnostic points from 1 A of coil circuit currents
grdv => Br at diagnostic points from 1 A of vessel circuit currents
grpd => Br at diagnostic points from 1 A of current in grid cells
gzdc => Bz at diagnostic points from 1 A of coil circuit currents
gzdv => Bz at diagnostic points from 1 A of vessel circuit currents
gzpd => Bz at diagnostic points from 1 A of current in grid cells
hr => Settings for high-resolution profiles
i16 => index differences to 16 grid points
iE => Flag nonzero elements of E
ics => Internal control settings for the GSevolve module (gsevolve_sfun.m)
il => Indices used for interpolation to rl, zl
irl => r-indices
iu => Contains indices to ic, ih, sp, sf in actuator vector u
ix => Contains indices to ic, iv, ih, sp, sf in state vector x
iy => Contains indices for diagnostics in output vector y
izl => z-indices
kl => Collapsed indices
limpoly => Variables used to determine if point is inside limiter
lstarx => Inductance in cables between power supplies and coils
mcc => Mutuals between coil circuits and coil circuits
mcv => Mutuals between coil circuits and vessel circuits
mdc => Mutual inductance between dp and coil circuits
mdv => Mutual inductance between dp and vessel circuits
mf0 => ffprim/mu0 coefficients f0 = mf0*sf
mf1 => ffprim/mu0 coefficients f1 = mf1*sf
mf2 => ffprim/mu0 coefficients f2 = mf2*sf
mf3 => ffprim/mu0 coefficients f3 = mf3*sf
mgg => Mutuals between grid points and grid cells
mhc => Mutuals between voltage loops and coils
mhv => Mutuals between voltage loops and vessel
mlc => Mutuals between flux loops and coils
mlv => Mutuals between flux loops and vessel
mp0 => pprime coefficients p0 = mp0*sp
mp1 => pprime coefficients p1 = mp1*sp
mp2 => pprime coefficients p2 = mp2*sp
mp3 => pprime coefficients p3 = mp3*sp
mpc => Mutuals between grid points and coil circuits
mpd => Mutual inductance between grid and dp
mph => Mutuals between voltage loops and grid
mpl => Mutuals between flux loops and grid
mpp => Mutuals between grid points and lowest row of grid cells
mps => Mutuals between saddle loops and grid
mpu => Mutuals from grid to usig
mpv => Mutuals between grid points and vessel circuits
msc => Mutuals between saddle loops and coils
msv => Mutuals between saddle loops and vessel
muc => Mutuals from coil circuits to usig
muv => Mutuals from vessel to usig
mvv => Mutuals between vessel circuits and vessel circuits
nb => Number of different points on a contour
nbdef => Rectangles that can contain x-points that are inside the plasma
nbdtest => Number of points checked for jump in boundary-defining point
nbmax => Fixed number of points returned by gscontour22 (padding is nans)
nbp => Number of magnetic probes
ndp => Number of user-specified diagnostic points (dp)
nfl => Number of flux loops
nfla => Number of angles from axis to fl points
nflc => Number of contours with fl points
nfltest => Number of fluxes checked for validity of linear approximation
ngg => Number of grid points = nz*nr
nic => Number of circuits with coils (and a power supply)
niv => Number of circuits with vessel elements (no power supply)
nkf => length(psikf)-1
nkp => length(psikp)-1
nl => length(rl)
nlv => Number of voltage loops
nn => Max number of previous plasma responses remembered by gsupdate
no_j_crit => jphi<jerr or [sf;sp]<[sferr;sperr] is considered no jphi
now => Time stamp used as identity for this configuration
np => Number of contours
nr => Number of grid points radially
nrog => Number of Rogowski loops
nsavebd => number of times to save breakdown information x,c to file gscreate_inputs_*
nsf => Number of states for ffprim/mu0 and boundary fpol
nsl => Number of saddle loops
nsp => Number of states for pprime and boundary pressure
nu => Number of actuators
nusig => Number of user-defined signals
nx => Total number of states = nic+niv+nsp+nsf
nxpoints => Number of x-points in outputs rx, zx
ny => Number of outputs
nz => Number of grid points in Z-direction
pb => pres at boundary = pb + vp*sp
plot_profiles => Flag to plot profiles at end of gsprofiles
plots => Settings for plotting such as rtplot
profile => Settings for outputs named profile.*
psibarp => Normalized poloidal flux for contours
psikf => ffprim/mu0 spline knot positions
psikp => pprime spline knot positions
rg => R of grid points [m]
rgg => R for all nz*nr grid points
rl => R of limiter at corners and lines of a doubly dense grid
rlc => Rogowski signals from coils = rlc*ic
rldata => Fraction of coil, vessel, plasma inside Rogowski loops
rlnames => names of rogowski coils
rlv => Rogowski signals from vessel = rlv*iv
rpl => Rogowski signals from plasma = rpl'*pcurrt(:)
rzero => R for a reference point inside the machine
sf0 => spline parameters that create default ffprim
sg0 => spline parameters that create default ffprim peaking
sldata => saddle loop geometry
slnames => names of saddle loops
sp0 => spline parameters that create default pprime
term => Settings for plasma termination
tokamak => Name of tokamak
trl => rl(1:nl-1)'-irl
turnin => Matrix to vector along limiter into machine
tzl => zl(1:nl-1)'-izl
vbus => R,L,M for voltage output vbus = R*x + d/dt(L*x+M*pcurrt(:))
vf => fpol at boundary = fb + vf*sf
vp => pres at boundary = pb + vp*sp
vvdata => Vessel element geometry
wl => Flux at rl,zl = sum(wl.*psizr(il),2)
wld1 => wld1(i,:) calculates derivative at i along vector from i to i+1
wld2 => wld2(i,:) calculates derivative at i+1 along vector from i to i+1
wlr => d(Flux)/d(R) at rl,zl = sum(wlr.*psizr(il),2)
wlz => d(Flux)/d(Z) at rl,zl = sum(wlz.*psizr(il),2)
xpointorder => 1,3=sorted according to normalized flux, 2,4=sorted according to zx
zg => Z of grid points [m]
zgg => Z for all nz*nr grid points
zl => Z of limiter at corners and lines of a doubly dense grid
zzero => Z for a reference point inside the machine
"""