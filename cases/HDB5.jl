function Parameters(::Type{Val{:HDB5}}; case)
    par = Parameters()
    par.general.casename = "HDB_$(case[:TOK])_$(case[:SHOT]))"
    par.general.init_from = :scalars

    # Equilibrium parameters
    par.equilibrium.B0 = abs(case[:BT])
    par.equilibrium.R0 = case[:RGEO]
    par.equilibrium.Z0 = 0.0
    par.equilibrium.ϵ = case[:AMIN]/case[:RGEO]
    par.equilibrium.κ = case[:KAPPA]
    par.equilibrium.δ = case[:DELTA] 
    par.equilibrium.βn = 1.0
    par.equilibrium.ip = abs(case[:IP])
    par.equilibrium.x_point = false
    par.equilibrium.symmetric = true
    
    # Core_profiles parameters
    par.core_profiles.ne_ped = case[:NEL] / 1.3
    par.core_profiles.n_peaking = 1.5
    par.core_profiles.T_shaping = 1.8
    par.core_profiles.w_ped = 0.03
    par.core_profiles.zeff = case[:ZEFF]
    par.core_profiles.rot_core = 50e3
    par.core_profiles.ngrid = 201
    par.core_profiles.bulk = :D
    par.core_profiles.impurity = :C
    
    # nbi
    if case[:PNBI] + case[:POHM] > 0.
        par.nbi.beam_power = case[:PNBI] + case[:POHM]
        par.nbi.beam_energy = case[:ENBI]
        par.nbi.beam_mass = 2
        par.nbi.toroidal_angle = 0.
    end
        
    if case[:PECRH] > 0
        par.ec.power_launched = case[:PECRH]
    end
       
    if case[:PICRH] > 0
        par.ic.power_launched = case[:PICRH]
    end
    
    return par
end