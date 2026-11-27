#= ========= =#
#  ActorGPEC  #
#= ========= =#
import GeneralizedPerturbedEquilibrium as GPEC

function gpec_step!(actor::ActorGPEC)
    dd = actor.dd
    par = actor.par

    eq_config = GPEC.Equilibrium.EquilibriumConfig(; eq_type="imas", jac_type=par.jac_type, psilow=par.psilow,
        psihigh=par.psihigh, mpsi=par.mpsi, mtheta=par.mtheta, imas_cocos=par.imas_cocos)
    equil = GPEC.Equilibrium.setup_equilibrium(eq_config, dd)

    wall = par.wall_shape == "nowall" ? GPEC.Vacuum.WallShapeSettings(; shape="nowall") :
           GPEC.Vacuum.WallShapeSettings(; shape=par.wall_shape, a=par.wall_a)
    problem = GPEC.EulerLagrangeProblem(equil; nn=par.nn_low:par.nn_high, wall, vac_flag=par.vac_flag,
        local_stability_flag=par.local_stability_flag, psiedge=par.psiedge, qlow=par.qlow,
        write_outputs_to_HDF5=false, verbose=par.verbose)
    ffs = GPEC.solve(problem, GPEC.Riccati())

    GPEC.write_imas(dd, (; ffs))

    return actor
end
