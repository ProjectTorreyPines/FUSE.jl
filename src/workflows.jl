function transport_workflow(dd::IMAS.dd, par::Parameters; do_plot = false, verbose = false, kw...)
    if do_plot
        plot(dd.core_profiles; color = :gray, label = "")
    end
    actor = TaueNNactor(dd; kw...)
    step(actor; verbose = verbose)
    finalize(actor)
    if do_plot
        display(plot!(dd.core_profiles))
    end
    return dd
end