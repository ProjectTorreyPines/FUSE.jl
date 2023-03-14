using Aqua, FUSE
using Pkg

pkg_specs = Pkg.Operations.load_all_deps(Pkg.API.EnvCache())
to_test = pkg_specs[findall(x->!isnothing(x.path), pkg_specs)]

pkg_ids = map(p -> Base.PkgId(p.uuid, p.name), to_test)
stale_deps = Aqua.analyze_stale_deps(pkg_ids)
# ps = to_test .=> stale_deps
map(x->x.name, pkg_ids) .=> stale_deps