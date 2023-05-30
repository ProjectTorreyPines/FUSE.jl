using Profile
using FUSE

@profile FUSE.warmup()

Profile.clear()

@profile FUSE.warmup()

using PProf

pprof()