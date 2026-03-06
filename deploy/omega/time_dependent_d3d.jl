using FUSE

for reconstruction in (false, true)
    sty = FUSE.study_parameters(:Postdictive)
    sty.device = :D3D
    sty.shots = [168830]
    sty.reconstruction = reconstruction
    sty.save_folder = mktempdir()
    sty.server = "localhost"
    sty.n_workers = 0
    FUSE.run(FUSE.StudyPostdictive(sty))
end
