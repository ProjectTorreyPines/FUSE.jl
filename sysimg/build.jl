using PackageCompiler

create_sysimage(["FUSE"],
                sysimage_path=joinpath(@__DIR__, "sysimg.so"),
                project=joinpath(@__DIR__, ".."),
                precompile_execution_file="precompile_script.jl",
                incremental=true, # set this to false to generate a smaller sysimg
                filter_stdlibs=false # set this to true to generate a smaller sysimg
                #uncomment this line to anonymize the sysimg for proprietary distributions
                #, sysimage_build_args=`-g0 --strip-ir --strip-metadata --compile=all`
                )
