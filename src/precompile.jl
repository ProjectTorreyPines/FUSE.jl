if VERSION > v"1.8"

    @precompile_setup begin
        @precompile_all_calls begin
            FUSE.warmup()
        end
    end

end