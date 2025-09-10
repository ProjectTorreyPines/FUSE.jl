using FUSE
using Test
import FUSE.CSV as CSV
import FUSE.HDF5 as HDF5
using FUSE.DataFrames

@testset "study_database" begin
    
    # =============================
    # Create a fake db for testing
    # =============================

    # --- Create three dd objects from sample ---
    dd1 = FUSE.IMAS.json2imas(joinpath(@__DIR__, "..", "sample", "CAT_eq_ods.json"))
    dd2 = deepcopy(dd1)
    dd2.equilibrium.vacuum_toroidal_field.b0 = [0.5]
    dd3 = deepcopy(dd1)
    dd3.equilibrium.vacuum_toroidal_field.b0 = [0.8]

    # --- Create study_database_item (absent fields -> nothing) ---
    item1 = FUSE.study_database_item(dd = dd1)
    item2 = FUSE.study_database_item(dd = dd2, ini = nothing, act = nothing, log = nothing, timer = nothing, error = nothing)
    item3 = FUSE.study_database_item(dd = dd3, ini = nothing, act = nothing, log = "this is log", timer = nothing, error = "this is error")
    items = [item1, item2, item3]

    # --- Simple DataFrame metadata (only columns used in the tests) ---
    df = FUSE.extract_dds_to_dataframe([dd1, dd2, dd3])
    df_additional = DataFrame(
        dir = ["case1", "case2", "case3"],
        case = [1, 2, 3],
        gparent = ["/case1", "/case2", "/case3"],  # must be valid HDF5 group paths
        status = ["success", "success", "fail"],
        worker_id = [1, 2, 3],
        elapsed_time = [0.1, 0.2, 0.3],
    )
    df = hcat(df, df_additional)   # assumes no name conflicts and same row count

    # --- Create study_database ---
    db = FUSE.study_database(df, items)

    
    # =============================
    #        Start testing
    # =============================
    
    # Basic shape checks
    @test length(db.items) == nrow(db.df)
    @test names(db.df) == names(df)

    @testset "filter" begin
        # filter(Function) -> new db
        db_success = filter(row -> row.status == "success", db)
        @test nrow(db_success.df) == 2
        @test length(db_success.items) == 2

        db_fail = filter(row -> row.status == "fail", db)
        @test nrow(db_fail.df) == 1

        # filter(mask::Vector{Bool}) -> new db
        mask_success = db.df.status .== "success"
        db_mask_success = filter(mask_success, db)
        @test nrow(db_mask_success.df) == 2

        # filter!(Function) -> in-place
        db_copy = FUSE.study_database(copy(db.df), copy(db.items))
        filter!(row -> row.case != 3, db_copy)
        @test nrow(db_copy.df) == 2
        @test all(db_copy.df.case .!= 3)

        # filter!(mask::Vector{Bool}) -> in-place by mask
        db_copy_mask = FUSE.study_database(copy(db.df), copy(db.items))
        mask = db_copy_mask.df.status .== "success"
        filter!(mask, db_copy_mask)
        @test nrow(db_copy_mask.df) == 2
        @test all(db_copy_mask.df.status .== "success")
    end

    @testset "subset" begin
        # subset (keep case==2)
        db_subset = DataFrames.subset(db, :case => DataFrames.ByRow(==(2)))
        @test nrow(db_subset.df) == 1
        @test db_subset.df.case[1] == 2

        # subset by status
        db_subset2 = DataFrames.subset(db, :status => DataFrames.ByRow(==("success")))
        @test nrow(db_subset2.df) == 2

        # subset! (keep case<=2)
        db_copy2 = FUSE.study_database(copy(db.df), copy(db.items))
        DataFrames.subset!(db_copy2, :case => DataFrames.ByRow(<=(2)))
        @test nrow(db_copy2.df) == 2
        @test all(db_copy2.df.case .<= 2)
    end

    @testset "getindex" begin
        # getindex variations (column-only)
        @test db[:, :case] == df[:, :case]
        @test db[:, [:case, :dir]] == df[:, [:case, :dir]]
        @test db[:, [:case, :dir, :status]] == df[:, [:case, :dir, :status]]

        # getindex variations (row and column)
        @test db[1, :case] == df[1, :case]
        @test db[2:3, :dir] == df[2:3, :dir]
        @test db[[1, 3], [:case, :status]] == df[[1, 3], [:case, :status]]
        @test db[db.df.status .== "success", [:case, :dir]] == df[df.status .== "success", [:case, :dir]]

        # getindex variations (row-only)
        @test isa(db[1], FUSE.study_database)
        @test nrow(db[1].df) == 1
        @test nrow(db[1:2].df) == 2
        @test nrow(db[[1, 2, 3]].df) == 3
        @test nrow(db[[true, true, true]].df) == 3
    
        # two-arg getindex delegates to DataFrame
        @test db[:, :case] == df[:, :case]
    end

    @testset "equality" begin
        # identical copy
        db_copy = FUSE.study_database(deepcopy(db.df), deepcopy(db.items))
        @test db == db_copy
        @test isequal(db, db_copy)

        # different order -> not equal
        db_reordered = FUSE.study_database(db.df[[2, 1, 3], :], db.items[[2, 1, 3]])
        @test db != db_reordered

        # modify one field in items -> not equal
        db_mod = FUSE.study_database(copy(db.df), deepcopy(db.items))
        db_mod.items[1].log = "changed"
        @test db != db_mod
    end

    @testset "File I/O" begin
        mktempdir() do save_dir
            h5file = joinpath(save_dir, "fake_db.h5")

            # Disable timer to keep parity on item.timer (remain nothing)
            FUSE.save_study_database(h5file, db; timer=false)
            @test isfile(h5file)

            # Load entire database
            db2 = FUSE.load_study_database(h5file)
            # Equality with normalization: loaded db must match saved db
            @test db2 == db
            @test isequal(db2, db)

            # Load a single group path
            db_case1 = FUSE.load_study_database(h5file, "/case1")
            @test nrow(db_case1.df) == 1

            # Load with a filtering function
            db_cond = FUSE.load_study_database(h5file, x -> x.case == 2)
            @test nrow(db_cond.df) == 1

            # Sample-and-write: by groups
            h5file_sampled = joinpath(save_dir, "fake_db_sampled.h5")
            df_sampled = FUSE.sample_and_write_study_database(h5file, h5file_sampled, ["/case2"])
            @test isfile(h5file_sampled)
            @test nrow(df_sampled) == 1

            # Sample-and-write: by ratio
            h5file_sampled2 = joinpath(save_dir, "fake_db_sampled2.h5")
            df_sampled2 = FUSE.sample_and_write_study_database(h5file, h5file_sampled2; sampling_ratio=1.0)
            @test isfile(h5file_sampled2)
            @test nrow(df_sampled2) == 3
        end
        # mktempdir do ... end ensures temporary dir cleanup automatically
    end
end
