using FUSE
using Test

@testset "nb_toroidal_velocity_fraction" begin
    # A straight beam keeps its angular momentum about the machine axis, so the lever arm R * ftor is the same
    # at every point of the path: tangency radius times cos(vertical tilt), wherever the source sits toroidally.
    source_r = 8.0
    tt = collect(range(0.0, 14.0, 201))
    for source_phi in deg2rad.([0.0, 30.0, 150.0, 210.0, 330.0]), tangency_radius in (0.76, 1.15), tilt in (0.0, 0.28), direction in (1, -1)
        # same construction as ActorSimpleNB
        angleh = direction * asin(tangency_radius / source_r)
        vx0 = -cos(angleh) * cos(tilt)
        vy0 = -sin(angleh) * cos(tilt)
        vz = cos(angleh) * sin(tilt)
        px = source_r * cos(source_phi)
        py = source_r * sin(source_phi)
        vx = vx0 * cos(source_phi) - vy0 * sin(source_phi)
        vy = vx0 * sin(source_phi) + vy0 * cos(source_phi)

        Xs, Ys, Zs, Rs = IMAS.pencil_beam([px, py, 0.0], [vx, vy, vz], tt)
        @test any(<(0.0), Xs) # the path visits X < 0, where asin(Y / R) is not the toroidal angle
        ftors = FUSE.nb_toroidal_velocity_fraction(Xs, Ys, vx, vy)
        @test all(f -> 0.0 <= f <= 1.0 + 1e-12, ftors)
        @test Rs .* ftors ≈ fill(tangency_radius * cos(tilt), length(tt)) rtol = 1e-10
    end
end
