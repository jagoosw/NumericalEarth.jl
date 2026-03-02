include("runtests_setup.jl")

using SeawaterPolynomials: TEOS10EquationOfState
using Oceananigans: location
using Oceananigans.Operators: Ay
using Oceananigans.Models: buoyancy_operation
using NumericalEarth.Diagnostics: MixedLayerDepthField, MixedLayerDepthOperand
using NumericalEarth.Diagnostics: Streamfunction, StreamfunctionField, StreamfunctionOperand

for arch in test_architectures, dataset in (ECCO4Monthly(),)
    A = typeof(arch)
    @info "Testing MixedLayerDepthField with $(typeof(dataset)) on $A"

    @testset "MixedLayerDepthField" begin
        grid = LatitudeLongitudeGrid(arch;
                                     size = (3, 3, 100),
                                     latitude  = (0, 30),
                                     longitude = (150, 180),
                                     z = (-1000, 0))

        bottom_height = regrid_bathymetry(grid;
                                          minimum_depth = 10,
                                          interpolation_passes = 5,
                                          major_basins = 1)

        grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

        start = DateTimeProlepticGregorian(1993, 1, 1)
        stop  = DateTimeProlepticGregorian(1993, 2, 1)
        dates = range(start; stop, step=Month(1))

        Tmeta = Metadata(:temperature; dataset, dates)
        Smeta = Metadata(:salinity; dataset, dates)

        Tt = FieldTimeSeries(Tmeta, grid; time_indices_in_memory=2)
        St = FieldTimeSeries(Smeta, grid; time_indices_in_memory=2)

        equation_of_state = TEOS10EquationOfState()
        sb = SeawaterBuoyancy(; equation_of_state)
        tracers = (T=Tt[1], S=St[1])
        h = MixedLayerDepthField(sb, grid, tracers)

        @test h isa Field
        @test location(h) == (Center, Center, Nothing)
        @test h.operand isa MixedLayerDepthOperand
        @test h.operand.buoyancy_perturbation isa KernelFunctionOperation

        compute!(h)
        if dataset isa ECCO4Monthly
            @test @allowscalar h[1, 1, 1] ≈ 16.2558363 # m
        end

        tracers = (T=Tt[2], S=St[2])
        h.operand.buoyancy_perturbation = buoyancy_operation(sb, grid, tracers)
        compute!(h)
        if dataset isa ECCO4Monthly
            @test @allowscalar h[1, 1, 1] ≈ 9.2957298 # m
        end
    end
end

for arch in test_architectures
    A = typeof(arch)
    @info "Testing Streamfunction diagnostic on $A"

    @testset "Streamfunction modes on $A" begin
        grid = RectilinearGrid(arch;
                               size = (8, 6, 5),
                               x = (0, 1),
                               y = (-20, 20),
                               z = (-500, 0),
                               topology = (Periodic, Bounded, Bounded),
                               halo = (3, 3, 3))

        eos = TEOS10EquationOfState()
        buoyancy = SeawaterBuoyancy(; equation_of_state = eos)
        model = HydrostaticFreeSurfaceModel(grid; buoyancy, tracers = (:T, :S))

        set!(model, v = (x, y, z) -> 1e-2, T = 10.0, S = 35.0)

        ψρy = Streamfunction(model; type = "rho-y", ρmin = 1018, ρmax = 1038, Nρ = 32)

        @test ψρy isa Field
        @test ψρy isa StreamfunctionField
        @test ψρy.operand isa StreamfunctionOperand

        compute!(ψρy)
        Ψρy = Array(interior(on_architecture(CPU(), ψρy)))

        @test all(isfinite, Ψρy)
        @test size(Ψρy, 1) == 1

        vAy = Field(model.velocities.v * Ay)
        compute!(vAy)
        transport = dropdims(sum(Array(interior(on_architecture(CPU(), vAy))), dims = (1, 3)); dims = (1, 3)) ./ 1e6

        @test Ψρy[1, :, 1] ≈ transport

        ψzρ = Streamfunction(model; type = "z-rho", ρmin = 1018, ρmax = 1038, Nρ = 32)
        compute!(ψzρ)
        Ψzρ = Array(interior(on_architecture(CPU(), ψzρ)))

        @test all(isfinite, Ψzρ)
        @test size(Ψzρ, 1) == 1
        @test size(Ψzρ, 2) == size(interior(vAy), 3)

        depth_transport = dropdims(sum(Array(interior(on_architecture(CPU(), vAy))), dims = (1, 2)); dims = (1, 2)) ./ 1e6
        @test Ψzρ[1, :, 1] ≈ depth_transport

        ψxy = Streamfunction(model; type = "lat-lon")
        compute!(ψxy)
        Ψxy = Array(interior(on_architecture(CPU(), ψxy)))

        @test all(isfinite, Ψxy)
        @test size(Ψxy, 3) == 1

        zonal_cumulative = dropdims(sum(Array(interior(on_architecture(CPU(), vAy))), dims = 3); dims = 3) ./ 1e6
        @test Ψxy[:, :, 1] ≈ cumsum(zonal_cumulative, dims = 1)

        # Alias coverage
        ψalias = Streamfunction(model; type = "depth-density", ρmin = 1018, ρmax = 1038, Nρ = 32)
        compute!(ψalias)
        @test Array(interior(on_architecture(CPU(), ψalias))) ≈ Ψzρ
    end
end
