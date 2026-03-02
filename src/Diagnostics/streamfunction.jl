using Oceananigans.Architectures: CPU, architecture, on_architecture
using Oceananigans.BoundaryConditions: FieldBoundaryConditions, fill_halo_regions!
using Oceananigans.Fields: FieldStatus, compute!, interior
using Oceananigans.Grids: Bounded, Flat
using Oceananigans.Models: seawater_density
using Oceananigans.Operators: Ay

mutable struct StreamfunctionOperand{V, R, FT}
    meridional_transport :: V
    density :: R
    ρmin :: FT
    ρmax :: FT
    Nρ :: Int
    in_sverdrups :: Bool
    type :: String
end

const StreamfunctionField = Field{<:Any, <:Any, <:Any, <:StreamfunctionOperand}

Base.summary(sfo::StreamfunctionOperand) = "StreamfunctionOperand($(sfo.type))"

const ρy_aliases = ("rho-y", "density-latitude")
const zrho_aliases = ("z-rho", "depth-density", "rho-z", "density-depth")
const horizontal_aliases = ("horizontal", "lat-lon", "barotropic")

"""
    Streamfunction(model; type="rho-y", ρmin=1020, ρmax=1032, Nρ=200, in_sverdrups=true,
                   geopotential_height=nothing, condition=nothing, mask=0)

Return a `Field` streamfunction diagnostic.

Supported `type` values:
- `"rho-y"` (aliases: `"density-latitude"`): output shape `(1, Nyᵥ, Nρ)`
- `"z-rho"` (aliases: `"depth-density"`, `"rho-z"`, `"density-depth"`): output shape `(1, Nz, Nρ)`
- `"horizontal"` (aliases: `"lat-lon"`, `"barotropic"`): output shape `(Nx, Nyᵥ, 1)`
"""
function Streamfunction(model;
                        type = "rho-y",
                        ρmin = 1020,
                        ρmax = 1032,
                        Nρ = 200,
                        in_sverdrups = true,
                        geopotential_height = nothing,
                        condition = nothing,
                        mask = 0)

    normalized_type = normalize_streamfunction_type(type)

    (condition === nothing && mask == 0) || throw(ArgumentError("`condition` and `mask` are not yet implemented for Streamfunction."))

    ρmin = convert(eltype(model.grid), ρmin)
    ρmax = convert(eltype(model.grid), ρmax)

    v_transport = Field(model.velocities.v * Ay)
    density = if normalized_type == "horizontal"
        nothing
    else
        Nρ > 0 || throw(ArgumentError("Nρ must be positive, got Nρ=$Nρ."))
        ρmax > ρmin || throw(ArgumentError("ρmax must be greater than ρmin, got ρmin=$ρmin and ρmax=$ρmax."))
        isnothing(geopotential_height) ? Field(seawater_density(model)) :
                                         Field(seawater_density(model; geopotential_height))
    end

    @show Nx, Nyᵥ, Nz = size(interior(v_transport))
    @show Nyᵥ = size(interior(v_transport), 2)
    @show arch = architecture(model.grid)
    output_grid = output_streamfunction_grid(arch, normalized_type, Nx, Nyᵥ, Nz, Nρ, ρmin, ρmax)

    operand = StreamfunctionOperand(v_transport, density, ρmin, ρmax, Nρ, in_sverdrups, normalized_type)
    loc = (Center(), Center(), Center())
    indices = (:, :, :)
    bcs = FieldBoundaryConditions(output_grid, loc)
    data = new_data(output_grid, loc, indices)
    status = FieldStatus()

    return Field(loc, output_grid, data, bcs, indices, operand, status)
end

function compute!(ψ::StreamfunctionField, time=nothing)
    compute_streamfunction!(ψ)
    fill_halo_regions!(ψ)
    return ψ
end

function compute_streamfunction!(ψ::StreamfunctionField)
    operand = ψ.operand

    compute!(operand.meridional_transport)
    v = Array(interior(on_architecture(CPU(), operand.meridional_transport)))

    ψ_out = if operand.type == "horizontal"
        compute_horizontal_streamfunction(v; in_sverdrups = operand.in_sverdrups)
    else
        compute!(operand.density)
        ρ = Array(interior(on_architecture(CPU(), operand.density)))

        if operand.type == "rho-y"
            compute_rhoy_streamfunction(v, ρ, operand.ρmin, operand.ρmax, operand.Nρ;
                                        in_sverdrups = operand.in_sverdrups)
        else
            compute_depth_density_streamfunction(v, ρ, operand.ρmin, operand.ρmax, operand.Nρ;
                                                 in_sverdrups = operand.in_sverdrups)
        end
    end

    interior(ψ) .= on_architecture(architecture(ψ), ψ_out)
    return ψ
end

normalize_streamfunction_type(type::AbstractString) = begin
    lowercase(type) in ρy_aliases && return "rho-y"
    lowercase(type) in zrho_aliases && return "z-rho"
    lowercase(type) in horizontal_aliases && return "horizontal"
    throw(ArgumentError("Unsupported Streamfunction type=\"$type\". Supported types are \"rho-y\", \"z-rho\", and \"horizontal\"."))
end

function output_streamfunction_grid(arch, type, Nx, Nyᵥ, Nz, Nρ, ρmin, ρmax)
    if type == "rho-y"
        return RectilinearGrid(arch;
                               size = (1, Nyᵥ, Nρ),
                               halo = (0, 0, 0),
                               x = (0, 1),
                               y = (1, Nyᵥ),
                               z = (ρmin, ρmax),
                               topology = (Flat, Bounded, Bounded))
    elseif type == "z-rho"
        return RectilinearGrid(arch;
                               size = (1, Nz, Nρ),
                               halo = (0, 0, 0),
                               x = (0, 1),
                               y = (1, Nz),
                               z = (ρmin, ρmax),
                               topology = (Flat, Bounded, Bounded))
    else
        return RectilinearGrid(arch;
                               size = (Nx, Nyᵥ, 1),
                               halo = (0, 0, 0),
                               x = (1, Nx),
                               y = (1, Nyᵥ),
                               z = (0, 1),
                               topology = (Bounded, Bounded, Flat))
    end
end

function face_density_from_center_density(ρᶜᶜᶜ, Nyᵥ)
    Nx, Nyᶜ, Nz = size(ρᶜᶜᶜ)
    ρᶜᶠᶜ = similar(ρᶜᶜᶜ, Nx, Nyᵥ, Nz)

    @inbounds for k = 1:Nz, j = 1:Nyᵥ, i = 1:Nx
        ρᶜᶠᶜ[i, j, k] = if Nyᵥ == Nyᶜ
            ρᶜᶜᶜ[i, j, k]
        elseif j == 1
            ρᶜᶜᶜ[i, 1, k]
        elseif j > Nyᶜ
            ρᶜᶜᶜ[i, Nyᶜ, k]
        else
            0.5 * (ρᶜᶜᶜ[i, j - 1, k] + ρᶜᶜᶜ[i, j, k])
        end
    end

    return ρᶜᶠᶜ
end

@inline function density_bin_index(ρ, ρmin, Δρ, Nρ)
    return clamp(Int(floor((ρ - ρmin) / Δρ)) + 1, 1, Nρ)
end

function compute_rhoy_streamfunction(v, ρᶜᶜᶜ, ρmin, ρmax, Nρ; in_sverdrups = true)
    Nx, Nyᵥ, Nz = size(v)
    ρᶜᶠᶜ = face_density_from_center_density(ρᶜᶜᶜ, Nyᵥ)
    Δρ = (ρmax - ρmin) / Nρ
    transports = zeros(eltype(v), Nyᵥ, Nρ)

    @inbounds for k = 1:Nz, j = 1:Nyᵥ, i = 1:Nx
        q = v[i, j, k]
        ρ = ρᶜᶠᶜ[i, j, k]
        (isfinite(q) && isfinite(ρ)) || continue
        r = density_bin_index(ρ, ρmin, Δρ, Nρ)
        transports[j, r] += q
    end

    ψρy = similar(transports)
    @inbounds for j = 1:Nyᵥ
        running = zero(eltype(transports))
        for r = Nρ:-1:1
            running += transports[j, r]
            ψρy[j, r] = running
        end
    end

    in_sverdrups && (ψρy ./= 1e6)
    return reshape(ψρy, 1, Nyᵥ, Nρ)
end

function compute_depth_density_streamfunction(v, ρᶜᶜᶜ, ρmin, ρmax, Nρ; in_sverdrups = true)
    Nx, Nyᵥ, Nz = size(v)
    ρᶜᶠᶜ = face_density_from_center_density(ρᶜᶜᶜ, Nyᵥ)
    Δρ = (ρmax - ρmin) / Nρ
    transports = zeros(eltype(v), Nz, Nρ)

    @inbounds for k = 1:Nz, j = 1:Nyᵥ, i = 1:Nx
        q = v[i, j, k]
        ρ = ρᶜᶠᶜ[i, j, k]
        (isfinite(q) && isfinite(ρ)) || continue
        r = density_bin_index(ρ, ρmin, Δρ, Nρ)
        transports[k, r] += q
    end

    ψzρ = similar(transports)
    @inbounds for k = 1:Nz
        running = zero(eltype(transports))
        for r = Nρ:-1:1
            running += transports[k, r]
            ψzρ[k, r] = running
        end
    end

    in_sverdrups && (ψzρ ./= 1e6)
    return reshape(ψzρ, 1, Nz, Nρ)
end

function compute_horizontal_streamfunction(v; in_sverdrups = true)
    Nx, Nyᵥ, _ = size(v)
    zonally_resolved_transport = dropdims(sum(v, dims = 3), dims = 3)
    ψxy = similar(zonally_resolved_transport)

    @inbounds for j = 1:Nyᵥ
        running = zero(eltype(zonally_resolved_transport))
        for i = 1:Nx
            running += zonally_resolved_transport[i, j]
            ψxy[i, j] = running
        end
    end

    in_sverdrups && (ψxy ./= 1e6)
    return reshape(ψxy, Nx, Nyᵥ, 1)
end
