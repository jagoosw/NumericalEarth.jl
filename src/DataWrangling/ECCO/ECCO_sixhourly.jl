using MeshArrays
using JLD2

# URL for 6-hourly input forcing files
const ECCO4_sixhourly_url = "https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/input_forcing/"

# Variable name mapping: symbol → file prefix (used in eccov4r4_{prefix}_{YYYY})
ECCO4_sixhourly_variable_names = Dict(
    :downwelling_longwave  => "dlw",
    :downwelling_shortwave => "dsw",
    :rain_freshwater_flux  => "rain",
    :air_specific_humidity => "spfh2m",
    :air_temperature       => "tmp2m_degC",
    :eastward_stress       => "ustr",
    :northward_stress      => "vstr",
    :wind_speed            => "wspeed",
    :sea_level_pressure    => "pres",
)

# Output grid size after MeshArrays regridding (same lat-lon resolution as ECCO4Monthly)
Base.size(::ECCO4Sixhourly, variable) = (720, 360, 1)

# Coordinate interfaces
longitude_interfaces(::ECCO4Sixhourly) = (-180, 180)
# latitude_interfaces inherited from ECCODataset: (-90, 90)
z_interfaces(::ECCO4Sixhourly) = [0, 1]

# Date range: 6-hourly from 1992 to 2017
all_dates(::ECCO4Sixhourly, variable) = DateTime(1992, 1, 1) : Hour(6) : DateTime(2017, 12, 31, 18)

# Variable mapping
dataset_variable_name(data::Metadata{<:ECCO4Sixhourly}) = ECCO4_sixhourly_variable_names[data.name]
available_variables(::ECCO4Sixhourly) = ECCO4_sixhourly_variable_names

# Data properties
is_three_dimensional(::Metadata{<:ECCO4Sixhourly}) = false
reversed_vertical_axis(::ECCO4Sixhourly) = false
default_mask_value(::ECCO4Sixhourly) = 0

# Grid and binary format (same LLC90 grid as ECCO Darwin)
binary_data_grid(::ECCO4Sixhourly) = GridSpec(ID=:LLC90)
binary_data_size(::ECCO4Sixhourly) = (90, 1170)  # 2D compact format

# Download directory
function default_download_directory(::ECCO4Sixhourly)
    path = joinpath(download_ECCO_cache, "v4", "sixhourly")
    return mkpath(path)
end

# File naming: yearly binary files (e.g., eccov4r4_dlw_1992)
function metadata_filename(metadata::Metadatum{<:ECCO4Sixhourly})
    shortname = dataset_variable_name(metadata)
    yearstr = string(Dates.year(metadata.dates))
    return "eccov4r4_" * shortname * "_" * yearstr
end

# URL construction
function metadata_url(m::Metadatum{<:ECCO4Sixhourly})
    return ECCO4_sixhourly_url * metadata_filename(m)
end

# Sign convention:
# dlw and dsw are described as "positive to decrease ocean temperature"
# PrescribedAtmosphere uses positive-downwelling (positive = warming ocean)
# TODO: verify sign by comparing monthly-averaged 6-hourly data against ECCO4Monthly EXFlwdn/oceQsw
function conversion_units(metadatum::Metadatum{<:ECCO4Sixhourly})
    varname = dataset_variable_name(metadatum)
    if varname in ("dlw", "dsw")
        return InverseSign()
    else
        return nothing
    end
end

"""
    download_dataset(metadata::Metadata{<:ECCO4Sixhourly})

Download 6-hourly ECCO forcing files. Since many metadatums (all 6-hourly timestamps
in a year) share the same yearly file, this deduplicates downloads to one per year.
"""
function download_dataset(metadata::Metadata{<:ECCO4Sixhourly})
    username = get(ENV, "ECCO_USERNAME", nothing)
    password = get(ENV, "ECCO_WEBDAV_PASSWORD", nothing)
    dir = metadata.dir

    # Collect unique yearly files to avoid redundant downloads
    seen = Set{String}()
    unique_metadatums = Metadatum{ECCO4Sixhourly}[]
    for metadatum in metadata
        fname = metadata_filename(metadatum)
        if !(fname in seen)
            push!(seen, fname)
            push!(unique_metadatums, metadatum)
        end
    end

    @root mktempdir(dir) do tmp
        downloader = netrc_downloader(username, password, "ecco.jpl.nasa.gov", tmp)
        ntasks = Threads.nthreads()

        asyncmap(unique_metadatums; ntasks) do metadatum
            fileurl  = metadata_url(metadatum)
            filepath = metadata_path(metadatum)

            if !isfile(filepath)
                instructions_msg = "\n See NumericalEarth.jl/src/DataWrangling/ECCO/README.md for instructions."
                if isnothing(username)
                    msg = "Could not find the ECCO_USERNAME environment variable." * instructions_msg
                    throw(ArgumentError(msg))
                elseif isnothing(password)
                    msg = "Could not find the ECCO_WEBDAV_PASSWORD environment variable." * instructions_msg
                    throw(ArgumentError(msg))
                end
                @info "Downloading ECCO 6-hourly data: $(metadatum.name) year $(Dates.year(metadatum.dates)) to $(metadatum.dir)..."
                Downloads.download(fileurl, filepath; downloader, progress=download_progress)
            end
        end
    end

    return nothing
end

"""
    retrieve_data(metadata::Metadatum{<:ECCO4Sixhourly})

Read a single 6-hourly record from a yearly LLC90 compact binary file
and regrid to a regular lat-lon grid using MeshArrays.

The yearly file contains all 6-hourly records for one year (1460 for non-leap,
1464 for leap years). The correct record is located by byte offset based on the date.
"""
function retrieve_data(metadata::Metadatum{<:ECCO4Sixhourly})
    native_size = binary_data_size(metadata.dataset)  # (90, 1170)
    native_grid_spec = binary_data_grid(metadata.dataset)

    # Calculate record index within the yearly file
    year_start = DateTime(Dates.year(metadata.dates), 1, 1)
    record_index = Int(Dates.value(metadata.dates - year_start) ÷ Dates.value(Hour(6))) + 1

    # Read the specific 2D record via byte offset
    record_numel = prod(native_size)  # 90 * 1170 = 105300
    byte_offset = (record_index - 1) * record_numel * sizeof(Float32)

    native_data = zeros(Float32, record_numel)
    filepath = metadata_path(metadata)
    open(filepath, "r") do io
        seek(io, byte_offset)
        read!(io, native_data)
    end
    native_data = bswap.(native_data)

    # Convert compact array to MeshArray (2D)
    meshed_data = MeshArrays.read(reshape(native_data, native_size...), native_grid_spec)

    Nx, Ny, _, _ = size(metadata)  # (720, 360, 1, 1)
    data = zeros(Float32, Nx, Ny)
    mask = zeros(Float32, Nx, Ny)

    # Load native grid coordinates (cached by MeshArrays)
    native_grid_coords = GridLoad(native_grid_spec; option="full")

    # Interpolation coefficients (cached to disk for reuse)
    interp_file = joinpath(dirname(filepath), "native_interp_coeffs.jld2")
    if !isfile(interp_file)
        resolution_X = 360 / Nx
        resolution_Y = 180 / Ny
        longitudes = longitude_interfaces(metadata.dataset)
        latitudes  = latitude_interfaces(metadata.dataset)
        lon = [i for i = longitudes[1]+resolution_X/2:resolution_X:longitudes[2]-resolution_X/2,
                     j = latitudes[1]+resolution_Y/2:resolution_Y:latitudes[2]-resolution_Y/2]
        lat = [j for i = longitudes[1]+resolution_X/2:resolution_X:longitudes[2]-resolution_X/2,
                     j = latitudes[1]+resolution_Y/2:resolution_Y:latitudes[2]-resolution_Y/2]
        coeffs = interpolation_setup(; Γ=native_grid_coords, lat, lon, filename=interp_file)
    else
        coeffs = interpolation_setup(interp_file)
    end

    # Surface ocean mask from native grid
    native_grid_fac_center = GridLoadVar("hFacC", native_grid_spec)

    # Interpolate 2D field (no depth loop — surface only)
    _, _, c = MeshArrays.Interpolate(meshed_data, coeffs)
    data .= c

    _, _, c = MeshArrays.Interpolate(land_mask(native_grid_fac_center[:, 1]), coeffs)
    mask .= c

    # Fill NaNs with default mask value
    data[isnan.(data)] .= default_mask_value(metadata.dataset)

    return data .* mask
end
