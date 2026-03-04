using NumericalEarth.DataWrangling: DatasetBackend
using Oceananigans.OutputReaders
using NumericalEarth.Atmospheres: PrescribedAtmosphere, TwoBandDownwellingRadiation

"""
    ECCOPrescribedAtmosphere([architecture = CPU(), FT = Float32];
                              dataset = ECCO4Monthly(),
                              start_date = first_date(dataset, :air_temperature),
                              end_date = last_date(dataset, :air_temperature),
                              dir = default_download_directory(dataset),
                              time_indices_in_memory = 10,
                              time_indexing = Cyclical(),
                              surface_layer_height = 2,  # meters
                              other_kw...)

Return a [`PrescribedAtmosphere`](@ref) representing ECCO state estimate data.
The atmospheric data will be held in `FieldTimeSeries` objects containing
- velocities: u, v
- air temperature and humidity: T, q
- surface pressure: p
- freshwater flux: rain
- downwelling radiation: ℐꜜˢʷ, ℐꜜˡʷ

When `dataset` is an `ECCO4Sixhourly`, the 6-hourly forcing data is used.
Wind velocity components are reconstructed from wind speed and wind stress direction:
`u = wspeed * τx / |τ|`, `v = wspeed * τy / |τ|`.
"""
function ECCOPrescribedAtmosphere(architecture = CPU(), FT = Float32;
                                  dataset = ECCO4Monthly(),
                                  start_date = first_date(dataset, :air_temperature),
                                  end_date = last_date(dataset, :air_temperature),
                                  dir = default_download_directory(dataset),
                                  time_indexing = Cyclical(),
                                  time_indices_in_memory = 10,
                                  surface_layer_height = 2,  # meters
                                  other_kw...)

    if dataset isa ECCO4Sixhourly
        return _ecco_sixhourly_prescribed_atmosphere(architecture, FT;
                   dataset, start_date, end_date, dir,
                   time_indexing, time_indices_in_memory,
                   surface_layer_height, other_kw...)
    end

    ua_meta = Metadata(:eastward_wind;         dataset, start_date, end_date, dir)
    va_meta = Metadata(:northward_wind;        dataset, start_date, end_date, dir)
    Ta_meta = Metadata(:air_temperature;       dataset, start_date, end_date, dir)
    qa_meta = Metadata(:air_specific_humidity; dataset, start_date, end_date, dir)
    pa_meta = Metadata(:sea_level_pressure;    dataset, start_date, end_date, dir)
    ℐꜜˡʷ_meta = Metadata(:downwelling_longwave;  dataset, start_date, end_date, dir)
    ℐꜜˢʷ_meta = Metadata(:downwelling_shortwave; dataset, start_date, end_date, dir)
    Fr_meta = Metadata(:rain_freshwater_flux;  dataset, start_date, end_date, dir)

    kw = (; time_indices_in_memory, time_indexing)
    kw = merge(kw, other_kw)

    ua = FieldTimeSeries(ua_meta, architecture; kw...)
    va = FieldTimeSeries(va_meta, architecture; kw...)
    Ta = FieldTimeSeries(Ta_meta, architecture; kw...)
    qa = FieldTimeSeries(qa_meta, architecture; kw...)
    pa = FieldTimeSeries(pa_meta, architecture; kw...)
    ℐꜜˡʷ = FieldTimeSeries(ℐꜜˡʷ_meta, architecture; kw...)
    ℐꜜˢʷ = FieldTimeSeries(ℐꜜˢʷ_meta, architecture; kw...)
    Fr = FieldTimeSeries(Fr_meta, architecture; kw...)

    auxiliary_freshwater_flux = nothing
    freshwater_flux = (; rain = Fr)

    times = ua.times
    grid  = ua.grid

    velocities = (u = ua, v = va)
    tracers = (T = Ta, q = qa)
    pressure = pa

    downwelling_radiation = TwoBandDownwellingRadiation(shortwave=ℐꜜˢʷ, longwave=ℐꜜˡʷ)

    FT = eltype(ua)
    surface_layer_height = convert(FT, surface_layer_height)

    atmosphere = PrescribedAtmosphere(grid, times;
                                      velocities,
                                      freshwater_flux,
                                      auxiliary_freshwater_flux,
                                      tracers,
                                      downwelling_radiation,
                                      surface_layer_height,
                                      pressure)

    return atmosphere
end

function _ecco_sixhourly_prescribed_atmosphere(architecture, FT;
            dataset, start_date, end_date, dir,
            time_indexing, time_indices_in_memory,
            surface_layer_height, other_kw...)

    # Metadata for all 6-hourly atmospheric variables
    Ta_meta  = Metadata(:air_temperature;       dataset, start_date, end_date, dir)
    qa_meta  = Metadata(:air_specific_humidity;  dataset, start_date, end_date, dir)
    pa_meta  = Metadata(:sea_level_pressure;     dataset, start_date, end_date, dir)
    ℐꜜˡʷ_meta = Metadata(:downwelling_longwave;  dataset, start_date, end_date, dir)
    ℐꜜˢʷ_meta = Metadata(:downwelling_shortwave; dataset, start_date, end_date, dir)
    Fr_meta  = Metadata(:rain_freshwater_flux;   dataset, start_date, end_date, dir)
    ws_meta  = Metadata(:wind_speed;             dataset, start_date, end_date, dir)
    τx_meta  = Metadata(:eastward_stress;        dataset, start_date, end_date, dir)
    τy_meta  = Metadata(:northward_stress;       dataset, start_date, end_date, dir)

    kw = (; time_indices_in_memory, time_indexing)
    kw = merge(kw, other_kw)

    # Build FieldTimeSeries for each variable
    Ta  = FieldTimeSeries(Ta_meta,  architecture; kw...)
    qa  = FieldTimeSeries(qa_meta,  architecture; kw...)
    pa  = FieldTimeSeries(pa_meta,  architecture; kw...)
    ℐꜜˡʷ = FieldTimeSeries(ℐꜜˡʷ_meta, architecture; kw...)
    ℐꜜˢʷ = FieldTimeSeries(ℐꜜˢʷ_meta, architecture; kw...)
    Fr  = FieldTimeSeries(Fr_meta,  architecture; kw...)
    ws  = FieldTimeSeries(ws_meta,  architecture; kw...)
    τx  = FieldTimeSeries(τx_meta,  architecture; kw...)
    τy  = FieldTimeSeries(τy_meta,  architecture; kw...)

    # Reconstruct wind velocity components from wind speed and stress direction:
    #   u = wspeed * τx / |τ|
    #   v = wspeed * τy / |τ|
    # We overwrite the stress FieldTimeSeries data in-place to hold velocity.
    # TODO: when time_indices_in_memory < total timesteps, the derived velocity
    # fields cannot be lazily reloaded by DatasetBackend. For now, ensure
    # time_indices_in_memory covers the full date range for wind fields.
    ua = τx  # reuse the FieldTimeSeries container
    va = τy

    for t in eachindex(ws.times)
        τx_data = parent(τx[t])
        τy_data = parent(τy[t])
        ws_data = parent(ws[t])

        # Compute stress magnitude with small epsilon to avoid division by zero
        τ_mag = @. sqrt(τx_data^2 + τy_data^2) + eps(FT)

        parent(ua[t]) .= @. ws_data * τx_data / τ_mag
        parent(va[t]) .= @. ws_data * τy_data / τ_mag
    end

    auxiliary_freshwater_flux = nothing
    freshwater_flux = (; rain = Fr)

    times = Ta.times
    grid  = Ta.grid

    velocities = (u = ua, v = va)
    tracers = (T = Ta, q = qa)
    pressure = pa

    downwelling_radiation = TwoBandDownwellingRadiation(shortwave=ℐꜜˢʷ, longwave=ℐꜜˡʷ)

    FT = eltype(Ta)
    surface_layer_height = convert(FT, surface_layer_height)

    atmosphere = PrescribedAtmosphere(grid, times;
                                      velocities,
                                      freshwater_flux,
                                      auxiliary_freshwater_flux,
                                      tracers,
                                      downwelling_radiation,
                                      surface_layer_height,
                                      pressure)

    return atmosphere
end
