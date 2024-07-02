"""
Here we define some structs to aid in image geometry exploration

1) to-sun azimuth angle (decimal degrees, clockwise from local north)
2) to-sun zenith angle (incidence angle in decimal degrees, zero at zenith)
3) to-sensor azimuth angle (decimal degrees, clockwise from local north)
4) to-sensor zenith angle (emission angle in decimal degrees, zero at zenith)
5) observation phase angle (decimal degrees, in plane of to-sun and to-sensor rays)
6) to-sun path length (decimal au with scene mean subtracted and noted in PDS label)
7) to-sensor path length (decimal meters)
8) surface slope from DEM (decimal degrees, zero at horizontal)
9) surface aspect from DEM (decimal degrees, clockwise from local north)
10) local cosine i (unitless, cosine of angle between to-sun and local DEM facet normal
vectors)
"""

using ArchGDAL, Rasters
import ArchGDAL as AG

@kwdef mutable struct M3Geometry
    """
    Struct defining all relevant geometrical parameters to a M3 image.
    """
    solaz :: Matrix{<:AbstractFloat}
    solze :: Matrix{<:AbstractFloat}
    m3az :: Matrix{<:AbstractFloat}
    m3ze :: Matrix{<:AbstractFloat}
    phase :: Matrix{<:AbstractFloat}
    solen :: Matrix{<:AbstractFloat}
    m3len :: Matrix{<:AbstractFloat}
    slope :: Matrix{<:AbstractFloat}
    aspec :: Matrix{<:AbstractFloat}
    cosi :: Matrix{<:AbstractFloat}
end

function convert_to_rad!(geo :: M3Geometry)
    fields = fieldnames(M3Geometry)

    for i = eachindex(fields)
        val = getfield(geo,fields[i])
        val = val * (π/180)
        setfield!(geo,fields[i],val)
    end

    return nothing
end

function calc_i(geo :: M3Geometry) :: Matrix{<:AbstractFloat}
    return (180/π) .* acos.(cos.(geo.solze) .* cos.(geo.slope) .+ sin.(geo.solze)  .* sin.(geo.slope) .* cos.((geo.solaz .- geo.aspec)))    
end

function calc_e(geo :: M3Geometry) :: Matrix{<:AbstractFloat}
    return (180/π) .* acos.(cos.(geo.m3ze) .* cos.(geo.slope) .+ sin.(geo.m3ze) .* sin.(geo.slope) .* cos.((geo.m3az .- geo.aspec)))
end

function modify_to_DEM(geo :: M3Geometry) :: Matrix{<:AbstractFloat}
    i = calc_i(geo)
    e = calc_e(geo)
    geo.phase = i.+e
    geo.cosi = cos.((π/180).*i)
end

function get_geom_fromfile(path::String)::M3Geometry
    obs = h5open(abspath(path)) do f
        f["Backplanes/ObsGeometry"][:,:,:]
    end

    return M3Geometry(
        solaz = obs[:,:,1],
        solze = obs[:,:,2],
        m3az = obs[:,:,3],
        m3ze = obs[:,:,4],
        phase = obs[:,:,5],
        solen = obs[:,:,6],
        m3len = obs[:,:,7],
        slope = obs[:,:,8],
        aspec = obs[:,:,9],
        cosi = obs[:,:,10]
    )
end

function get_geom_fromDEM(m3path::String,slopedatapath::String,aspecdatapath::String)
    m3 = AG.readraster(m3path)
    slop = AG.readraster(slopedatapath)
    aspe = AG.readraster(aspecdatapath)

    return M3Geometry(
        solaz = m3[:,:,1],
        solze = m3[:,:,2],
        m3az = m3[:,:,3],
        m3ze = m3[:,:,4],
        phase = m3[:,:,5],
        solen = m3[:,:,6],
        m3len = m3[:,:,7],
        slope = slop[:,:],
        aspec = aspe[:,:],
        cosi = m3[:,:,10]
    )

end

get_geom_fromDEM("C:/Lunar_Imagery_Data/M3_data/obsgeom/obsgeom_wall.tif","C:/Lunar_Imagery_Data/LOLA_Slope/wall_slope.tif","C:/Lunar_Imagery_Data/LOLA_Slope/wall_aspect.tif")
GC.gc()