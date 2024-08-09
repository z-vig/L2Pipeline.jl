"""
    L2Pipeline.jl

A pipeline for processing Moon Mineralogy Mapper data from L1 to L2 in Julia.
"""

module L2Pipeline

using HDF5
using Dates
using ProgressMeter

GOOD_BANDS = 9:255

@kwdef mutable struct L1CalData
    """
    Struct containing all of the data that is needed to run the processing pipeline
    """
    wvl :: Vector{<:AbstractFloat}
    rdn :: Array{<:AbstractFloat,3}
    obs :: Array{<:AbstractFloat,3}
    solspec :: Tuple{Vector{<:AbstractFloat},Vector{<:AbstractFloat},Float64} #(wvl,irr,distance)
    statpol :: Tuple{Vector{<:AbstractFloat},Vector{<:AbstractFloat},Dates.DateTime} 
    falpha :: Matrix{<:AbstractFloat}

    current_step :: Array{<:AbstractFloat,3} = rdn
end

function make3d(im::Array{Vector{Float64},2})
    return permutedims([im[I][k] for k=eachindex(im[1,1]),I=CartesianIndices(im)],(2,3,1))
end

function caldata_from_url(url::String,imageid::String)
    """
    Function for creating a caldata object from a directory of datasets. The directory structure must be:
    --> xx_Caldata // where xx is the image ID
        --> xx_rdn.hdr // ENVI header file
        --> xx_spectral_data.hdf5 // HDF5 file containing spectral datasets, see below for structure
        --> xx_solspec.tab // Solar spectrum values
        --> xx_L1B.lbl // PDS label file
        --> xx_statpol1.tab // First statistical polishing coefficients
        --> xx_statpol2.tab // Second statistical polishing coefficients
        --> xx_falpha.tab // Phase function coefficients

    The spectral data hdf5 file should be structured as:
    --> Backplanes
        --> LatLongElev
        --> ObsGeometry
        ...various terrain and phase angle maps
    --> ScalarDatasets
        ...derived scalar datasets
    --> ShadowMaps
        --> facet_shadows
        --> lowsignal_shadows
    --> VectorDatasets
        --> Radiance
        --> Radiance_Smooth
        --> Reflectance
        --> Reflectance_GNDTRU
        --> Reflectance_GNDTRU_Smooth
        --> 2pContinuum
        --> 2pContinuum_GNDTRU_Smooth
    """

    function filt_fxn(c::Char)
        if isdigit(c)
            return true
        elseif isequal(c,'.')
            return true
        else
            return false
        end
    end

    #Getting wvl
    wvl = open(joinpath(url,"$(imageid)_rdn.hdr")) do f
        raw = readlines(f)
        nbands = [parse(Int,filter(isdigit,i)) for i in raw if occursin("bands",i)][1]
        wvl_idx = [n for (n,i) in enumerate(raw) if occursin("wavelength",i)][1]

        return raw[wvl_idx+1:wvl_idx+nbands] .|> 
               x->filter(filt_fxn,x) .|> 
               x->parse(Float64,x)
    end

    #Getting rdn
    rdn = h5open(joinpath(url,"$(imageid)_spectral_data.hdf5"),"r") do f
        return read(f["VectorDatasets/Radiance"])
    end

    #Getting obs
    obs = h5open(joinpath(url,"$(imageid)_spectral_data.hdf5"), "r") do f
        return read(f["Backplanes/ObsGeometry"])
    end

    #Getting solspec
    solwvl,solirr = open(joinpath(url,"$(imageid)_solspec.tab")) do f
        raw = readlines(f)
        raw_split = [split(i," ") for i in raw]
        nums = [parse.(Float64,i[i .!= ""]) for i in raw_split]
        wvl = [i[1] for i in nums]
        irr = [i[2] for i in nums]
        return wvl,irr
    end
    
    soldist = open(joinpath(url,"$(imageid)_L1B.lbl")) do f
        return [parse(Float64,filter(filt_fxn,i)) for i in readlines(f) if occursin("SOLAR_DISTANCE",i)][1]
    end

    #Getting statistical polishing factors
    sp1 = open(joinpath(url,"$(imageid)_statpol1.tab")) do f
        strs = [split(i,"  ") for i in readlines(f)]
        num_strs = [i[i .!= ""] for i in strs]
        nums = [parse.(Float64,i) for i in num_strs]
        g_sp = [i[3] for i in nums]
        return g_sp
    end

    sp2 = open(joinpath(url,"$(imageid)_statpol2.tab")) do f
        strs = [split(i,"  ") for i in readlines(f)]
        num_strs = [i[i .!= ""] for i in strs]
        nums = [parse.(Float64,i) for i in num_strs]
        g_sp = [i[3] for i in nums]
        return g_sp
    end

    starttime = open(joinpath(url,"$(imageid)_L1B.lbl")) do f
        str = [filter(isdigit,i) for i in readlines(f) if occursin("START_TIME",i)][1]
        return DateTime(str,dateformat"yyyymmddHHMMSS")
    end

    #Getting falpha
    falpha_vecs = open(joinpath(url,"$(imageid)_falpha.tab")) do f
        raw = readlines(f)
        splt = [split(i," ") for i in raw]
        return [parse.(Float64,i[i .!= ""]) for i in splt]
    end

    f_alpha = zeros(113,256)
    for (n,i) in enumerate(falpha_vecs[2:114])
        f_alpha[n,:] = i[2:end]
    end

    return L1CalData(
        wvl = wvl[GOOD_BANDS],
        rdn = rdn[:,:,GOOD_BANDS],
        obs = obs[:,:,:],
        solspec = (solwvl[GOOD_BANDS],solirr[GOOD_BANDS],soldist),
        statpol = (sp1[GOOD_BANDS],sp2[GOOD_BANDS],starttime),
        falpha = f_alpha[:,GOOD_BANDS]
    )

end

include("geometry_objects.jl")
export M3Geometry,convert_to_rad!,calc_e,calc_i,get_geom_fromfile,get_geom_fromDEM

include("spectral_utilities.jl")
export movingavg

include("solspec_removal.jl")
export rem_solspec

include("statistical_polishing.jl")
export statistical_polish

include("photometric_correction.jl")
export photometric_correction

include("thermal_corrections.jl")
export clark_etal,li_milliken,B

export L1CalData,
       caldata_from_url

end #L2Pipeline