"""
    L2Pipeline.jl

A pipeline for processing Moon Mineralogy Mapper data from L1 to L2 in Julia.
"""

module L2Pipeline

using HDF5

@kwdef mutable struct L1CalData
    """
    Struct containing all of the data that is needed to run the processing pipeline
    """
    wvl :: Vector{<:AbstractFloat}
    rdn :: Array{<:AbstractFloat,3}
    solspec :: Tuple{Vector{<:AbstractFloat},Vector{<:AbstractFloat}}
end

function caldata_from_url(url::String)

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
    wvl = open(joinpath(url,"rdn.hdr")) do f
        raw = readlines(f)
        nbands = [parse(Int,filter(isdigit,i)) for i in raw if occursin("bands",i)][1]
        wvl_idx = [n for (n,i) in enumerate(raw) if occursin("wavelength",i)][1]

        return raw[wvl_idx+1:wvl_idx+nbands] .|> 
               x->filter(filt_fxn,x) .|> 
               x->parse(Float64,x)
    end

    #Getting rdn
    rdn = h5open(joinpath(url,"spectral_data.hdf5"),"r") do f
        return read(f["VectorDatasets/Radiance"])
    end

    #Getting solspec
    solwvl,solirr = open(joinpath(url,"solspec.tab")) do f
        raw = readlines(f)
        raw_split = [split(i," ") for i in raw]
        nums = [parse.(Float64,i[i .!= ""]) for i in raw_split]
        wvl = [i[1] for i in nums]
        irr = [i[2] for i in nums]
        return wvl,irr
    end

    return L1CalData(
        wvl = wvl,
        rdn = rdn,
        solspec = (solwvl,solirr)
    )

end

include("solspec_removal.jl")

export L1CalData,
       caldata_from_url

end #L2Pipeline