"""
Solar Spectrum Removal Step
"""

function rem_solspec(dat::L1CalData)
    ax = axes(dat.rdn)
    solspec_removed = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        return (10^6 .* dat.current_step[x,y,:] .* π) ./ (10^6 .* dat.solspec[2] ./ dat.solspec[3]^2)
    end

    dat.current_step = make3d(solspec_removed)
    println("I/F Correction Complete!")
    return nothing
end