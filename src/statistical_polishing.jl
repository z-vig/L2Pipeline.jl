"""
Statistical Polishing Step
"""

function statistical_polish(dat::L1CalData)
    starttime = dat.statpol[3]
    ax = axes(dat.rdn)

    if (starttime>=DateTime("2009-01-19T00:00:00") && 
        starttime<DateTime("2009-02-15T00:00:00")) ||
       (starttime>=DateTime("2009-04-15T00:00:00") &&
        starttime<DateTime("2009-04-28T00:00:00")) ||
       (starttime>=DateTime("2009-07-12T00:00:00") &&
        starttime<DateTime("2009-08-17T00:00:00"))

        statpol = map(CartesianIndices(ax[1:2])) do i
            x,y = Tuple(i)
            return dat.current_step[x,y,:] .* dat.statpol[1]
        end
        print("Cold instrument polishing applied...")

    else

        statpol = map(CartesianIndices(ax[1:2])) do i
            x,y = Tuple(i)
            return dat.current_step[x,y,:] .* dat.statpol[2]
        end
        print("Warm instrument polishing applied...")

    end

    dat.current_step = make3d(statpol)
    println("Statistical Polishing Complete!")
end