"""
Running photometric correction with the defined photometric correction coefficients found within L1 calibration directory
"""

function photometric_correction(dat::L1CalData, geom::M3Geometry)
    """
    Reads in calibration data for the falpha coefficients and the geometry of the M3 scene to calculate the photometric correction
    """
    
    convert_to_rad!(geom)
    i_topo = calc_i(geom)
    e_topo = calc_e(geom)
    
    function XL(i,e)
        return cos(i*(π/180))/(cos(e*(π/180))+cos(i*(π/180)))
    end

    photometric_corrected = @showprogress map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        i = i_topo[x,y]
        e = e_topo[x,y]
        if round(Int,i+e) > 100#size(dat.falpha,1)
            alph = 100#size(dat.falpha,1)
        else
            alph = round(Int,i+e)
        end
        xl = XL(30,0)/XL(i,e)
        falph = dat.falpha[30,:]./dat.falpha[alph,:]
        return xl.*falph
    end

    dat.current_step = make3d(photometric_corrected)
    println("Photometric Correction Complete!")
    return nothing

end