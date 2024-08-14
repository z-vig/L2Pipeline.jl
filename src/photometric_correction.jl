"""
Running photometric correction with the defined photometric correction coefficients found within L1 calibration directory
"""

function XL(i,e)
    """
    Lomell-Seeliger Photmoetric Limb Darkening Factor
    """
    return cos(i*(π/180))/(cos(e*(π/180))+cos(i*(π/180)))
end

function photometric_coef(dat::L1CalData)
    """
    Reads in calibration data for the falpha coefficients and the geometry of the M3 scene to calculate the photometric correction
    """
    
    i_topo = dat.illum[:,:,1]
    e_topo = dat.illum[:,:,2]
    g_topo = dat.illum[:,:,3]

    ax = axes(dat.current_step)
    correction_coefficients = @showprogress map(CartesianIndices(ax[1:2])) do ind
        x,y = Tuple(ind)
        i = i_topo[x,y]
        e = e_topo[x,y]
        g = g_topo[x,y]
        if isfinite(i)

            if round(Int,g) > 100#size(dat.falpha,1)
                alph = 100#size(dat.falpha,1)
            else
                alph = round(Int,g)
            end

            if cos(i*(π/180)) < 0.05
                xl = XL(30,0)/XL(acos(0.05),e)
            elseif cos(i*(π/180)) > 0.05
                xl = XL(30,0)/XL(i,e)
            end

            falph = dat.falpha[30,:]./dat.falpha[alph,:]
            return xl.*falph
        elseif isnan(i)
            return NaN * ones(size(dat.falpha,2))
        end
    end

end

function photometric_correction!(dat::L1CalData; return_corrected_values :: Bool = false)
    """
    Reads in calibration data for the falpha coefficients and the geometry of the M3 scene to calculate the photometric correction
    """
    
    i_topo = dat.illum[:,:,1]
    e_topo = dat.illum[:,:,2]
    g_topo = dat.illum[:,:,3]

    ax = axes(dat.current_step)
    correction_coefficients = @showprogress map(CartesianIndices(ax[1:2])) do ind
        x,y = Tuple(ind)
        i = i_topo[x,y]
        e = e_topo[x,y]
        g = g_topo[x,y]
        if isfinite(i)
            if round(Int,g) > 100#size(dat.falpha,1)
                alph = 100#size(dat.falpha,1)
            else
                alph = round(Int,g)
            end
            xl = XL(30,0)/XL(i,e)
            falph = dat.falpha[30,:]./dat.falpha[alph,:]
            return xl.*falph
        elseif isnan(i)
            return NaN
        end
    end

    corrected = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        return dat.current_step[x,y,:] .* correction_coefficients[x,y]
    end

    dat.current_step = make3d(corrected)
    println("Photometric Correction Complete!")

    if return_corrected_values == false
        return nothing
    elseif return_corrected_values == true
        return corrected
    end

end