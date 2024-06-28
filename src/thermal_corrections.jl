"""
Here we try some different thermal corrections to M3 data.

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

h = 6.626*10^-34 #J*s
kᵦ = 1.381*10^-23 #J/K
c = 2.998*10^8 #m/s

function B(λ::Vector{Float64},T::Float64,ϵ::Float64)
    return ((2*h*c^2) ./ (λ .^ 5)) .* (ϵ ./ (exp.((h*c)./(λ*kᵦ*T)).-1))
end

function get_temp(B,ϵ,λ,F) :: Float64
    return (h*c/(λ*kᵦ)) * (log((2*h*c^2*ϵ/(F*B*λ^5))-1))^-1
end

function clark_etal(dat::L1CalData)
    #Undoing Solar Distance Correction (Step1)
    dat.current_step = dat.current_step .* dat.solspec[3]^2

    #Optional smoothing Step (Step1.5)
    dat.current_step,avg_λ = movingavg(dat.current_step,dat.wvl,9)

    #Linearly projecting I/F value (Step2)
    wvlA = 1550; idxA = argmin(abs.(dat.wvl.-wvlA))
    wvlB = 2350; idxB = argmin(abs.(dat.wvl.-wvlB))
    wvlC = 2700; idxC = argmin(abs.(dat.wvl.-wvlC))
    wvlD = 2280; idxD = argmin(abs.(dat.wvl.-wvlD))
    wvlE = 2590; idxE = argmin(abs.(dat.wvl.-wvlE))

    ax = axes(dat.current_step)
    projIF = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        m = (dat.current_step[x,y,idxB]-dat.current_step[x,y,idxA]) / (wvlB-wvlA)
        projIF = m*(wvlC-wvlA)+dat.current_step[x,y,idxA]
        return projIF
    end

    #Determining initial thermal component (Step3)
    T1 = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        T1 = dat.current_step[x,y,idxC] - projIF[x,y]
        if T1>0
            return T1
        elseif T1<0
            return -999 #Error value when temperature is undetermined
        end
    end

    #Determining initial emissivity and temperature (Step4 + extra)
    ϵ = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        return 1 - dat.current_step[x,y,idxA]
    end

    temp_derived = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        if T1[x,y] != -999
            Fidx = argmin(abs.(dat.solspec[1].-wvlC))
            F = 10^6 .* dat.solspec[2][Fidx] ./ π
            return get_temp(T1[x,y],ϵ[x,y],wvlC*10^-9,F)
        else
            return -999.
        end
    end

    #Removing initial thermal emission (Step5)
    IOF1 = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        if temp_derived[x,y] != -999.
            emiss1 = B(dat.wvl.*10^-9,temp_derived[x,y],ϵ[x,y])
            return dat.current_step[x,y,:] .- (emiss1./(10^6*dat.solspec[2]))
        else
            return dat.current_step[x,y,:]
        end
    end
    
    #Phase Angle Correction (Step6)
    i_topo = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        saz = dat.obs[x,y,1]*(π/180)
        sze = dat.obs[x,y,2]*(π/180)
        slo = dat.obs[x,y,8]*(π/180)
        asp = dat.obs[x,y,9]*(π/180)
        return (180/π)*acos(cos(sze)*cos(slo)+sin(sze)*sin(slo)*cos((saz-asp)))
    end

    e_topo = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        maz = dat.obs[x,y,3]*(π/180)
        mze = dat.obs[x,y,4]*(π/180)
        slo = dat.obs[x,y,8]*(π/180)
        asp = dat.obs[x,y,9]*(π/180)
        return (180/π)*acos(cos(mze)*cos(slo)+sin(mze)*sin(slo)*cos((maz-asp)))
    end
    
    function XL(i,e)
        return cos(i*(π/180))/(cos(e*(π/180))+cos(i*(π/180)))
    end

    photo_correction = @showprogress map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        i = i_topo[x,y]
        e = e_topo[x,y]
        if round(Int,i+e) > 100#size(dat.falpha,1)
            α = 100#size(dat.falpha,1)
        else
            α = round(Int,i+e)
        end
        xl = XL(30,0)/XL(i,e)
        fα = dat.falpha[30,:]./dat.falpha[α,:]
        return xl.*fα
    end

    IOF1c = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        return IOF1[x,y] .* photo_correction[x,y]
    end

    IOF_notherm = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        return dat.current_step[x,y,:] .* photo_correction[x,y]
    end

    #Wavelength-dependent emissivity (Step7)
    ϵ_λ = map(CartesianIndices(ax[1:2])) do i
        x,y=Tuple(i)
        return 1 .- IOF1c[x,y]
    end

    #Second projection (Step8)
    projIF2 = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        m = (IOF1[x,y][idxE]-IOF1[x,y][idxD]) / (wvlE-wvlD)
        projIF2 = m*(wvlC-wvlD)+IOF1[x,y][idxD]
        return projIF2
    end

    T2 = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        T2 = dat.current_step[x,y,idxC] - projIF2[x,y]
        if T2>0
            return T2
        elseif T2<0
            return -999 #Error value when temperature is undetermined
        end
    end

    #Second blackbody estimation
    temp_derived2 = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        if T2[x,y] != -999
            Fidx = argmin(abs.(dat.solspec[1].-wvlC))
            F = 10^6 .* dat.solspec[2][Fidx] ./ π
            try
                return get_temp(T2[x,y],ϵ_λ[x,y][idxC],wvlC*10^-9,F)
            catch
                println(ϵ_λ[x,y][idxC])
            end
        else
            return -999.
        end
    end

    pnts(x,y) = [(dat.wvl[idxD],IOF1[x,y][idxD]),(dat.wvl[idxE],IOF1[x,y][idxE]),(dat.wvl[idxC],projIF2[x,y])]

    ret_dict = Dict(
        "temp" => temp_derived,
        "T1" => T1,
        "ep" => ϵ,
        "IOF1" => make3d(IOF1),
        "i_topo" => i_topo,
        "e_topo" => e_topo,
        "photo" => make3d(photo_correction),
        "IOF1c" => make3d(IOF1c),
        "IOF_notherm" => make3d(IOF_notherm),
        "ep_wvl" => make3d(ϵ_λ),
        "p2" => pnts,
        "T2" => T2
        # "temp2" => temp_derived2
    )

    return ret_dict
    
end

function li_milliken(dat::L1CalData)
    println("Test")
end

