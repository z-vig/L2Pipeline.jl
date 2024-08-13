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
    # println("$B, $ϵ, $λ, $F")
    return (h*c/(λ*kᵦ)) * (log((2*h*c^2*ϵ/(F*B*λ^5))-1))^-1
end

function clark_etal!(dat::L1CalData)
    #Undoing Solar Distance Correction (Step1)
    step1 = dat.current_step .* dat.solspec[3]^2
    println("Step 1: Solar distance correction removed...")
    
    #Optional smoothing Step (Step1.5)
    step1,avg_λ = movingavg(step1,dat.wvl,9)
    println("Step 1.5: Spectra smoothed...")

    #Linearly projecting I/F value (Step2)
    wvlA = 1550; idxA = argmin(abs.(dat.wvl.-wvlA))
    wvlB = 2350; idxB = argmin(abs.(dat.wvl.-wvlB))
    wvlC = 2700; idxC = argmin(abs.(dat.wvl.-wvlC))
    wvlD = 2280; idxD = argmin(abs.(dat.wvl.-wvlD))
    wvlE = 2590; idxE = argmin(abs.(dat.wvl.-wvlE))

    ax = axes(step1)
    projIF = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        m = (step1[x,y,idxB]-step1[x,y,idxA]) / (wvlB-wvlA)
        projIF = m*(wvlC-wvlA)+step1[x,y,idxA]
        return projIF #Returns the projected value at wvlC,\. Scalar array
    end
    println("Step 2: Spectra projected to 2.7μm")

    #Determining initial thermal component and emissivity (Step3)
    T1 = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        T1 = step1[x,y,idxC] - projIF[x,y]
        if T1>0
            return T1
        elseif T1<0
            return NaN #Error value when temperature is undetermined
        elseif isnan(T1)
            return NaN
        end
    end
    println("Step 3: First thermal component obtained...")

    
    ϵ = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        if isfinite(step1[x,y,idxA])
            return 1 - step1[x,y,idxA]
        elseif isnan(step1[x,y,idxA])
            return NaN
        end
    end

    #Determining initial temperature and planck spectrum (Step4 + extra)
    temp_derived = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        if isfinite(T1[x,y])
            Fidx = argmin(abs.(dat.solspec[1].-wvlC))
            F = 10^6 .* dat.solspec[2][Fidx] ./ π
            return get_temp(T1[x,y],ϵ[x,y],wvlC*10^-9,F)
        elseif isnan(T1[x,y])
            return NaN
        end
    end

    planck1 = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        if isfinite(temp_derived[x,y])
            return B(dat.wvl .* 10^-9, temp_derived[x,y], ϵ[x,y])./(10^6*dat.solspec[2])
        elseif isnan(temp_derived[x,y])
            return NaN .* ones(size(dat.wvl))
        end
    end
    println("Step 4: initial emissivity determined and temperature derived...")

    #Removing initial thermal emission (Step5)
    IOF1 = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        if isfinite(temp_derived[x,y])
            emiss1 = B(dat.wvl .*10^-9,temp_derived[x,y],ϵ[x,y])
            return step1[x,y,:] .- (emiss1./(10^6*dat.solspec[2]))
        elseif isnan(temp_derived[x,y])
            return step1
        end
    end
    println("Step 5: Initial thermal emission removed...")
    
    #Phase Angle Correction (Step6)
    photo_correction = photometric_correction(
        IOF1,
        dat.illum...,
        dat.falpha,
        return_corrected_values=true
    )

    # println("Step 6: Photometric correction applied...")

    # #Wavelength-dependent emissivity (Step7)
    # ϵ_λ = map(CartesianIndices(ax[1:2])) do i
    #     x,y=Tuple(i)
    #     return 1 .- IOF1c[x,y]
    # end
    # println("Step 7: Wavelength-dependent emissivity derived...")

    # #Second projection (Step8)
    # projIF2 = map(CartesianIndices(ax[1:2])) do i
    #     x,y = Tuple(i)
    #     m = (IOF1[x,y][idxE]-IOF1[x,y][idxD]) / (wvlE-wvlD)
    #     projIF2 = m*(wvlC-wvlD)+IOF1[x,y][idxD]
    #     return projIF2
    # end

    # T2 = map(CartesianIndices(ax[1:2])) do i
    #     x,y = Tuple(i)
    #     T2 = dat.current_step[x,y,idxC] - projIF2[x,y]
    #     if T2>0
    #         return T2
    #     elseif T2<0
    #         return -999 #Error value when temperature is undetermined
    #     end
    # end
    # println("Step 8: Projecting to 2.59μm and obtaining residual thermal component...")

    # #Second blackbody estimation (step9)
    # temp_derived2 = map(CartesianIndices(ax[1:2])) do i
    #     x,y = Tuple(i)
    #     if T2[x,y] != -999
    #         Fidx = argmin(abs.(dat.solspec[1].-wvlC))
    #         F = 10^6 .* dat.solspec[2][Fidx] ./ π
    #         try
    #             return get_temp(T2[x,y],ϵ_λ[x,y][idxC],wvlC*10^-9,F)
    #         catch
    #             println(ϵ_λ[x,y][idxC])
    #         end
    #     else
    #         return -999.
    #     end
    # end

    debug_dict = Dict(
        "wvl_used" => (wvlA,wvlB,wvlC,wvlD,wvlE),
        "smooth_wvl" => avg_λ,
        "idx_used" => (idxA,idxB,idxC,idxD,idxE),
        "step1" => step1, #smoothed, if active
        "step2" => projIF, #scalar array of projected values
        "step3" => (T1,ϵ), #(Thermal component, derived emissivity)
        "step4" => (temp_derived,planck1), #Derived temperature and corresponding planck spectrum
        "step5" => IOF1, #Thermal spectrum removed
        # "step6" => photo_correction #Photometrically corrected data
        # "IOF1" => make3d(IOF1),
        # "i_topo" => i_topo,
        # "e_topo" => e_topo,
        # "photo" => make3d(photo_correction),
        # "IOF1c" => make3d(IOF1c),
        # "ep_wvl" => make3d(ϵ_λ),
        # "p2" => pnts,
        # "T2" => T2
        # "temp2" => temp_derived2
    )

    return debug_dict
    
end

function li_milliken(dat::L1CalData)
    println("Test")
end

