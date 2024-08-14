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

function B(λ::Vector{Float64},T::Float64,ϵ::Union{Float64,Vector{Float64}})
    return ((2*h*c^2) ./ (λ .^ 5)) .* (ϵ ./ (exp.((h*c)./(λ*kᵦ*T)).-1))
end

function get_temp(B,ϵ,λ,F) :: Float64
    # println("$B, $ϵ, $λ, $F")
    return (h*c/(λ*kᵦ)) * (log((2*h*c^2*ϵ/(F*B*λ^5))+1))^-1
end

function get_temp_photometric(B,ϵ,λ,F,Φ,x,y) :: Float64
    # println("$B, $ϵ, $λ, $F, $Φ, $x, $y")
    return (h*c/(λ*kᵦ)) * (log((2*h*c^2*ϵ*Φ/(F*B*λ^5))+1))^-1
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
            return B(dat.wvl .* 10^-9, temp_derived[x,y], ϵ[x,y])./(10^6*dat.solspec[2]) .* dat.solspec[3]^2
        elseif isnan(temp_derived[x,y])
            return NaN .* ones(size(dat.wvl))
        end
    end
    println("Step 4: initial emissivity determined and temperature derived...")

    #Removing initial thermal emission (Step5)
    IOF1 = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        if isfinite(temp_derived[x,y])
            emiss1 = B(dat.wvl .*10^-9,temp_derived[x,y],ϵ[x,y]) .* dat.solspec[3]^2
            return step1[x,y,:] .- (emiss1./(10^6*dat.solspec[2]))
        elseif isnan(temp_derived[x,y])
            return step1[x,y,:]
        end
    end
    println("Step 5: Initial thermal emission removed...")
    
    #Phase Angle Correction (Step6)
    photo_coef = photometric_coef(dat)
    photo_corrected = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        if isfinite(photo_coef[x,y][1])
            corrected = photo_coef[x,y] .* IOF1[x,y]
            # return corrected
            if maximum(corrected) > 0.6
                return 0.6 .* corrected ./ maximum(corrected) #scaling large reflectance values!
            elseif maximum(corrected) < 0.6
                return corrected
            end
        elseif isnan(photo_coef[x,y][1])
            return NaN .* ones(size(dat.wvl))
        end
    end
    println("Step 6: Photometric correction applied...")

    #Iterative thermal removal function

    #Second projection (Step7)
    projIF2 = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        m = (IOF1[x,y][idxE]-IOF1[x,y][idxD]) / (wvlE-wvlD)
        projIF2 = m*(wvlC-wvlD)+IOF1[x,y][idxD]
        return projIF2
    end


    println("Step 7: Projecting to 2.59μm...")

    #Wavelength-dependent emissivity and second thermal component (Step8)
    ϵ_λ = map(CartesianIndices(ax[1:2])) do i
        x,y=Tuple(i)
        return 1 .- photo_corrected[x,y]
    end

    T2 = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        T2 = step1[x,y] - projIF2[x,y]
        if T2>0
            return T2
        elseif T2<0
            return NaN #Error value when temperature is undetermined
        elseif isnan(T2)
            return NaN
        end
    end
    println("Step 8: Second thermal component obtained...")

    #2nd Derived temperature and thermal emission spectrum (step9)
    temp_derived2 = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        if isfinite(T2[x,y])
            Fidx = argmin(abs.(dat.solspec[1].-wvlC))
            F = 10^6 .* dat.solspec[2][Fidx] ./ π
            return get_temp_photometric(T2[x,y],ϵ_λ[x,y][idxC],wvlC*10^-9,F,photo_coef[x,y][idxC],x,y)
        elseif isnan(T2[x,y])
            return NaN
        end
    end

    planck2 = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        if isfinite(temp_derived2[x,y])
            therm_spectrum = (dat.solspec[3]^2 .* B(dat.wvl .* 10^-9, temp_derived2[x,y], ϵ_λ[x,y])) ./ (10^6*dat.solspec[2])
            return photo_coef[x,y] .* therm_spectrum
        elseif isnan(temp_derived2[x,y])
            return NaN .* ones(size(dat.wvl))
        end
    end
    println("Step 9: Second temperature estimate obtained...")

    #Removing 2nd thermal spectrum and photometrically correcting (step10)
    IOF2 = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        if isfinite(temp_derived2[x,y])
            return step1[x,y,:] .- planck2[x,y]
        elseif isnan(temp_derived2[x,y])
            return step1[x,y,:]
        end
    end

    #photometrically correcting IOF2 (step10.5)
    IOF2_photo = map(CartesianIndices(ax[1:2])) do i
        x,y = Tuple(i)
        if isfinite(photo_coef[x,y][1])
            corrected = photo_coef[x,y] .* IOF2[x,y]
            # return corrected
            if maximum(corrected) > 0.6
                return 0.6 .* corrected ./ maximum(corrected) #scaling large reflectance values!
            elseif maximum(corrected) < 0.6
                return corrected
            end
        elseif isnan(photo_coef[x,y][1])
            return NaN .* ones(size(dat.wvl))
        end
    end


    function iter_thermal(IOF::Matrix{Vector{Float64}},photo_coef::Matrix{Vector{Float64}},orig_rfl::Array{Float64,3},bad_temp_mask::Matrix{Bool})
        projIF = map(CartesianIndices(ax[1:2])) do i
            x,y = Tuple(i)
            m = (IOF[x,y][idxE]-IOF[x,y][idxD]) / (wvlE-wvlD)
            projIF = m*(wvlC-wvlD)+IOF[x,y][idxD]
            return projIF
        end

        #Wavelength-dependent emissivity and second thermal component (Step8)
        ϵ_λ = map(CartesianIndices(ax[1:2])) do i
            x,y=Tuple(i)
            return 1 .- (photo_coef[x,y] .* IOF[x,y]) 
        end

        T = map(CartesianIndices(ax[1:2])) do i
            x,y = Tuple(i)
            T = step1[x,y] - projIF[x,y]
            if T>0
                return T
            elseif T<0
                return NaN #Error value when temperature is undetermined
            elseif isnan(T2)
                return NaN
            end
        end

        #2nd Derived temperature and thermal emission spectrum (step9)
        temp_derived = map(CartesianIndices(ax[1:2])) do i
            x,y = Tuple(i)
            if isfinite(T[x,y])
                Fidx = argmin(abs.(dat.solspec[1].-wvlC))
                F = 10^6 .* dat.solspec[2][Fidx] ./ π
                return get_temp_photometric(T[x,y],ϵ_λ[x,y][idxC],wvlC*10^-9,F,photo_coef[x,y][idxC],x,y)
            elseif isnan(T[x,y])
                return NaN
            end
        end

        planck = map(CartesianIndices(ax[1:2])) do i
            x,y = Tuple(i)
            if isfinite(temp_derived[x,y])
                therm_spectrum = (dat.solspec[3]^2 .* B(dat.wvl .* 10^-9, temp_derived2[x,y], ϵ_λ[x,y])) ./ (10^6*dat.solspec[2])
                return photo_coef[x,y] .* therm_spectrum
            elseif isnan(temp_derived[x,y])
                return NaN .* ones(size(dat.wvl))
            end
        end

        #Removing 2nd thermal spectrum and photometrically correcting (step10)
        IOF2 = map(CartesianIndices(ax[1:2])) do i
            x,y = Tuple(i)
            if isfinite(temp_derived2[x,y])
                return step1[x,y,:] .- planck[x,y]
            elseif isnan(temp_derived2[x,y])
                return step1[x,y,:]
            end
        end

        #photometrically correcting IOF2 (step10.5)
        IOF2_photo = map(CartesianIndices(ax[1:2])) do i
            x,y = Tuple(i)
            if isfinite(photo_coef[x,y][1])
                corrected = photo_coef[x,y] .* IOF2[x,y]
                # return corrected
                if maximum(corrected) > 0.6
                    return 0.6 .* corrected ./ maximum(corrected) #scaling large reflectance values!
                elseif maximum(corrected) < 0.6
                    return corrected
                end
            elseif isnan(photo_coef[x,y][1])
                return NaN .* ones(size(dat.wvl))
            end
        end

        return projIF2,ϵ_λ,T2,temp_derived2,planck2,IOF2,IOF2_photo
    end

    projIF2,ϵ_λ,T2,temp_derived2,planck2,IOF2,IOF2_photo = iter_thermal(IOF1)


    debug_dict = Dict(
        "wvl_used" => (wvlA,wvlB,wvlC,wvlD,wvlE),
        "smooth_wvl" => avg_λ,
        "idx_used" => (idxA,idxB,idxC,idxD,idxE),
        "step1" => step1, #smoothed, if active
        "step2" => projIF, #scalar array of projected values
        "step3" => (T1,ϵ), #(Thermal component, derived emissivity)
        "step4" => (temp_derived,planck1), #Derived temperature and corresponding planck spectrum
        "step5" => make3d(IOF1), #Thermal spectrum removed
        "step6" => (make3d(photo_coef),make3d(photo_corrected)), #Photometrically corrected data
        "step7" => projIF2, #scalar array of projected values
        "step8" => (T2,make3d(ϵ_λ)), #wavelength dependent emissivity
        "step9" => (temp_derived2,planck2), #2nd temperature and thermal emission
        "step10" => make3d(IOF2),
        # "step10.5" => make3d(IOF2_photo)
    )

    return debug_dict
    
end

function li_milliken(dat::L1CalData)
    println("Test")
end

