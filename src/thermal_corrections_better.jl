h = 6.626*10^-34 #J*s
kᵦ = 1.381*10^-23 #J/K
c = 2.998*10^8 #m/s

function B(λ::Vector{Float64},T::Float64)
    return ((2*h*c^2) ./ (λ .^ 5)) .* (1 ./ (exp.((h*c)./(λ*kᵦ*T)).-1))
end

function get_temp(B,ϵ,λ,F,r₀) :: Float64
    # println("$B, $ϵ, $λ, $F")
    # for i in [B,ϵ,λ,F,r₀]
    #     println(typeof(i))
    # end
    return (h*c/(λ*kᵦ)) * (log(((2*h*c^2*ϵ*r₀^2)/(F*B*λ^5)) .+ 1))^-1
end

function get_temp_photometric(B,ϵ,λ,F,r₀,Φ) :: Float64
    # println("$B, $ϵ, $λ, $F, $Φ, $x, $y")
    return (h*c/(λ*kᵦ)) * (log((2*h*c^2*ϵ*Φ/(F*B*λ^5))+1))^-1
end

function clark_etal!(dat::L1CalData)

    function project_spectrum(IOF::Array{Float64,3},idx1::Int,idx2::Int,wvl1::Float64,wvl2::Float64,wvl_proj::Float64)
        """
        Draws a line through wvl1 and wvl2 and returns the projected reflectance values at wvl_proj
        """
        ax = axes(IOF)
        projected = map(CartesianIndices(ax[1:2])) do i
            x,y = Tuple(i)
            m = (IOF[x,y,idx2]-IOF[x,y,idx1]) / (wvl2-wvl1)
            return m*(wvl_proj-wvl1)+IOF[x,y,idx1]
        end
        return projected
    end

    function derive_temp(IOF::Array{Float64,3},IOF_proj::Matrix{Float64},idx_proj::Int, λ_proj::Float64; emiss_constant::Bool=false,idx_emiss::Union{Nothing,Int}=nothing, include_photometric::Bool = false)
        """
        Derives temperature from the estimated thermal component of a spectrum
        """
        if emiss_constant == true && isnothing(idx_emiss)
            println("Please specify index for constant emissivity wavelength!")
        end

        ax = axes(IOF)
        therm_est = map(CartesianIndices(ax[1:2])) do i
            x,y = Tuple(i)
            diff = IOF[x,y,idx_proj] - IOF_proj[x,y]
            if diff > 0
                return diff
            else
                return NaN
            end
        end

        ϵ = map(CartesianIndices(ax[1:2])) do i
            x,y = Tuple(i)
            if emiss_constant == true
                return ones(size(IOF,3)) .* (1 - IOF[x,y,idx_emiss])
            elseif emiss_constant == false
                return 1 .- IOF[x,y,:]
            end
        end

        temp = map(CartesianIndices(ax[1:2])) do i 
            x,y = Tuple(i)
            if isfinite(therm_est[x,y])
                Fidx = argmin(abs.(dat.solspec[1].-wvlC))
                F = 10^6 .* dat.solspec[2][Fidx] ./ π
                return get_temp(therm_est[x,y],ϵ[x,y][idx_proj],λ_proj*10^-9,F,dat.solspec[3])
            elseif isnan(therm_est[x,y])
                return NaN
            end
        end

        planck = map(CartesianIndices(ax[1:2])) do i 
            x,y = Tuple(i)
            if isfinite(temp[x,y])
                BT = B(dat.wvl .* 10^-9, temp[x,y])
                F = 10^6 .* dat.solspec[2] ./ π
                return (dat.solspec[3]^2 .* ϵ[x,y] .* BT) ./ F
            elseif isnan(temp[x,y])
                return NaN .* ones(size(dat.wvl))
            end
        end

        thermal_removed = map(CartesianIndices(ax[1:2])) do i
            x,y = Tuple(i)
            if isfinite(planck[x,y][1])
                return IOF[x,y,:] .- planck[x,y]
            elseif isnan(planck[x,y][1])
                return IOF[x,y,:]
            end
        end

        return temp,planck,thermal_removed

    end

    N = 2
    function iterative_thermal(IOF::Array{Float64,3},tempN0::Matrix{Float64},λ1::Float64,λ2::Float64,λ3::Float64,idx1::Int,idx2::Int,idx3::Int)
        final_thermal_removal = Array{Float64,3}(undef,size(IOF))
        final_proj = Matrix{Float64}(undef,size(tempN0))
        final_temp = Matrix{Float64}(undef,size(tempN0))
        final_planck = Array{Float64,3}(undef,size(IOF))

        while true
            proj = project_spectrum(IOF,idx1,idx2,λ1,λ2,λ3)

            tempN1,planck,thermal_removed = derive_temp(IOF,proj,idx3,λ3,emiss_constant=false)

            temp_diff = abs.(tempN1 .- tempN0)
            tempN1[temp_diff.<2] .= NaN
            bad_num = count(temp_diff[isfinite.(temp_diff)].>2)
            
            if all(temp_diff[isfinite.(temp_diff)].<2) || N == 5
                println("Iteration: $N, $bad_num, $(extrema(tempN0[isfinite.(tempN0)]))")
                println("Iteration Complete!")
                final_thermal_removal .= make3d(thermal_removed)
                final_planck .= make3d(planck)
                final_proj .= proj
                final_temp .= tempN1
                break
            else
                println("Iteration: $N, $bad_num, $(extrema(tempN0[isfinite.(tempN0)]))")
                N+=1
                IOF .= make3d(thermal_removed)
                tempN0 .= tempN1
            end
        end

        return final_proj,final_temp,final_planck,final_thermal_removal
        
    end

    #Linearly projecting I/F value
    wvlA = 1550; idxA = argmin(abs.(dat.wvl.-wvlA))
    wvlB = 2350; idxB = argmin(abs.(dat.wvl.-wvlB))
    wvlC = 2700; idxC = argmin(abs.(dat.wvl.-wvlC))
    wvlD = 2280; idxD = argmin(abs.(dat.wvl.-wvlD))
    wvlE = 2590; idxE = argmin(abs.(dat.wvl.-wvlE))

    #Real Wavelength Values, from those defined in Clark et al., 2011
    λA = dat.wvl[idxA]
    λB = dat.wvl[idxB]
    λC = dat.wvl[idxC]
    λD = dat.wvl[idxD]
    λE = dat.wvl[idxE]

    #Undoing Solar Distance Correction
    IOF0 = dat.current_step .* dat.solspec[3]^2

    #Optional smoothing
    IOF0,avg_λ = movingavg(IOF0,dat.wvl,9)

    #Projecting from wvlA and wvlB to wvlC
    proj1 = project_spectrum(IOF0,idxA,idxB,λA,λB,λC)

    #Deriving initial temp estimate
    temp1,planck1,thermal_removed1 = derive_temp(IOF0,proj1,idxC,λC,emiss_constant=true,idx_emiss=idxA)

    Φ = photometric_coef(dat)

    IOF1c = map(CartesianIndices(axes(temp1))) do i
        x,y = Tuple(i)
        corrected = Φ[x,y] .* thermal_removed1[x,y]
        if maximum(corrected) > 0.6
            return 0.6 .* corrected ./ maximum(corrected)
        elseif maximum(corrected) < 0.6
            return corrected
        elseif any(isnan.(corrected))
            return NaN .* ones(size(dat.wvl))
        end
    end

    final_proj, final_temp, final_planck, final_thermal_removed = iterative_thermal(make3d(IOF1c),temp1,λD,λE,λC,idxD,idxE,idxC)

    # println(typeof(IOF1c))

    debug_dict = Dict(
        "λ" => [λA,λB,λC,λD,λE],
        "idx" => [idxA,idxB,idxC,idxD,idxE],
        "IOF0" => IOF0,
        "proj1" => proj1,
        "temp1" => temp1,
        "planck1" => planck1,
        "notherm" => make3d(thermal_removed1),
        "photometric" => make3d(IOF1c),
        "projf" => final_proj,
        "tempf" => final_temp,
        "planckf" => final_planck,
        "nothermf" => final_thermal_removed
    )
    return debug_dict


end