"""
Several utility functions for manipulating spectra
"""

using Statistics

function movingavg(input_image::Array{<:AbstractFloat,3},input_λvector::Vector{<:AbstractFloat},box_size::Int)

    if box_size%2==0
        throw(DomainError(box_size,"Box Size must be odd!"))
    end

    split_index::Int = (box_size-1)/2
    avg_im_size = (size(input_image)[1:2]...,size(input_image)[3]-(2*split_index))
    avg_im = zeros(avg_im_size)

    for band ∈ 1:size(avg_im)[3]
        subset_img = input_image[:,:,band:band+(2*split_index)]
        av_subset = mean(subset_img,dims=3)
        sd_subset = std(subset_img,dims=3)
        upperlim_subset = av_subset.+(2*sd_subset)
        lowerlim_subset = av_subset.-(2*sd_subset)
        
        subset_img[(subset_img.<lowerlim_subset).||(subset_img.>upperlim_subset)].=0.0
        wiseav_missingvals = convert(Array{Union{Float64,Missing}},subset_img)
        wiseav_missingvals[wiseav_missingvals.==0.0].=missing

        wiseav_denom = size(wiseav_missingvals)[3].-sum(ismissing.(wiseav_missingvals),dims=3)

        avg_im[:,:,band] = sum(subset_img,dims=3)./wiseav_denom
        avg_im = convert(Array{Float32},avg_im)
        #println("$(avg_im[20,20,band])...$band")
    end

    avg_im = cat(input_image[:,:,1:Int((box_size-1)/2)],avg_im,dims=3)
    avg_im = cat(avg_im,input_image[:,:,end-Int((box_size-1)/2)+1:end],dims=3)

    avg_λvector = input_λvector[split_index+1:size(input_image)[3]-split_index]

    return (avg_im,avg_λvector)
end