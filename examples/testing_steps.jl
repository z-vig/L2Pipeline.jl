#Testing struct
using L2Pipeline
using GLMakie

caldata = caldata_from_url("C:/Lunar_Imagery_Data/020644_caldata")

#Running the pipeline
rem_solspec(caldata,"C:/Lunar_Imagery_Data/020644_caldata")

xtest = rand(axes(caldata.rdn,1))
ytest = rand(axes(caldata.rdn,2))

f = Figure()
ax1 = Axis(f[1,1],ylabel="Radiance")
ax2 = Axis(f[2,1],ylabel="I/F Corrrection")
ax3 = Axis(f[:,2])

lines!(ax2,test[xtest,ytest,:],label="SolSpec Removed")
lines!(ax1,caldata.rdn[xtest,ytest,:],label="Original")
image!(ax3,caldata.rdn[:,:,240])
scatter!(ax3,xtest,ytest)
hidedecorations!(ax3)

Label(f[0,:],"Displaying: ($(xtest),$(ytest))")

display(GLMakie.Screen(),f)