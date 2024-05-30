#Testing struct
using L2Pipeline
using GLMakie

caldata = caldata_from_url("C:/Lunar_Imagery_Data/020644_caldata")

#Running the pipeline
rem_solspec(caldata)

f = Figure()
ax = Axis(f[1,1])
image!(caldata.rdn[:,:,10])
display(GLMakie.Screen(),f)