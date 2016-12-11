#= Implements particle based solution of the inviscid Burgers' equation
=#

include("particleset.jl")
using Particles

"""
Computes the allowable time step for each particle based on current positions and speeds for Burgers' flux.

Returns the minimum time step.
Note: dt[i] > 0 if the ith particle is faster than the next particle.
"""
#=function burgersTimeSteps!(p::ParticleSet)
	# compute time it takes for particle i to get to within dmin distance of particle i+1
	mdt = 1.0
	for i = 1:n-1
		p.dt[i] = -(p.x[i+1]-p.dmin-p.x[i])/(p.u[i+1]-p.u[i])
		mdt = min(mdt,p.dt[i])
	end
	return mdt
end=#

plist = ParticleList(5,0.0,0.0)
xp = linspace(plist.xstart,plist.xend,plist.n)
vp = ones(plist.n)
up = zeros(plist.n)
initialize(plist,xp,vp,up)
Particles.printList(plist)
