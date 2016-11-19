"""
Module for particle-based method to solve Burgers' equation.
@author Aditya Kashi
@date October 2016
"""

"""
Module for Hyparticle functionality for Burgers' equation
"""
module BurgersParticles

# Stuff to be exported to code that import this module
export real, ParticleSet, burgersTimeSteps!

""" The floating-point type to be used """
real = Float64

""" A particle type for constructing a linked list of particles."""
type Particle
	x::real
	v::real
	u::real
	dt::real
	next::Particle
end

type ParticleList
	first::Particle
	n::UInt32				# number of particles
	xstart::real			# left boundary of domain
	xend::real				# right boundary of domain
	dmax::real				# max distance between adjacent particles
	dmin::real				# min distance between adjacent particles
	hi::real				# initial sampling resolution

	function ParticleList(n, xstart, xend, dmaxfactor=4.0/3.0, dmin=eps(real))
		# set up initial list
		first1 = Particle(0.0,0.0,0.0,0.0,0)
		cur = first1
		for i = 1:n
			cur.next = Particle(0.0,0.0,0.0,0.0,0)
			cur = cur.next
		end
		cur = 0
		hi = (xend-xstart)/n
		dmax = h*dmaxfactor
		return new(first1,n,xstart,xend,dmax.dmin,hi)
	end
end

"""
A type representing the set of particles and their properties.
"""
type ParticleSet
	n::UInt32				# number of particles
	xstart::real			# left boundary of domain
	xend::real				# right boundary of domain
	x::Array{real}			# positions
	v::Array{real}			# velocities
	u::Array{real}			# values of conserved variable
	dt::Array{real}			# time before collision
	dmax::real				# max distance between adjacent particles
	dmin::real				# min distance between adjacent particles
	hi::real				# initial sampling resolution

	function ParticleSet(n, xstart, xend, dmaxfactor=4.0/3.0, dmin=eps(real))
		x = zeros(n)
		v = zeros(n)
		u = zeros(n)
		dt = zeros(n)
		hi = (xend-xstart)/n
		dmax = h*dmaxfactor
		return new(n,xstart,xend,x,v,u,dt,dmax,dmin,hi)
	end
end

"""
Computes the allowable time step for each particle based on current positions and speeds for Burgers' flux.

Returns the minimum time step.
Note: dt[i] > 0 if the ith particle is faster than the next particle.
"""
function burgersTimeSteps!(p::ParticleSet)
	# compute time it takes for particle i to get to within dmin distance of particle i+1
	mdt = 1.0
	for i = 1:n-1
		p.dt[i] = -(p.x[i+1]-p.dmin-p.x[i])/(p.u[i+1]-p.u[i])
		mdt = min(mdt,p.dt[i])
	end
	return mdt
end

end
