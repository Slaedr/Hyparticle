"""
Module for particle-based method to solve Burgers' equation.
@author Aditya Kashi
@date October 2016
"""


# Module for base Hyparticle data structures.
module Particles

# Stuff to be exported to code that import this module
export real, Particle, ParticleList, deleteParticleAfter, insertParticleAfter, initialize

""" The floating-point type to be used """
real = Float64

""" A particle type for constructing a linked list of particles."""
type Particle
	x::real
	v::real
	u::real
	dt::real
	next::Particle

	function Particle()
		return new(0.0,0.0,0.0,0.0)
	end
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
		# set up initial list; note that one extra particle is allocated for convenience
		first1 = Particle()
		cur = first1
		for i = 1:n
			cur.next = Particle()
			cur = cur.next
		end
		
		hi = (xend-xstart)/n
		dmax = hi*dmaxfactor
		return new(first1,n,xstart,xend,dmax,dmin,hi)
	end
end

""" Deletes the next particle after par, ie the one referred to by par.next. 
 Note that the memory is released only when the GC sees fit."""
function deleteParticleAfter(plist::ParticleList, par::Particle)
	nexttonext = par.next.next
	par.next = nexttonext
	n -= 1
end

""" Add a particle newpar to the list after particle par."""
function insertParticleAfter(plist::ParticleList, newpar::Particle, par::Particle)
	temp = par.next
	par.next = newpar
	newpar.next = temp
	n += 1
end

""" Initialize the list with, for instalce, an initial condition."""
function initialize(plist::ParticleList, pos, vel, uval)
	cur = plist.first
	for i = 1:plist.n
		cur.x = pos[i]
		cur.v = vel[i]
		cur.u = uval[i]
		cur = cur.next
	end
end

""" For debugging"""
function printList(plist::ParticleList)
	cur = plist.first
	for i = 1:plist.n
		println(cur.x, " ", cur.v, " ", cur.u, " ", cur.dt)
		cur = cur.next
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

end
