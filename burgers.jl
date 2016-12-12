#= Implements particle based solution of the inviscid Burgers' equation
=#

include("particleset.jl")
using Particles

"""
Computes the allowable time step for each particle based on current positions and speeds for Burgers' flux.

Returns the minimum time step.
Note: dt[i] > 0 if the ith particle is faster than the next particle.
"""
function burgersTimeSteps!(p::ParticleList)
	# compute time it takes for particle i to get to within dmin distance of particle i+1
	cur = p.first
	p.gdt = 1.0
	for i = 1:p.n-1
		nxt = cur.next
		cur.dt = -(nxt.x-p.dmin-cur.x)/(nxt.u-cur.u)
		p.gdt = min(p.gdt,cur.dt)
		cur = cur.next
	end
end

"""
One step of particle management. Note that boundaries are not treated.

When distance between adjacent particles is greater than dmax, a new particle is inserted.
when distance between adjacent particles is less than dmin, the two particles are deleted and a new one is inserted.
Value of the inserted particle is computed from local conservation.
"""
function burgersParticleManagement!(p::ParticleList)
	cur = p.first.next
	for i = 2:n-2
		nxt = cur.next
		dist = nxt.x - cur.x
		if dist < p.dmin
			# TODO: remove the two particles and insert a new one, update n
		elseif dist > dmax
			# TODO: insert a new particle, update n
		end
		cur = cur.next
	end
	# TODO: handle boundaries
end

"""
Move particles using computed time step.
"""
function moveParticles!(p::ParticleList)
	cur = p.first
	for i = 1:n
		cur.x += cur.v*p.gdt
		cur = cur.next
	end
end

plist = ParticleList(5,0.0,0.0)
xp = linspace(plist.xstart,plist.xend,plist.n)
vp = ones(plist.n)
up = zeros(plist.n)
initialize(plist,xp,vp,up)
Particles.printList(plist)
