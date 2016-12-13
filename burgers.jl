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
		cur.dt = -(nxt.x-p.dmin/2.0-cur.x)/(nxt.u-cur.u)
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
	prev = p.first
	cur = p.first.next
	for i = 2:p.n-2
		nxt = cur.next
		dist = nxt.x - cur.x
		if dist < p.dmin
			# remove the two particles and insert a new one
			newp = Particle()
			# compute the new particle's converved variable
			x1 = prev.x; x2 = cur.x; x23 = (nxt.x+cur.x)/2.0; x3 = nxt.x; x4 = nxt.next.x
			u1 = prev.u; u2 = cur.u; u3 = nxt.u; u4 = nxt.next.u
			tworexp = (x2-x1)*(u2+u1) + (x3-x2)*(u3+u2) + (x4-x3)*(u4+u3)
			newp.u = (tworexp - (x23-x1)*u1 - (x4-x23)*u4)/(x4-x1)
			newp.x = x23
			newp.v = newp.u
			# add the new particle to the list
			newp.next = nxt.next
			prev.next = newp
			p.n -= 1
			cur = newp
		elseif dist > p.dmax
			# insert a new particle, update n
			newp = Particle()
			x2 = cur.x; x3 = nxt.x; x23 = (nxt.x+cur.x)/2.0
			u2 = cur.u; u3 = nxt.u
			twor = (x3-x2)*(u3+u2)
			newp.u = (twor - (x23-x2)*u2 - (x3-x23)*u3)/(x3-x2)
			newp.x = x23
			newp.v = newp.u

			newp.next = nxt
			cur.next = newp
			cur = nxt
			p.n += 1
		else
			# nothing to do, move on
			prev = cur
			cur = cur.next
		end
	end
end

"""
Apply Dirichlet BC at left. Also remove extra particles at right.

We simply insert an element exactly at the left boundary having the boundary value bvalue.
"""
function applyDirichletBC!(p::ParticleList, bvalue::real)
	if p.first.x - p.xstart > p.dmax
		newp = Particle()
		newp.x = xstart
		newp.u = bvalue
		newp.v = bvalue
		newp.next = p.first
		p.first = newp
		p.n += 1
	end

	# check for particles that have overshot the right boundary
	cur = p.first.next
	prev = p.first
	for i = 3:p.n
		cur = cur.next
		prev = prev.next
	end
	if cur.x > p.xend
		prev.next = nothing
		n -= 1
	end
end

"""
Sets a sinosoidal initial condition.
"""
function initsin(xarr)
	a = xarr[1]; b = xarr[end]
	wl = (b-a)/5.0
	ws = a + (b-a)/10.0
	ic = zeros(xarr)
	for i = 1:length(xarr)
		if xarr[i] >= ws && xarr[i] <= ws+wl
			ic[i] = sin(2*pi*wl*(xarr[i]-ws))
		end
	end
	return ic
end

"""
Main time-stepping loop for Burgers' equation.
"""
function burgersLoop(N, xstart, xend, ttime)

	plist = ParticleList(N, xstart, xend)

	# set initial conditions
	pos = linspace(xstart,xend,N)
	uval = initsin(pos)
	vval = zeros(uval); vval[:] = uval[:]
	initialize!(plist,pos,vval,uval)

	t = 0; step = 0
	while t < ttime
		burgersTimeSteps!(plist)
		moveParticles!(plist)
		burgersParticleManagement!(plist)
		if step % 20 == 0
			println("Step ", step, ", time = ", t)
		end
		t += plist.gdt
		step += 1
	end

	# get solution arrays
	xsol,usol = outputToArrays(plist)
end

plist = ParticleList(5,0.0,0.0)
xp = linspace(plist.xstart,plist.xend,plist.n)
vp = ones(plist.n)
up = zeros(plist.n)
initialize!(plist,xp,vp,up)
Particles.printList(plist)
