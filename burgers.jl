#= Implements particle based solution of the inviscid Burgers' equation
=#

include("particleset.jl")
using Particles
using PyPlot

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
		if abs(nxt.v-cur.v) > eps(Particles.real)
			cur.dt = -(nxt.x-p.dmin/2.0-cur.x)/(nxt.v-cur.v)
		else
			cur.dt = 1.0
		end
		#print(cur.dt, " ")
		if cur.dt > 0.0
			p.gdt = min(p.gdt,cur.dt)
		end
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
	i = 2

	# because of the way cur moves, we need a while-limit that's independent of p.n
	nlim = p.n-2

	while i <= nlim
		nxt = cur.next
		dist = nxt.x - cur.x
		if dist < p.dmin
			# remove the two particles and insert a new one

			newp = Particle()
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
			nlim -= 1 
			# note that nlim mirrors p.n for this, cur moves 1 node forward

			cur = newp

		elseif dist > p.dmax
			# insert a new particle
			
			newp = Particle()
			x2 = cur.x; x3 = nxt.x; x23 = (nxt.x+cur.x)/2.0
			u2 = cur.u; u3 = nxt.u
			twor = (x3-x2)*(u3+u2)
			newp.u = (twor - (x23-x2)*u2 - (x3-x23)*u3)/(x3-x2)
			newp.x = x23
			newp.v = newp.u

			newp.next = nxt
			cur.next = newp
			p.n += 1 
			# note that nlim is not updated as cur moves 2 nodes forward

			cur = nxt
			prev = newp

		else
			prev = cur
			cur = cur.next

		end
		i += 1
	end
end

"""
Apply Dirichlet BC at left. Also remove extra particles at right.

We simply insert an element exactly at the left boundary having the boundary value bvalue.
"""
function applyDirichletBC!(p::ParticleList, bvalue::Particles.real)
	# if the fisrt particle is too far from the left boundary, insert a new particle
	if p.first.x - p.xstart > p.dmax
		newp = Particle()
		newp.x = xstart
		newp.u = bvalue
		newp.v = newp.u
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
		prev.next = 0
		p.n -= 1
	end
end

"""
Sets a sinosoidal initial condition.
"""
function initsin(xarr,bval,amp)
	a = xarr[1]; b = xarr[end]
	wl = (b-a)/5.0
	ws = a + (b-a)/10.0
	ic = bval*ones(xarr)
	for i = 1:length(xarr)
		if xarr[i] >= ws && xarr[i] <= ws+wl
			ic[i] = bval+amp*sin(2*pi*1.0/wl*(xarr[i]-ws))
		end
	end
	return ic
end

"""
Main time-stepping loop for Burgers' equation.
"""
function burgersLoop(N, xstart, xend, bvalue, initamp, ttime, maxiter)

	plist = ParticleList(N, xstart, xend)

	# set initial conditions
	pos = linspace(xstart,xend,N)
	uval = initsin(pos,bvalue,initamp)
	vval = zeros(uval); vval[:] = uval[:]
	initialize!(plist,pos,vval,uval)

	t = 0; step = 0
	while t < ttime && step < maxiter
		burgersTimeSteps!(plist)
		moveParticles!(plist)
		burgersParticleManagement!(plist)
		applyDirichletBC!(plist,bvalue)
		if step % 10 == 0
			println("Step ", step, ", time = ", t, ", time step = ", plist.gdt, ", n = ", plist.n)
		end
		t += plist.gdt
		step += 1
	end
	println("Time loop exited. Steps = ", step, ", time = ", t)

	# get solution arrays
	xsol,usol = outputToArrays(plist)
	return (xsol,usol)
end

# main

N = 200
xstart = 0.0
xend = 4*pi
bvalue = 1.5
amplitude = 1.0
finaltime = 3.0
maxtimesteps = 1000

println("Hyparticle for Burgers' equation - Initial data:")
println("  N = ", N, ", h = ", (xend-xstart)/N, ", final time = ", finaltime)

xso,uso=burgersLoop(N,xstart,xend,bvalue,amplitude,finaltime,maxtimesteps)

plot(xso,uso,"o-")
xlabel("x")
ylabel("u")
grid("on")
show()
println()
