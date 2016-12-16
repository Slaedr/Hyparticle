include("burgers.jl")

function gridConvergence(Nstart,nmesh,xstart,xend,bvalue,amplitude,finaltime,maxtimesteps,ngraphplot)
	imesh = 1
	n = Nstart
	xsols = []
	vsols = []
	for imesh = 0:nmesh-1
		n = imesh*2^(2*imesh)
		xso,uso=burgersLoop(n,xstart,xend,bvalue,amplitude,finaltime,maxtimesteps,ngraphplot)
		push!(xsols,xso)
		push!(vsols,uso)
	end
end

# main

N = 20
xstart = 0.0
xend = 4*pi
bvalue = 1.5
amplitude = 0.5
finaltime = 0.5
maxtimesteps = 1000
ngraphplot = 200

println("Hyparticle for Burgers' equation - Initial data:")
println("  N = ", N, ", h = ", (xend-xstart)/N, ", final time = ", finaltime)

xso,uso=burgersLoop(N,xstart,xend,bvalue,amplitude,finaltime,maxtimesteps,ngraphplot)
#println(Array(xso)); println(uso)

plot(xso,uso,"o-")
xlabel("x")
ylabel("u")
grid("on")
show()
println()
