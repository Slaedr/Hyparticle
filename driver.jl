include("Burgers.jl")
using Burgers
using PyPlot

function gridConvergence(Nstart,nmesh,Nref,xstart,xend,bvalue,amplitude,finaltime,maxtimesteps,ngraphplot)
	imesh = 1
	n = Nstart
	vsols = []
	
	# reference solution
	println("Computing reference solution")
	xsols,uref=burgersLoop(Nref,xstart,xend,bvalue,amplitude,finaltime,maxtimesteps,ngraphplot)
	println("Reference solution computed")
	errs::Array{Particles.real} = []
	hs::Array{Particles.real} = []

	for imesh = 0:nmesh-1
		n = n*2
		xso,uso,nfin = burgersLoop(n,xstart,xend,bvalue,amplitude,finaltime,maxtimesteps,ngraphplot)
		push!(vsols,uso)
		push!(hs, (xend-xstart)/(nfin-1))
		push!(errs, sum(abs(uso .- uref))*hs[end])
	end
	hsl, errsl = log10(hs), log10(errs)
	plot(hsl,errsl,"o-")
	plot(hsl, 2.*(hsl.-hsl[1]).+errsl[1], "--",label="Slope 2.0")
	legend(loc="best")
	grid("on")
	finslope = (errsl[1:end-1].-errsl[2:end])./(hsl[1:end-1].-hsl[2:end])
	println("Slope = ", finslope)
	show()
end

# main

nstart = 20
nref = 20000
nmesh = 8
ngraphplot = 10000
xstart = 0.0
xend = 4*pi
bvalue = 1.5
amplitude = 0.5
finaltime = 0.1
maxtimesteps = 1000000

println("Hyparticle for Burgers' equation - Initial data:")
println("  Nstart = ", nstart, ", Nref = ", nref, ", final time = ", finaltime)

#gridConvergence(nstart,nmesh,nref,xstart,xend,bvalue,amplitude,finaltime,maxtimesteps,ngraphplot)

xso,uso = burgersLoop(200, xstart,xend,bvalue,amplitude,finaltime,maxtimesteps,200)

plot(xso,uso,"-")
xlabel("x")
ylabel("u")
grid("on")
show()
println()
