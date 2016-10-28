""" 
Module for particle-based method to solve hyperbolic PDEs.
@author Aditya Kashi
@date October 2016
"""

"""
Import for Hyparticle functionality
"""
module Particles

# Stuff to be exported to code that import this module
export HyparticleSimulation

"""
An immutable type representing the set of particles and their properties.
"""
immutable ParticleSet
	x::Array{Float64}
	v::Array{Float64}
end

immutable HyparticleSimulation
	pset::ParticleSet				# Set of particles in the simulation
	dt::Float32						# time step
end

end
