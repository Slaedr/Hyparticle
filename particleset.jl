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
An type representing the set of particles and their properties.
"""
type ParticleSet
	n::UInt32				# number of particles
	x::Array{Float64}		# positions
	v::Array{Float64}		# velocities
	u::Array{Float64}		# values of conserved variable
	dt::Array{Float64}		# time before collision
end

type HyparticleSimulation
	pset::ParticleSet					# Set of particles in the simulation
	dts::Float64						# time step
end

end
