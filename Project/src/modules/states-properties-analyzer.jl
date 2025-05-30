#!/usr/bin/julia

using DelimitedFiles

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../setup/graphic-setup.jl")
include(PROJECT_ROOT * "/dmrg-routine.jl")
include(PROJECT_ROOT * "/methods-plotting.jl")
PROJECT_ROOT *= "/../.."

@doc raw"""
function GetStateProperties(
		FilePathOut::String,
		ModelParameters::Vector{Float64},
		ConserveNumber::Bool,
		ConserveParity::Bool
	)
	
Returns: E, nMean, nVariance, DensityFluctuations, LocalE, P, S

This functions computes, given some `ModelParameters` and a logical input about
conserving or not the particles number (`ConserveNumber`) and the particles 
number parity (`ConserveParity`), energy, mean density and variance density,
the density fluctuation matrix, the local energy per site, the single-particle
state populations and the bipartite entropy.
"""
function GetStateProperties(
		FilePathOut::String,
		ModelParameters::Vector{Float64},
		ConserveNumber::Bool,
		ConserveParity::Bool
	)

	# ------------------------------- Simulation ------------------------------- 

	nsweep = 10
	maxlinkdim = [10,50,75,200,500]
	cutoff = [1E-8]
	DMRGParameters = [nsweep, maxlinkdim, cutoff]
	
	Observables = RunDMRGAlgorithm(
		ModelParameters,
		DMRGParameters,
		"StateAnalyzer"; #UserMode
		verbose=true,
		FixedN=ConserveNumber,
		FixedParity=ConserveParity
	)
	
	
	# GO ON FROM HERE
	# Extract compressibility and charge stiffness of the state: to be understood
	# if to be done as in Eq 5.4.1 conserving the particles number or not
	
	for N in L/2:2:L-2
		# Get compressibility and charge stiffnesss
		@sync begin
			# Create tasks for each DMRG call
			task1 = @spawn RunDMRGAlgorithm(
				[L, N, 1.0, V, μ0, 0.0],		# Central
				DMRGParameters,
				"Fast"; # UserMode
				verbose=false,
				FixedN=true,
				FixedParity=true
			)
			task2 = @spawn RunDMRGAlgorithm(
				[L, N, 1.0, V, μ0, ε], 		# Rotate clockwise
				DMRGParameters,
				"Fast"; # UserMode
				verbose=false,
				FixedN=true,
				FixedParity=true
			)
			task3 = @spawn RunDMRGAlgorithm(
				[L, N, 1.0, V, μ0, -ε], 		# Rotate counter-clockwise
				DMRGParameters,
				"Fast"; # UserMode
				verbose=false,
				FixedN=true,
				FixedParity=true
			)
			task4 = @spawn RunDMRGAlgorithm(
				[L, N+2, 1.0, V, μ0, 0.0], 	# Add two particles
				DMRGParameters,
				"Fast"; # UserMode
				verbose=false,
				FixedN=true,
				FixedParity=true
			)
			task5 = @spawn RunDMRGAlgorithm(
				[L, N-2, 1.0, V, μ0, 0.0], 	# Remove two particles
				DMRGParameters,
				"Fast"; # UserMode
				verbose=false,
				FixedN=true,
				FixedParity=true
			)

			# Wait for all tasks to complete and collect results
			E, psi = fetch(task1)
			EClock, _ = fetch(task2)
			ECounterClock, _ = fetch(task3)
			EAdd, _ = fetch(task4)
			EPop, _ = fetch(task5)
			
			D = pi * L * (EClock+ECounterClock-2*E)/(4ε^2)
			k = 4/(L*(EAdd + EPop - 2*E))
		end
	end
	
	# Energy, local energy, state, local denisty, density-density correlator, entropy
	E, e, psi, n, Cnn, S, k, D = Observables

	# --------------------------------- Write ----------------------------------

	DataFile = open(FilePathOut, "a")
		write(DataFile, "$ConserveNumber; $E; $e; $n; $nC; $S\n")
	close(DataFile)
	
	return Observables
end

# ----------------------------------- Main -------------------------------------
	
function main()

	global Counter = 0
	L = 50
	N = 25
	η = 0.0
	
	print("Should I run the simulations? (y/n) ")
	Compute = readline()
	
	for Phase in ["XY", "IAF", "IF"]
	
		DirPathOut = PROJECT_ROOT * "/simulations/states-properties/"
		mkpath(DirPathOut)
		
		if Phase=="XY"			# XY state (Jordan-Wigner mapping)  
			V = 0.1
			μ = 0.0
			FilePathOut = DirPathOut * "XY_V=$(V)_μ=$(μ).txt"
		elseif Phase=="IAF"		# IAF state (Jordan-Wigner mapping), Mott insulating
			V = 2.5
			μ = 0.0
			FilePathOut = DirPathOut * "IAF_V=$(V)_μ=$(μ).txt"	
		elseif Phase=="IF"		# IF state (Jordan-Wigner mapping), Mott insulating
			V = -2.5
			μ = 0.0
			FilePathOut = DirPathOut * "IF_V=$(V)_μ=$(μ).txt"
		end
		
		@info "State parameters" L N V μ
		
		ModelParameters = [L, N, 1.0, V, μ, 0.0]
		ConserveNumber = true
		ConserveParity = true
		
		if Compute=="y"
			DataFile = open(FilePathOut, "w")
			write(DataFile, "# FixedN; E; LocalE; nMean; nVariance; ",
				"DensityFluctuations; Populations; Entropy\n")
			close(DataFile)
		end
		
		DirPathOut = PROJECT_ROOT * "/analysis/states-properties/"
		mkpath(DirPathOut)
		
		Observables = GetStateProperties(
			FilePathOut,
			ModelParameters,
			ConserveNumber,
			ConserveParity
		)
		E, e, psi, n, Cnn, S = Observables
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
