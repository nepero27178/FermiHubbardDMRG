#!/usr/bin/julia

using DelimitedFiles
using Base.Threads

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../setup/simulations-setup.jl")
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
	
	ε = 0.1     # Arbitrary
	Δμ = 0.5    # Arbitrary
	L, N, _, V, μ, _ = ModelParameters
	
	# ------------------------------- Simulation ------------------------------- 

	DMRGParameters = global DMRGParametersXY # Always use XY parameters
	
	Observables = RunDMRGAlgorithm(
		ModelParameters,
		DMRGParameters,
		"StateAnalyzer"; #UserMode
		verbose=false,
		FixedN=ConserveNumber,
		FixedParity=ConserveParity
	)
	
	# Energy, local energy, state, local denisty, density-density correlator, entropy
	E, lE, psi, n, Cnn, S = Observables
	ρ = sum(n)/L
	δn2M = GetBlockVariance(n, Cnn)z
	k = 0 # Initialize
	D = 0 # Initialize
	
	# Get compressibility and charge stiffnesss
	@sync begin
		# Create tasks for each DMRG call
		task1 = @spawn RunDMRGAlgorithm(
			[L, N, 1.0, V, μ+Δμ, 0.0],		# Central
			DMRGParameters,
			"Fast"; # UserMode
			verbose=false,
			FixedN=ConserveNumber,
			FixedParity=ConserveParity
		)
		task2 = @spawn RunDMRGAlgorithm(
			[L, N, 1.0, V, μ, ε], 			# Rotate clockwise
			DMRGParameters,
			"Fast"; # UserMode
			verbose=false,
			FixedN=ConserveNumber,
			FixedParity=ConserveParity
		)
		task3 = @spawn RunDMRGAlgorithm(
			[L, N, 1.0, V, μ, -ε], 		# Rotate counter-clockwise
			DMRGParameters,
			"Fast"; # UserMode
			verbose=false,
			FixedN=ConserveNumber,
			FixedParity=ConserveParity
		)

		# Wait for all tasks to complete and collect results
		EShift, psiShift = fetch(task1)
		EClock, _ = fetch(task2)
		ECounterClock, _ = fetch(task3)
		
        ρShift = sum(expect(psiShift, "n"))/L
        Δρ = ρShift - ρ
    
		D = pi * L * (EClock+ECounterClock-2*E)/(4ε^2)
		k = Δρ/Δμ
	end

	# --------------------------------- Write ----------------------------------

	DataFile = open(FilePathOut, "a")
		write(DataFile, "$ConserveNumber; $ConserveParity; $L; $N; $V; $μ; $E; $lE; $n; $(δn2M); $S; $k; $D\n")
	close(DataFile)
	
	return E, psi, ρ, δn2M, S, k, D
end

# ----------------------------------- Main -------------------------------------
	
function main()

	L = global StatePropertiesL
	η = 0.0
	
	print("Should I run the simulations? (y/n) ")
	Compute = readline()
	
	for Phase in ["XY", "IAF", "IF"]		
	
		printstyled("Simulating phase: $(Phase)", color=:yellow)

		DirPathOut = PROJECT_ROOT * "/simulations/states-properties/"
		mkpath(DirPathOut)
		
		if Phase=="XY"			# XY state (Jordan-Wigner mapping)  
			V = 0.5
			μ = 1.0
			FilePathOut = DirPathOut * "XY_V=$(V)_μ=$(μ).txt"
		elseif Phase=="IAF"		# IAF state (Jordan-Wigner mapping), Mott insulating
			V = 5.0
			μ = 4.0
			FilePathOut = DirPathOut * "IAF_V=$(V)_μ=$(μ).txt"	
		elseif Phase=="IF"		# IF state (Jordan-Wigner mapping), Mott insulating
			V = -2.5
			μ = 1.0
			FilePathOut = DirPathOut * "IF_V=$(V)_μ=$(μ).txt"
		end

	        if Compute=="y"
			DataFile = open(FilePathOut, "w")
				write(DataFile, "# FixedN; FixedParity; L; N; V; μ; E; e; n;"
					* " δn_M^2; S; k; D [calculated at $(now())]\n")
			close(DataFile)
		end

		for ConserveNumber in [true, false]

			if ConserveNumber
				# If particles number is conserved also parity is.
				BoolArray = [true]
			elseif !ConserveNumber
				BoolArray = [true, false]
			end

			for ConserveParity in BoolArray
				
				println()
				@info "State parameters" ConserveNumber ConserveParity L N V μ
				ModelParameters = [L, N, 1.0, V, μ, 0.0]
				
				DirPathOut = PROJECT_ROOT * "/analysis/states-properties/"
				mkpath(DirPathOut)
				
				Observables = GetStateProperties(
					FilePathOut,
					ModelParameters,
					ConserveNumber,
					ConserveParity
				)
				E, lE, ρ, δn2M, S, k, D = Observables
				
				@info Phase ConserveNumber ConserveParity ρ k D
				println()
			end

		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
