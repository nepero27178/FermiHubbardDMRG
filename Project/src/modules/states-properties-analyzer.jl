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
	L, N = Int64.(ModelParameters[1:2])
    _, V, μ, _ = ModelParameters[3:6]
	
	# ------------------------------- Simulation ------------------------------- 

#	DMRGParameters = global DMRGParametersXY # Always use XY parameters
	
	Observables = RunDMRGAlgorithm(
		ModelParameters,
		DMRGParametersXY,	# Overkill
		"StateAnalyzer"; 	# UserMode
		verbose=false,
		FixedN=ConserveNumber,
		FixedParity=ConserveParity
	)
	
	# Energy, local energy, state, local density, density-density correlator, entropy
	E, lE, psi, n, CnnMatrix, S = Observables
	
	ρ = sum(n)/L
	δn2M = GetBlockVariance(n, CnnMatrix)
	
	Cnn = zeros(Int64(L/2)+1)
	Csu = zeros(Int64(L/2)+1)	
	CsuMatrix = GetSuperconductingCorrelator(psi)
#	eCsu = std(Csu, dims=2)[:] # [:] Matrix{Float64} -> Vector{Float64}
#	Csu = mean(Csu, dims=2)[:] # [:] Matrix{Float64} -> Vector{Float64}
	for r in 0:Int64(L/2)
		TmpCnn = 0
		TmpCsu = 0
		for j in 1:L
			TmpCnn += CnnMatrix[j, mod1(j+r, L)]
			TmpCsu += CsuMatrix[j, mod1(j+r, L)]
		end
		Cnn[r+1] = TmpCnn/L
		Csu[r+1] = TmpCsu/L
	end
	
#	eCsu = std(Csu, dims=2)[:] # [:] Matrix{Float64} -> Vector{Float64}
#	Csu = mean(Csu, dims=2)[:] # [:] Matrix{Float64} -> Vector{Float64}
	
	k = 0 # Initialize
	D = 0 # Initialize
	
	# Get compressibility and charge stiffnesss
	@sync begin
		# Create tasks for each DMRG call
		task1 = @spawn RunDMRGAlgorithm(
			[L, N, 1.0, V, μ+Δμ, 0.0],		# Central
			DMRGParametersXY,				# Overkill
			"Fast"; # UserMode
			verbose=false,
			FixedN=ConserveNumber,
			FixedParity=ConserveParity
		)
		task2 = @spawn RunDMRGAlgorithm(
			[L, N, 1.0, V, μ, ε], 			# Rotate clockwise
			DMRGParametersXY,				# Overkill
			"Fast"; # UserMode
			verbose=false,
			FixedN=ConserveNumber,
			FixedParity=ConserveParity
		)
		task3 = @spawn RunDMRGAlgorithm(
			[L, N, 1.0, V, μ, -ε],	 		# Rotate counter-clockwise
			DMRGParametersXY,				# Overkill
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

	Bools = "$ConserveNumber; $ConserveParity; "
	Macro = "$L; $N; $V; $μ; $E; "							# Macroscopic
	Chain = "$(lE); $(n); $(S); $(δn2M); $(Cnn); $(Csu); "	# Microscopic
	Bos = "$k; $D"											# Bosonization
	
	DataFile = open(FilePathOut, "a")
		write(DataFile, Bools * Macro * Chain * Bos * "\n")
	close(DataFile)
	
	return E, psi, ρ, S, δn2M, Cnn, Csu, k, D
end

# ----------------------------------- Main -------------------------------------
	
function main()

	L = StatePropertiesL
	N = StatePropertiesN
#	global IFPoint
#	global IAFPoint
#	global XYPoint1
#	global XYPoint2
	
	PhasesDict = Dict([
		("XY-Up", XYPoint1),	# Global
		("XY-Down", XYPoint2),	# Global
		("IF", IFPoint),		# Global
		("IAF", IAFPoint)		# Global
	])
	
	η = 0.0
	
	print("Should I run the simulations? (y/n) ")
	Compute = readline()
	
	for Phase in ["XY-Up", "XY-Down", "IAF", "IF"]		
	
		printstyled("Simulating phase: $(Phase)", color=:yellow)

		DirPathOut = PROJECT_ROOT * "/simulations/states-properties/"
		mkpath(DirPathOut)

		Point = PhasesDict[Phase]
		V = Point[1]
		μ = Point[2]
		FilePathOut = DirPathOut * Phase * "_V=$(V)_μ=$(μ)_L=$(L).txt"
		
	    if Compute=="y"
			DataFile = open(FilePathOut, "w")
				write(DataFile, "# FixedN; FixedParity; L; N; V; μ; E; e; n;"
					* " S; δn_M^2; <Cnn>; <Csu>; k; D; [calculated at $(now())]\n")
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
				E, lE, ρ, S, δn2M, Cnn, Csu, k, D = Observables
				
				@info Phase ConserveNumber ConserveParity ρ k D
				println()
			end

		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
