#!/usr/bin/julia

using DelimitedFiles

PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/../setup/graphic-setup.jl")
include(PROJECT_ROOT * "/dmrg-routine.jl")
include(PROJECT_ROOT * "/methods-plotting.jl")
PROJECT_ROOT *= "/../.."

function GetStateProperties(
		FilePathOut::String,
		ModelParameters::Vector{Float64},
		ConserveNumber::Bool
	)

	# ------------------------------- Simulation ------------------------------- 

	nsweep = 10
	maxlinkdim = [10,50,75,200,500]
	cutoff = [1E-8]
	DMRGParameters = [nsweep, maxlinkdim, cutoff]
	
	Observables = RunDMRGAlgorithm(
		ModelParameters,
		DMRGParameters,
		"Debug",
		ConserveNumber
	)
									
	E, nMean, nVariance, DensityFluctuations, LocalE, P, S = Observables

	# --------------------------------- Write ----------------------------------

	DataFile = open(FilePathOut, "a")
	write(DataFile, "$ConserveNumber; $E; $LocalE; $nMean; $nVariance; ",
		"$DensityFluctuations; $P; $S\n")
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
		
		if Compute=="y"
			DataFile = open(FilePathOut, "w")
			write(DataFile, "# FixedN; E; LocalE; nMean; nVariance; ",
				"DensityFluctuations; Populations; Entropy\n")
			close(DataFile)
		end
		
		DirPathOut = PROJECT_ROOT * "/analysis/states-properties/"
		mkpath(DirPathOut)
		
		# Run both for fixed and variable fermions number
		for ConserveNumber in [true, false]
	
			global Counter += 1
			println("Simulation ($Counter/4): XY=$(XY) and ",
				"FixedN=$(ConserveNumber).")
	
			if Compute=="y"
				Observables = GetStateProperties(
					FilePathOut,
					ModelParameters,
					ConserveNumber
				)
				_, _, _, _, Γ, P, S = Observables
			elseif Compute=="n"
				ObservablesData = readdlm(FilePathOut, ';', Any, comments=true)
				for Index in 1:2
					if ObservablesData[Index,1]==ConserveNumber
						Γ, P, S = ObservablesData[6:8]
					end
				end			
			end
	
#			PlotPopulations(
#				DirPathOut,
#				P,
#				ModelParameters,
#				XY,
#				ConserveNumber
#			)
#			PlotBipartiteEntropy(
#				DirPathOut,
#				S,
#				ModelParameters,
#				XY,
#				ConserveNumber
#			)
#			printstyled("Plots ready!\n", color=:green)
		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
