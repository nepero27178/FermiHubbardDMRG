#!/usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../FermiHubbardDMRG/src

# Include setup
include(PROJECT_ROOT * "/setup/graphic-setup.jl")
include(PROJECT_ROOT * "/setup/simulations-setup.jl")

# Include modules
include(PROJECT_ROOT * "/modules/dmrg-routine.jl")
include(PROJECT_ROOT * "/modules/methods-plotting.jl")
include(PROJECT_ROOT * "/modules/methods-sweeping.jl")

# Move back
PROJECT_ROOT *= "/.."	# Absloute path up to .../FermiHubbardDMRG/

function Test()

	LL = [x for x in 38:16:150]
	VV = [5, 10, 20, 50]
		
	nSweeps = 5
	MaxLinkDim = [10,50,75,200,500]
	Cutoff = [1E-8]
	DMRGParameters = [IAFnSweeps, IAFMaxLinkDim, IAFCutoff]

	# LaTeX format
	Header = "& "
	for V in VV[1:end-1]
		Header *= "\$ $V \$ & "
	end
	Header *= "\$ $(VV[end]) \$ \\\\"
	println(Header)
	println("\\midrule")

	for (l,L) in enumerate(LL)
	
		Projections = zeros(length(VV))
		for (v,V) in enumerate(VV)
		
			N = Int64(L/2)
			E, psi = RunDMRGAlgorithm(
					[L, N, 1.0, V, V, 0.0],
					DMRGParameters,
					"Fast";
					verbose=false,
					FixedN=false,
					FixedParity=false
				)
				
			sites = siteinds(psi)

			HP = GetHalfMIProjector(sites)
			Projection = inner(psi', HP, psi)
			Projections[v] = Projection
		end
	
		# LaTeX format
		Line = "\$ $L \$ & "	
		for p in Projections[1:end-1]
			Line *= "\$ $(round(p,digits=3)) \$ & "
		end
		Line *= "\$ $(round(Projections[end],digits=3)) \$ \\\\"
	
		println(Line)

	end
end
