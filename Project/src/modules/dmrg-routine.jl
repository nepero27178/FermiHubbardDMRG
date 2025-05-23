#!/usr/bin/julia

using ITensors, ITensorMPS
using Statistics
using Dates
using DelimitedFiles

include("physics-definitions.jl")

# ------------------------------------------------------------------------------
# -------------------------------- Initializers -------------------------------- 
# ------------------------------------------------------------------------------

@doc raw"""
function SetUniformState(
		sites::Any
	)::MPS
	
Returns: uniformly filled fermionic state.	
	
Create an initial MPS with `N` particles on the fermionic chain; in the spinless
system, the state is at unitary filling; in the spinful system, it is at unitary
filling.
"""
function SetUniformState(	# Ferromagnetic state
		sites::Any
	)::MPS

    L = length(sites)
    states = ["1" for _ in 1:L]

    return MPS(sites, states)
end

@doc raw"""
function SetStartingState(
		sites::Any,
        N::Int64;
        ForceCenter=false
	)::MPS
	
Returns: starting fermionic state, filled only at center.	
	
Create an initial MPS with `N` particles on the fermionic chain, with `N≤L`. Two
modes are implemented:
- If `N≠L/2`, all the fermions are created at the center of the chain. 
For periodic and twisted boundary conditions this is irrelevant, but can become 
useful for open boundary conditions.
- If `N=L/2`, the fermions are alternated (`1`/`0`/`1`/`0`/...), except the 
parameter `ForceCenter` is set to `true`, in which case the above condition
applies.
"""
function SetStartingState(
		sites::Any,
		N::Int64;
        ForceCenter=false
	)::MPS

	L = length(sites)
	if N>L || N<1
		error("Invalid filling! Enter 1≤N≤L.")
	end 
	
    if N==L/2 && !ForceCenter
        states = ["0" for _ in 1:L]
        for j in 1:N
            states[2*j-1] = "1"
        end
    else
        states = vcat(["1" for _ in 1:N], ["0" for _ in N+1:L])
	    circshift!(states, floor(Int64,(L-N)/2))
    end

    return MPS(sites, states)
end

# ------------------------------------------------------------------------------
# ------------------------------------ DMRG ------------------------------------ 
# ------------------------------------------------------------------------------

@doc raw"""
function RunDMRGAlgorithm(
		ModelParameters::Vector{Float64},
		DMRGParameters::Vector{Any},
		UserMode::String,
        FixedN::Bool;
		verbose=false,
	)::Any

Returns: results of DMRG optimization and relevant observables as Float(s)64

Run DMRG algorithm with chosen parameters, and return results. The on-site
interaction is absent in the fermionic polarized case.

Input:

    - ModelParameters: array of [L::Int64, N::Int64, t::Float64, V::Float64, μ::Float64, η::Float64]
    - DMRGParameters: array of [nsweeps::Int64, maxdim::Int64, cutoff::Vector{Float64}]
    							
Parametric input:

    - ComputeAllObservables: boolean variable, if false only E, nVariance and Γ are computed (default: false)
    - ComputeGamma: boolean variable, if false all positional correalators are not extracted (default: false)
    - ComputeC: boolean variable, if false all positional correalators are not extracted (default: false)
"""
function RunDMRGAlgorithm(
		ModelParameters::Vector{Float64},		# Model parameters
		DMRGParameters::Vector{Any},			# Parameters for DMRG
		UserMode::String,						# Parametric input
		FixedN::Bool;							# Let N vary
		verbose=false							# Do not print at line
	)::Any
    
    L, N = Int64.(ModelParameters[1:2])
    t, V, μ, η = ModelParameters[3:6]
    nsweeps, maxdim, cutoff = DMRGParameters[1:3]

	if η>0.5 || η<-0.5
    	error("Invalid phase! Enter -1/2≤η≤1/2.")
	end 

	# Evaluate UserMode
	ModeErrorMsg = "Input error: use as argument \"OrderParameters\", \"Correlator\", \"Debug\" or \"Fast\":
- \"OrderParameters\": return E, nVariance, aAvg;
- \"Correlator\": return E, Γ, eΓ;
- \"Debug\": return E, LocalE, nMean, nVariance;
- \"Fast\": return E."
	
	OrderParameters=false
	Correlators=false
	Debug=false
	Fast=false
	
	if UserMode=="OrderParameters"
		
		@warn "Mode OrderParameters under construction."
	
#		OrderParameters=true
#		nVariance = 0		 					# Variance on central site
#   	aAvg = 0		 						# <a> on central site
	
	elseif UserMode=="Correlators"
		
		Correlators=true
		
	elseif UserMode=="Debug"
	
		Debug=true
		nMean = zeros(L)				# Mean number of particles per site
		nVariance = zeros(L) 			# Variance on n_i, for all sites i
    	LocalE = zeros(L)				# Local contribution to the energy
    	Populations = zeros(L,2)	 	# Single site populations
    	Entropy = zeros(L-1)			# Bipartite entropy
    
    elseif UserMode == "Fast"
    	Fast = true
    		
	else
		error(ModeErrorMsg)
	end

    if verbose
    	if !Fast
    		println("Selected mode: ", UserMode)
        end
        @info "Model parameters" L N t V μ η
        printstyled("Starting simulation...\n", color=:yellow)
    end

	# Initialize lattice 
    sites = siteinds("Fermion", L, conserve_nf=FixedN)
    
    # Compute hamiltonian (two-cases separation is needed in order to simplify
    # initialization).
    if η==0
		H = GetHamiltonianMPO(sites, t, V, μ)
	elseif η!==0 								
	    H = GetHamiltonianMPO(sites, t, V, μ; Φ=2*pi*η)
	end
    
    # Set starting state
    psi0 = SetStartingState(sites, N)
	# Set Ferromagnetic state
	# psi0 = SetUniformState(sites)

    if verbose
        # Print observables   
        N0 = GetTotalFermionNumber(psi0)
        E0 = inner(psi0', H, psi0)
        @info "Expectation values on the initial state" N0 E0
        ShowDensityProfile(psi0)
    end

    # Run DMRG algorithm and print results
    E, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, outputlevel=verbose)
    FinalAmplitude = dot(psi0, psi)
    @info "Superposition of final state with the initializer" FinalAmplitude

    # Sanity checks: calculate whether Ntot has been conserved, and the found 
    # ground state is actually an eigenstate of H.
    if verbose
    	ParticlesNumber = GetTotalFermionNumber(psi)
        VarE = inner(H, psi, H, psi) - E^2
		EnergyRelativeError = sqrt(VarE)/E
		@info "Expectation values on the final state" ParticlesNumber EnergyRelativeError 
		ShowDensityProfile(psi)
    end

    if OrderParameters

		@warn "Mode OrderParameters under construction."
#    	Index = ceil(Int64, L/2)
#    	
#		if !verbose
#			NTotAvg = inner(psi', NTot, psi)
#		end    	
#		
#		nVariance = GetNumberVariance(psi, sites, Index)
#	    aAvg = expect(psi, "a"; sites=Index)
#		
#		return E, NTotAvg, nVariance, aAvg
    
    elseif Correlators
    
    	# TODO Extend to Green's function

		l = floor(Int64, L/2)		# Correlations up to half the closed chain
		DensityCorrelator = GetDensityCorrelator(psi)
		Γ = zeros(l)
		eΓ = zeros(l)
		
		SuperconductingCorrelator = GetSuperconductingCorrelator(psi)
		O = zeros(l)
		eO = zeros(l)
		
		for r in 1:l
		
			TmpDensity = zeros(L)
			TmpSC = zeros(L)
			for i in 1:L
				TmpDensity[i] = DensityCorrelator[i,mod1(i+r,L)]
				TmpSC[i] = SuperconductingCorrelator[i,mod1(i+r,L)]
			end
			
			Γ[r] = mean(TmpDensity)
			eΓ[r] = std(TmpDensity)
			
			O[r] = mean(TmpSC)
			eO[r] = std(TmpSC)
			
		end
		
		return E, Γ, eΓ, O, eO
    
    elseif Debug
		
		DensityCorrelator = GetDensityCorrelator(psi)
		R = floor(Int64, L/2)	# Max distance
		DensityFluctuations = zeros(R)
		for r in 1:R
			Tmp = 0
			for j in 1:L
				Tmp += DensityCorrelator[j,mod1(j+r,L)]
			end
			Tmp /= R
			DensityFluctuations[r] = Tmp
		end
		
	    for i in 1:L
	    	# Locally defined
	    	nMean[i] = expect(psi, "n"; sites=i)
			nVariance[i] = DensityCorrelator[i,i]
			LocalE[i] = inner(psi', GetLocalHamiltonianMPO(sites, i, t, V, μ), psi)
    	end
    	Populations = GetLocalPopulation(psi)
    	
    	for b in 1:L-1
    		# Compute entropy for all possible bipartitions
    		Entropy[b] = GetVonNeumannEntropy(psi, b)
    	end
    	
    	ShowDensityProfile(psi)    	
    	return E, nMean, nVariance, DensityFluctuations, LocalE, Populations, Entropy
    	
    elseif Fast
    	return E
    end
end

# ------------------------------------------------------------------------------
# ------------------------------------ Main ------------------------------------ 
# ------------------------------------------------------------------------------

@doc raw"""
function main()

If the script is called directly from terminal, run one full DMRG routine
with the parameters defined inline.
"""
function main()
    
    L = 15
    N = 5
   	t = 1.0 # Keep it fixed
   	V = 1.0
    μ = 1.0
    η = 0

    nsweep = 10
    maxlinkdim = [10,50,75,200,500]
    cutoff = [1E-8]
    DMRGParameters = [nsweep, maxlinkdim, cutoff]
    
    # Choose one: "OrderParameters" / "Correlators" / "Debug" / "Fast"
    UserMode = "Correlators"

	if UserMode=="Fast"
		for i in 1:80
    		print("\r")
    		for j in 1:i
    			printstyled("-", color=:cyan)
    		end
    		printstyled(" >>> Speedin' >>>", color=:cyan)
    		for _ in 1:80-i
    			printstyled("-", color=:cyan)
    		end
    		sleep(0.01)
    	end
    	print("\n")
	end
	
    # ------------------------ Results for a single run ------------------------
    
    Observables = RunDMRGAlgorithm(
		[L, N, t, V, μ, η],
		DMRGParameters,
		UserMode,
		FixedN; 
		verbose=true,
	)
    								
	if UserMode=="OrderParameters"
		
		@warn "Mode OrderParameters under construction."
	
#		E, nVariance, aAvg = Observables
#		println("Results of the simulation:
#Energy of ground state: $(round.(E, digits=4))
#Number variance on central site: $(round.(nVariance, digits=4))
#Average of <a> on central site: $(round.(aAvg, digits=4))")
	
	elseif UserMode=="Correlators"
		E, Γ, eΓ, O, eO = Observables
		M = vcat(["r" "Γ" "eΓ" "O" "eO"], hcat([x for x in 1:floor(Int64, L/2)], [Γ eΓ O eO]))
		@info "Results of the simulation:" E M
		
	elseif UserMode=="Debug"
	
		E, nMean, nVariance, LocalE, Populations, Entropy = Observables
		println("Results of the simulation:
Energy of ground state: $(round.(E, digits=4))
\"Local\" part of energy: $(round.(LocalE, digits=4))
Mean number of particles: $(round.(nMean, digits=4))
Variance number of particles: $(round.(nVariance, digits=4))
Relative fluctuation: $(round.(sqrt.(nVariance)./nMean, digits=4))
Populations: $(round.(Populations, digits=4))
Bipartite entropy: $(round.(Entropy, digits=4))")
    	
	else
		E = Observables
		println("Results of the simulation:
Energy of ground state: $(round.(E, digits=4))")

		for i in 1:80
    		print("\r")
    		for _ in 1:80-i
    			printstyled("-", color=:cyan)
    		end
    		printstyled("<<< I'm speed <<< ", color=:cyan)
    		for _ in 1:i
				printstyled("-", color=:cyan)
    		end
    		sleep(0.01)
    	end
    	print("\n")
	end

end

if abspath(PROGRAM_FILE) == @__FILE__ # equivalent to if __name__ == "__main__"
	
	"""
	If this file is directly compiled, run a single DMRG simulation by the
    parameters defined in function main().
    """
    
    main()
end
