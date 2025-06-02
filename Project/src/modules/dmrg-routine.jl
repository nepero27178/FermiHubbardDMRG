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
		ModelParameters::Vector{Float64},		# Model parameters
		DMRGParameters::Vector{Any},			# Parameters for DMRG
		UserMode::String;						# Parametric input
		verbose=false,							# Do not print at line
		FixedN=false							# Let N vary
		FixedParity=false						# Let parity vary
	)::Any

Returns: results of DMRG optimization and relevant observables as Float(s)64

Run DMRG algorithm with chosen parameters, and return results. The on-site
interaction is absent in the fermionic polarized case.

Input:

    - ModelParameters: array of [L::Int64, N::Int64, t::Float64, V::Float64, μ::Float64, η::Float64]
    - DMRGParameters: array of [nsweeps::Int64, maxdim::Int64, cutoff::Vector{Float64}]
    - UserMode: string from the set [\"Debug\", \"Fast\", ...]
    							
Parametric input:

    - verbose=false (print at line)
    - FixedN=false (fix particles number)
    - FixedParity=false (fix particles number parity)
"""
function RunDMRGAlgorithm(
		ModelParameters::Vector{Float64},		# Model parameters
		DMRGParameters::Vector{Any},			# Parameters for DMRG
		UserMode::String;						# Parametric input
		verbose=false,							# Do not print at line
		FixedN=false,							# Let N vary
		FixedParity=false						# Let parity vary
	)::Any
    
    L, N = Int64.(ModelParameters[1:2])
    t, V, μ, η = ModelParameters[3:6]
    nsweeps, maxdim, cutoff = DMRGParameters[1:3]

	if η>0.5 || η<-0.5
    	error("Invalid phase! Enter -1/2≤η≤1/2.")
	end 

	# Evaluate UserMode
	ModeErrorMsg = "Input error: use as argument \"OrderParameters\", \"Correlator\", \"Fast\", \"StateAnalyzer\"  or \"Debug\":
- \"OrderParameters\": return E, nVariance, aAvg;
- \"Correlator\": return E, Γ, eΓ;
- \"Fast\": return E;
- \"StateAnalyzer\": return [...]
- \"Debug\": return E, LocalE, nMean, nVariance;"
	
	OrderParameters=false
	Correlators=false
	Fast=false
	StateAnalyzer=false
	Debug=false
	
	if UserMode=="OrderParameters"
		
		@warn "Mode OrderParameters under construction."
	
#		OrderParameters=true
#		nVariance = 0		 					# Variance on central site
#   	aAvg = 0		 						# <a> on central site
	
	elseif UserMode=="Correlators"
		Correlators=true
    
    elseif UserMode == "Fast"
    	Fast = true
    	
    elseif UserMode == "StateAnalyzer"
    	StateAnalyzer=true
    	SiteDensity = zeros(L)			# Mean number of particles per site
    	LocalE = zeros(L)				# Local contribution to the energy
#    	Populations = zeros(L,2)	 	# Single site populations
    	Entropy = zeros(L-1)			# Bipartite entropy (extend to L sites)
    	
    elseif UserMode=="Debug"

		@warn "Mode Debug under construction."
		Debug=true
    		
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
    sites = siteinds("Fermion", L, conserve_nf=FixedN, conserve_nfparity=FixedParity)
    
    # Compute hamiltonian (two-cases separation is needed in order to simplify
    # initialization).
    if η==0
		H = GetHamiltonianMPO(sites, t, V, μ)
	elseif η!==0 								
	    H = GetHamiltonianMPO(sites, t, V, μ; Φ=2*pi*η)
	end
    
    #TODO Finalize state initialization phase-wise
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

    # Sanity checks: calculate whether Ntot has been conserved, and the found 
    # ground state is actually an eigenstate of H.
    if verbose
		FinalAmplitude = dot(psi0, psi)
	    @info "Superposition of final state with the initializer" FinalAmplitude
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
    
	elseif Fast
    	return E, psi
    
    elseif StateAnalyzer
    
#		R = floor(Int64, L/2)	# Max distance
#		DensityFluctuations = zeros(R)
#		for r in 1:R
#			Tmp = 0
#			for j in 1:L
#				Tmp += DensityCorrelator[j,mod1(j+r,L)]
#			end
#			Tmp /= R
#			DensityFluctuations[r] = Tmp
#		end
		
		DensityCorrelator = GetDensityCorrelator(psi)
	    for i in 1:L
	    	# Locally defined
	    	SiteDensity[i] = expect(psi, "n"; sites=i)
			LocalE[i] = inner(psi', GetLocalHamiltonianMPO(sites, i, t, V, μ), psi)
    	end
    	Populations = GetLocalPopulation(psi)
    	
    	for b in 1:L-1
    		# Compute entropy for all possible bipartitions
    		Entropy[b] = GetVonNeumannEntropy(psi, b)
    	end
    	
    	# ShowDensityProfile(psi)    	
		return E, LocalE, psi, SiteDensity, DensityCorrelator, Entropy

	elseif Debug
		@warn "Mode Debug under construction."
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
    UserMode = "Fast"

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
    	
	elseif UserMode=="Fast"
		E, _ = Observables
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
