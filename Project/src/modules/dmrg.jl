#!/usr/bin/julia

using ITensors, ITensorMPS
using Statistics
using Dates
using DelimitedFiles

# ------------------------------------------------------------------------------
# --------------------------------- Physics ------------------------------------
# ------------------------------------------------------------------------------

# Full hamiltonian

H = GetHamiltonian(sites,
				   t::Float64,
				   V::Float64,
				   μ::Float64;
				   Φ=0)
    """
    Construct the 1D Fermi-Hubbard chain Hamiltonian as an MPO for given sites,
    hopping `t`, interaction `V`, and chemical potential `μ`. The parameter Φ
    controls the boundary conditions: if Φ=0, periodic boundary conditions are
    used. Otherwise twisted boundary conditions are implemented via a 
    pseudo-gauge map
    	c_k -> exp(- iΦ * k/L) c_k
    which algebraically gives the hamiltonian below.
    """
    os = OpSum()
    L = length(sites)
    
    for j=1:L
    	# Separate: initialize Float64 or ComplexFloat64 
    	if Φ!==0
            os += -t * (cos(Φ/L) + im*sin(Φ/L)),"Cdag",j,"C",mod1(j+1,L)
	        os += -t * (cos(Φ/L) - im*sin(Φ/L)),"Cdag",mod1(j+1,L),"C",j
		elseif Φ==0
			os += -t,"Cdag",j,"C",mod1(j+1,L)
	        os += -t,"Cdag",mod1(j+1,L),"C",j
		end
		os += V,"N",j,"N",mod1(j+1,L)
		os += -μ,"N",j
    end
    
    return MPO(os,sites)
end

# Local hamiltonian

function GetLocalH(sites,
				   j::Int64,
				   t::Float64,
				   V::Float64,
				   μ::Float64;
				   Φ=0)
    """
    Construct the local Hamiltonian term MPO for site `i` in the 1D 
    Fermi-Hubbard model with hopping `t`, NN interaction `V`, and chemical
    potential `μ`.
    """
    
    if j<1 || j>L
    	error("Invalid site! Enter 1≤j≤L.")
    	break
    end
    
    os = OpSum()
    L = length(sites)
    
	# Separate: initialize Float64 or ComplexFloat64 
	if Φ!==0
		os += -t/2 * (cos(Φ/L) + im*sin(Φ/L)),"Cdag",mod1(j-1,L),"C",j
        os += -t/2 * (cos(Φ/L) - im*sin(Φ/L)),"Cdag",j,"C",mod1(j-1,L)
        os += -t/2 * (cos(Φ/L) + im*sin(Φ/L)),"Cdag",j,"C",mod1(j+1,L)
        os += -t/2 * (cos(Φ/L) - im*sin(Φ/L)),"Cdag",mod1(j+1,L),"C",j
	elseif Φ==0
		os += -t/2,"Cdag",mod1(j-1,L),"C",j
        os += -t/2,"Cdag",j,"C",mod1(j-1,L)
		os += -t/2,"Cdag",j,"C",mod1(j+1,L)
        os += -t/2,"Cdag",mod1(j+1,L),"C",j
	end
	
	os += V/2,"N",mod1(j-1,L),"N",j
	os += V/2,"N",j,"N",mod1(j+1,L)
    os += -μ,"N",j
end

# Total number operator

function GetNumber(sites)
    """
    Calculate the total number operator `N` as a MPO for the given sites.
    """
    os = OpSum()
    for j=1:length(sites)
        os += "N",j
    end

    return MPO(os, sites)
end

function GetSuperConductingPairing(sites,
								   j::Int64,
								   Φ=0)

	if j<1 || j>L
    	error("Invalid site! Enter 1≤j≤L.")
    	break
    end								   

	os = OpSum()
	L = length(sites)
	if Φ==0
		os += "C",j,"C",mod1(j+1,L)
	elseif Φ!==0
		os += cos(Φ * (2*j+1)/L) + im*sin(Φ * (2*j+1)/L),"C",j,"C",mod1(j+1,L)
	end
	
	return MPO(os, sites)
end

function GetSuperConductingCorrelator(sites,
									  psi::MPS,
									  j::Int64,
									  k::Int64,
									  Φ=0)
	
	if j<1 || j>L
    	error("Invalid site! Enter 1≤j≤L.")
    	break
    end
    
    if k<=1 || k>=L
    	error("Invalid distance! Enter 1<k<L.")
    	break
    end
	
	# TODO Extend to non-zero flux threading the ring.
	
	A = GetSuperConductingCorrelator(sites,j)
	B = GetSuperConductingCorrelator(sites,mod1(j+k,L))
	inner(A, psi, B, psi)
	
	return 
									  
end

function GetNumberCorrelator(psi::MPS)
	
	"""
	Trivial function definition, just for notational coherence.
	"""
	
	return correlation_matrix(psi,"N","N")
									  
end

# Two points correlator # TODO

function GetEqualTimeGreensFunction(psi::MPS)
    
    """
	Trivial function definition, just for notational coherence.
	"""
	
	EqualTimeGreensFunction = - im .* correlation_matrix(psi,"C","Cdag") 
	return EqualTimeGreensFunction
end

# Von Neumann entropy

function GetVonNeumannEntropy(psi::MPS, sites, b::Int64)
    """
    Calculate Von Neumann entropy across bond b; which is, the bipartite entropy
    for the bipartition (1,...,b) -- (b+1,...,L).
    """
    # Perform SVD to all states except i, to prepare for calculating a local observable (end of Lect 4 page 4)
    orthogonalize!(psi, b)

    # First argument: the tensor on which to perform the SVD, i.e. psi[i]
    # Second argument: the indices on which to perform the SVD, i.e. linkind(psi, i-1) and sites[i] (left link and vertical link)
    
    if b==1
	  _, S, _ = svd(psi[b], (sites[b],))
	else
	  _, S, _ = svd(psi[b], (linkind(psi, b-1), sites[b]))
	end

    # S is the diagonal matrix containing the singular values
    SvN = 0.0 				# von Neumann entropy
    for n in 1:dim(S, 1)
      p = S[n,n]^2
      SvN -= p * log(p)
    end
    return SvN
end

# Local population # TODO

function GetLocalPopulation(psi::MPS)
	"""
	Get the local population for each site and for each boson state, save it
	into a matrix.
	"""
	L = length(psi)
	Populations=zeros(L,2)
	KMatrix = zeros(2,2)
	for k in 1:2
		KMatrix[k,k] = 1
		Populations[:,k] = expect(psi,KMatrix)
		KMatrix[k,k] = 0 
	end
	
	return Populations
end

# TODO GO ON FROM HERE

# ------------------------------------------------------------------------------
# ------------------------------------ DMRG ------------------------------------ 
# ------------------------------------------------------------------------------

function SetStartingState(sites)
    """
    Create an initial MPS with `N` particles on the fermionic chain; this move
    initializes the system at half-filling.
    """
    L = length(sites)
    states = ["1" for _ in 1:L]

    return MPS(sites, states)
end

function RunDMRGAlgorithm(ModelParameters::Vector{Float64},
						  DMRGParameters::Vector{Any},
						  UserMode::String;
						  verbose=false,				# Do not print at line
						  FixedN=false,					# Let N vary
						  pbc=true)	                    # Periodic boundary conditions
    
    """
    Run DMRG algorithm with chosen parameters, and return results.
    The on-site interaction is absent in the fermionic polarized case.
    Input:
        - ModelParameters: array of [L::Int64,
        							 N::Int64,
        							 t::Float64,
        							 V::Float64,
        							 μ::Float64]
        - DMRGParameters: array of [nsweeps::Int64, maxdim::Int64, 
        							cutoff::Vector{Float64}]
    Parametric input:
        - ComputeAllObservables: boolean variable, if false only E, 
          nVariance and Γ are computed (default: false)
        - ComputeGamma: boolean variable, if false all positional correalators
          are not extracted (default: false)
        - ComputeC: boolean variable, if false all positional correalators are
          not extracted (default: false)
    Output:
        - Results of DMRG optimization and relevant observables
    """
    
    L, N = Int64.(ModelParameters[1:2])
    t, Φ, V, μ = ModelParameters[3:6]
    nsweeps, maxdim, cutoff = DMRGParameters[1:3]

	# Evaluate UserMode
	ModeErrorMsg = "Input error: use as argument \"OrderParameters\"
- \"Correlator\", \"Debug\" or \"Fast\":
- \"OrderParameters\": return E, nVariance, aAvg;
- \"Correlator\": return E, Γ, eΓ;
- \"Debug\": return E, LocalE, nMean, nVariance;
- \"Fast\": return E."
	
	OrderParameters=false
	Correlator=false
	Debug=false
	Fast=false
	
	if UserMode=="OrderParameters"
		OrderParameters=true
		nVariance = 0		 		# Variance on central site
   		aAvg = 0		 			# <a> on central site
	
	elseif UserMode=="Correlator"
		Correlator=true
		# Initialize later
		
	elseif UserMode=="Debug"
		Debug=true
		nMean = zeros(L)				# Mean number of particles per site
		nVariance = zeros(L) 			# Variance on n_i, for all sites i
    	LocalE = zeros(L)				# Local contribution to the energy
    	Populations = zeros(L,nmax+1) 	# Single site populations
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
        println("Model parameters: L=$L, N=$N, nmax=$nmax, J=$J, U=$U, μ=$μ.")
        println("Starting simulation...\n")
    end

    # Calculate hamiltonian and number operators
    sites = siteinds("Fermion", L)
    H = GetHamiltonian(sites, t, Φ, V, μ; pbc)
    NTot = GetNumber(sites)
    
    # Set starting state and print observables
    if RandomPsi0
        psi0 = randomMPS(sites)					# Initialize random state
    else
        psi0 = SetStartingState(sites, N, d)	# Initialize Mott state
    end

    if verbose
        N0 = inner(psi0', NTot, psi0)
        E0 = inner(psi0', H, psi0)
        println("Expectation values on the initial state: N=$N0, E=$E0\n")
    end

    # Run DMRG algorithm and print results
    E, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, outputlevel=verbose)

    # Sanity checks: calculate whether Ntot has been conserved, and 
    # the found ground state is actually an eigenstate of H.
    if verbose
        VarE = inner(H, psi, H, psi) - E^2
		NTotAvg = inner(psi', NTot, psi)
        println("\nFinal N=$(round(NTotAvg, digits=3)), ", 
        		"E=$(round(E,digits=4)), ",
        		"VarE=$(round(VarE,digits=4))\n")
    end

    if OrderParameters
    
    	Index = ceil(Int64, L/2)
    	
		if !verbose
			NTotAvg = inner(psi', NTot, psi)
		end    	
		
		nVariance = GetNumberVariance(psi, sites, Index)
	    aAvg = expect(psi, "a"; sites=Index)
		
		return E, NTotAvg, nVariance, aAvg
    
    elseif Correlator

		"""
		Procedure: discard half of the lattice (the first quarter and the last);
		what remains is used as a domain to sweep across, collecting 
		symmetrically two-points and density-density correlators; these last
		are then avereaged; error is taken as the statistical standard 
		deviation.
		"""
		
		Trash = floor(Int64, L/4)	# How many to discard from each side?
		Start = Trash+1				# Start of useful segment
		Stop = L-Trash				# End of useful segment
		
		Segment = L-2*Trash			# Segment length
		Spacings = collect(Int64, 2:2:Segment)

		Γ = zeros(Float64, length(Spacings))
		eΓ = zeros(Float64, length(Spacings))
		
		# (Center) symmetric sweep for Γ
		for r in Spacings

		    ΓTmpArray = [] # stores the different Γ(r)

		    i = Start
		    j = Start+r-1

		    while j<=Stop
		        ΓOp = GetTwoPointCorrelator(sites, i, j)
		        ΓTmp = inner(psi', ΓOp, psi)
		        push!(ΓTmpArray, ΓTmp)	

		        i += 1
		        j += 1   
		    end

		    Γ[Int64(r/2)] = mean(ΓTmpArray)
		    eΓ[Int64(r/2)] = std(ΓTmpArray)
        end
        
        eΓ[end] = eΓ[end-1] # Correct NaN
		
        return E, Γ, eΓ
    
    elseif Debug
		
	    for i in 1:L
	    	# Locally defined
	    	nMean[i] = expect(psi, "n"; sites=i)
			nVariance[i] = GetNumberVariance(psi, sites, i)
			LocalE[i] = inner(psi', GetLocalH(sites, i, J, U, μ), psi)
    	end
    	Populations = GetLocalPopulation(psi, d)
    	
    	for b in 1:L-1
    		# Compute entropy for all possible bipartitions
    		Entropy[b] = GetVonNeumannEntropy(psi, sites, b)
    	end
    	
    	return E, nMean, nVariance, LocalE, Populations, Entropy
    	
    elseif Fast
    	return E
    end
end

# ------------------------------------------------------------------------------
# ------------------------------------ Main ------------------------------------ 
# ------------------------------------------------------------------------------

function main()

	"""
    If the script is called directly from terminal, run one full DMRG routineù
    with the following parameters.
    """
    
    L = 30
    N = 30
    nmax = 4
    J = 0.06		# MI: 0.06 / SF: 0.33
    μ = 0.4			# MI: 0.4 / SF: 0.8

    nsweep = 10
    maxlinkdim = [10,50,75,200,500]
    cutoff = [1E-8]
    DMRGParameters = [nsweep, maxlinkdim, cutoff]
    
    UserMode = "Debug" # "OrderParameters" / "Correlator" / "Debug" / "Fast"

	if UserMode=="Fast"
		for i in 1:100
    		print("\r")
    		for j in 1:i
    			print(" ")
    		end
    		printstyled(">>> Speedin' >>>", color=:cyan)
    		sleep(0.05)
    	end
    	print("\n")
	end
	
    # ------------------------ Results for a single run ------------------------
    
    Observables = RunDMRGAlgorithm([L, N, nmax, J, μ],
    							    DMRGParameters,
    							    UserMode; 
    								verbose=true,
    								FixedN=false,
    								RandomPsi0=false)
    								
	if UserMode=="OrderParameters"
		E, nVariance, aAvg = Observables
		println("Results of the simulation:
Energy of ground state: $(round.(E, digits=4))
Number variance on central site: $(round.(nVariance, digits=4))
Average of <a> on central site: $(round.(aAvg, digits=4))")
	
	elseif UserMode=="Correlator"
		E, Γ, eΓ = Observables
		println("Results of the simulation:
Energy of ground state: $(round.(E, digits=4))
Green's function: $(round.(Γ, digits=4))
Error on Green's function: $(round.(eΓ, digits=4))")
		
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

		for i in 1:100
    		print("\r")
    		for j in 1:i
    			print(" ")
    		end
    		printstyled(">>> I'm speed >>>", color=:cyan)
    		sleep(0.05)
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
