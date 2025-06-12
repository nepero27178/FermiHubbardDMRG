#!/usr/bin/julia

using ITensors, ITensorMPS
using Statistics
using Dates
using DelimitedFiles

# ------------------------------------------------------------------------------
# ------------------------------- Hamiltonian ----------------------------------
# ------------------------------------------------------------------------------

# Full hamiltonian

@doc raw"""
function GetHamiltonianMPO(
		sites::Any,
		t::Float64,
		V::Float64,
		μ::Float64;
		Φ=0
	)::MPO
		   
Returns: H as an MPO.

Build the 1D Fermi-Hubbard chain Hamiltonian as an MPO for given sites, with 
hopping `t`, interaction `V`, and chemical potential `μ`. The flux parameter `Φ`
controls the boundary conditions: if Φ=0, periodic boundary conditions are used.
Otherwise, twisted boundary conditions are implemented via a pseudo-gauge map
	```c_k -> exp(- iΦ * k/L) c_k```
which algebraically gives the hamiltonian as defined in code. The twisted
boundary conditions give a total Φ phase twist in the ground-state wavefunction
at the end of the chain.
"""
function GetHamiltonianMPO(
		sites::Any,
		t::Float64,
		V::Float64,
		μ::Float64;
		Φ=0
	)::MPO

    os = OpSum()
    L = length(sites)
    
    for j=1:L
    	# Separate: initialize Float64 or ComplexFloat64 
    	if Φ!==0
            os += -t * (cos(Φ/L) - im*sin(Φ/L)),"Cdag",j,"C",mod1(j+1,L)
	        os += -t * (cos(Φ/L) + im*sin(Φ/L)),"Cdag",mod1(j+1,L),"C",j
		elseif Φ==0
			os += -t,"Cdag",j,"C",mod1(j+1,L)
	        os += -t,"Cdag",mod1(j+1,L),"C",j
		end
		os += V,"N",j,"N",mod1(j+1,L)
		os += -μ,"N",j
    end

	return MPO(os, sites)
end

# Local hamiltonian

@doc raw"""
function GetLocalHamiltonianMPO(
		sites::Any,
		j::Int64,
		t::Float64,
		V::Float64,
		μ::Float64;
		Φ=0
	)::MPO
				   
Returns: Hj (the local hamiltonian) as an MPO.

Build the 1D Fermi-Hubbard chain local Hamiltonian at site `i` as an MPO, with 
hopping `t`, interaction `V`, and chemical potential `μ`. The flux parameter `Φ`
controls the boundary conditions: if Φ=0, periodic boundary conditions are used.
Otherwise, twisted boundary conditions are implemented via a pseudo-gauge map
	c_k -> exp(- iΦ * k/L) c_k
which algebraically gives the hamiltonian as defined in code. The twisted
boundary conditions give a total Φ phase twist in the ground-state wavefunction
at the end of the chain.
"""
function GetLocalHamiltonianMPO(
		sites::Any,
		j::Int64,
		t::Float64,
		V::Float64,
		μ::Float64;
		Φ=0
	)::MPO
    
    L = length(sites)
    if j<1 || j>L
    	error("Invalid site! Enter 1 ≤ j ≤ L.")
    end
    
    os = OpSum()
    
	# Separate: initialize Float64 or ComplexFloat64 
	if Φ!==0
		os += -t/2 * (cos(Φ/L) - im*sin(Φ/L)),"Cdag",mod1(j-1,L),"C",j
        os += -t/2 * (cos(Φ/L) + im*sin(Φ/L)),"Cdag",j,"C",mod1(j-1,L)
        os += -t/2 * (cos(Φ/L) - im*sin(Φ/L)),"Cdag",j,"C",mod1(j+1,L)
        os += -t/2 * (cos(Φ/L) + im*sin(Φ/L)),"Cdag",mod1(j+1,L),"C",j
	elseif Φ==0
		os += -t/2,"Cdag",mod1(j-1,L),"C",j
        os += -t/2,"Cdag",j,"C",mod1(j-1,L)
		os += -t/2,"Cdag",j,"C",mod1(j+1,L)
        os += -t/2,"Cdag",mod1(j+1,L),"C",j
	end
	
	os += V/2,"N",mod1(j-1,L),"N",j
	os += V/2,"N",j,"N",mod1(j+1,L)
    os += -μ,"N",j
    
    return MPO(os, sites)
end

# ------------------------------------------------------------------------------
# ------------------------------- Correlators ----------------------------------
# ------------------------------------------------------------------------------

# Total number of fermions and density-density correlator.

@doc raw"""
function GetTotalFermionNumber(
		psi::MPS
	)::Float64

Returns: expected particle number as a Float64.

This function calculates the expected total number of fermions on the MPS ground
state `psi`. 
"""
function GetTotalFermionNumber(
		psi::MPS
	)::Float64
	
	# Trivial function definition, just for notational coherence.
    return sum(expect(psi,"N"))
    
end

@doc raw"""
function GetDensityCorrelator(
		psi::MPS
	)::Matrix{Float64}

Returns: expected particle number as a Matrix{Float64}.

This function calculates the expected total number of fermions on the MPS ground
state `psi`. 
"""
function GetDensityCorrelator(
		psi::MPS
	)::Matrix{Float64}
	
	# Trivial function definition, just for notational coherence.
	return correlation_matrix(psi,"N","N")
									  
end

@doc raw"""
function GetBlockVariance(
		n::Vector{Float64},
		Cnn::Matrix{Float64};
		kVector=[]
	)::Vector{Float64}
	
Returns: `δn2M::Vector{Float64}`, the mean density variance over blocks. 
"""
function GetBlockVariance(
		n::Vector{Float64},		# Density vector
		Cnn::Matrix{Float64};	# Density-density correlation matrix
		kVector=[]				# Chosen block lengths
	)::Vector{Float64}
	
	# Connected density correlation matrix
	cCnn = Cnn .- n*n'
	
	# Ghost matrix to employ matricial PBC
	Ghost = [
		cCnn cCnn;
		cCnn zeros(size(cCnn));
	]

	L = size(Cnn,1)
	δn2M = Float64[]	# Unknown size
	
	if kVector!=[]
	
		if any(kVector .< 0 .|| kVector .> Int64(L/2))
			error("Invalid kVector! Enter 1 .≤ kVector .≤ L.")
		end
		
		δn2M = zeros(length(kVector))
		for (j,k) in enumerate(kVector)
			if k>0 && k<=Int64(L/2)
			
				Vectorδn2M = zeros(L)
				for s in 1:L
					Vectorδn2M[s] = sum(
						Ghost[
							s:s+k,
							s:s+k
						]
					)
				end
				
				push!(δn2M, mean(Vectorδn2M))
			
			end		
		end
		return δn2M
		
	elseif kVector==[]
		
		δn2M = zeros(Int64(L/2)+1)
		Matrixδn2M = zeros(Int64(L/2)+1,L)
		for (m,M) in enumerate(0:1:Int64(L/2)), s in 1:L
			Matrixδn2M[m,s] = sum(
				Ghost[
					s:s+M,
					s:s+M
				]
			)
		end
		
		δn2M .= mean(Matrixδn2M, dims=2)[:] # [:] Matrix{Float64} -> Vector{Float64}
		return δn2M
	
	end
end


# Superconducting order parameter and fluctuations correlator.

@doc raw"""
function GetSuperconductingPairingMPO(
		sites::Any,
		j::Int64;
		Φ=0
	)::MPO
	
Returns: superconducting order parameters as an MPO.

This function calculates for a given site `j` the superconducting pairing
operator on neighboring sites,
	C_j C_{j+1}
(note: usually the definition is hermitian-conjugate to this one).
"""
function GetSuperconductingPairingMPO(
		sites::Any,
		j::Int64;
		Φ=0
	)::MPO

	L = length(sites)
	if j<1 || j>L
    	error("Invalid site! Enter 1 ≤ j ≤ L.")
    end								   

	os = OpSum()
	if Φ==0
		os += "C",j,"C",mod1(j+1,L)
	elseif Φ!==0
		os += cos(Φ * (2*j+1)/L) - im*sin(Φ * (2*j+1)/L),"C",j,"C",mod1(j+1,L)
	end
	
	return MPO(os, sites)
end

@doc raw"""
function GetSuperconductingCorrelator(
		psi::MPS;
		Φ=0
	)::Matrix{Float64}
	
Returns: the superconducting pairing correlator as a Matrix{Float64}.

This function computes the pairing correlator for the superconducting order
parameter at all sites. 
"""
function GetSuperconductingCorrelator(
		psi::MPS;
		Φ=0
	)::Matrix{Float64}
	
	sites = siteinds(psi)
	L = length(sites)
	SuperconductingCorrelator = zeros(L,L)
	
	for j in 1:L, k in 1:L
		A = GetSuperconductingPairingMPO(sites,j)
		B = GetSuperconductingPairingMPO(sites,k)
		
		# Unusual definition, can get useful later in calculating mean correlators
		# SuperconductingCorrelator[j,mod1(j+k,L)] = inner(A, psi, B, psi)	
		
		SuperconductingCorrelator[j,k] = inner(A, psi, B, psi)
	end
	
	return SuperconductingCorrelator
									  
end

# ------------------------------------------------------------------------------
# -------------------------------- Operators -----------------------------------
# ------------------------------------------------------------------------------

# Density profile line printer

@doc raw"""
function ShowDensityProfile(
		psi::MPS
	)
	
Returns: nothing.

Prints on command line the density profile for the MPS `psi`, obtained
evaluating the density operator `N` on each site. The format output is a matrix
whose columns are the site indices (first) and the local denisty (second).
"""
function ShowDensityProfile(
		psi::MPS
	)
	
	DensityProfile = expect(psi, "N")
	println("Density profile:")
	for i in 1:length(DensityProfile)
		println("$i \t $(round(DensityProfile[i], digits=3))")
	end

	return
	
end

# Equal time Green's function

@doc raw"""
function GetEqualTimeGreensFunction(
		psi::MPS
	)::Matrix{Complex{Float64}}
	
Returns: G(i-j,0) as a ComplexF64.

This function calcualtes the equal time Green's function as a simple correlation
matrix for the operators `C` and `Cdag` on the MPS `psi`.
"""
function GetEqualTimeGreensFunction(
		psi::MPS
	)::Matrix{Complex{Float64}}
    
	# Trivial function definition, just for notational coherence.
	EqualTimeGreensFunction = - im .* correlation_matrix(psi,"C","Cdag") 
	return EqualTimeGreensFunction
end

# Von Neumann entropy

@doc raw"""
function GetVonNeumannEntropy(
		psi::MPS,
		b::Int64
	)::Float64
	
Returns: bipartite entropy at chain link `b` as a Float64.

Calculate Von Neumann entropy across bond `b` for the MPS `psi`; which is, the 
bipartite entropy for the bipartition (1,...,b) -- (b+1,...,L).
"""
function GetVonNeumannEntropy(
		psi::MPS,
		b::Int64
	)::Float64
	
	sites = siteinds(psi)
   
    # Perform SVD to all states except i.
    orthogonalize!(psi, b)

    # First argument: the tensor on which to perform the SVD, i.e. psi[i].
    # Second argument: the indices on which to perform the SVD, i.e. 
    # linkind(psi, i-1) and sites[i] (left link and vertical link).
    
    if b==1
	  _, S, _ = svd(psi[b], (sites[b],))
	else
	  _, S, _ = svd(psi[b], (linkind(psi, b-1), sites[b]))
	end

    # S is the diagonal matrix containing the singular values; SvN is the Von
    # Neumann entropy.
    
    SvN = 0.0
    for n in 1:dim(S, 1)
      p = S[n,n]^2
      SvN -= p * log(p)
    end
    return SvN
end

# Local population

@doc raw"""
function GetLocalPopulation(
		psi::MPS
	)::Matrix{Float64}
	
Returns: single-particle state populations as a Matrix{Float64}.

Gets the local population for each site and for each single-particle fermion 
state, and saves it into a matrix.
"""
function GetLocalPopulation(
		psi::MPS
	)::Matrix{Float64}
	
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

# Projectors
@doc raw"""
function GetHalfMIProjector(
		sites::Any
	)::MPO
	
Returns: Projector on the half-filled MI subspace.

This function computes the MPO projector onto the half-filled Mott-Insulating
Hilbert subspace.
"""
function GetHalfMIProjector(
		sites::Any
	)::MPO
	
	L = length(sites)
	
	states = [isodd(j) ? "0" : "1" for j in 1:L]	# 0101..
	e = MPS(sites, states)
	
	states = [isodd(j) ? "1" : "0" for j in 1:L]	# 1010..
	o = MPS(sites, states)
	
	E = projector(e) # Normalized projector
	O = projector(o) # Normalized projector
	HP = E+O
	return HP
	
end

@doc raw"""
function GetUnitaryMIProjector(
		sites::Any
	)::MPO
	
Returns: Projector on the unitary-filled MI subspace.

This function computes the MPO projector onto the unitary-filled Mott-Insulating
Hilbert subspace.
"""
function GetUnitaryMIProjector(
		sites::Any
	)::MPO
	
	L = length(sites)
	
	states = ["1" for j in 1:L]						# 1111..
	u = MPS(sites, states)
	
	UP = projector(u)
	return UP
	
end
