#!/usr/bin/julia

using ITensors, ITensorMPS

"""
Hellooo.
"""
function GetHamiltonian(Sites)::MPO
	os = OpSum()
	L = length(Sites)
	for j in 1:L
		os += (-1)^j,"N",j
	end
	return MPO(os, Sites)
end

function GetLocalPopulation(psi::MPS, d::Int64)
	"""
	Get the local population for each site and for each boson state, save it
	into a matrix.
	"""
	L = length(psi)
	Populations=zeros(L,d)
	KMatrix = zeros(d,d)
	for k in 1:d
		KMatrix[k,k] = 1
		Populations[:,k] = expect(psi,KMatrix)
		KMatrix[k,k] = 0 
	end
	
	return Populations
end

function GetLocalCoefficients(psiUniform::MPS, psi::MPS, d::Int64, Sites)
	"""
	Get the local population for each site and for each boson state, save it
	into a matrix.
	"""
	L = length(psi)
	Coefficients=zeros(L,d)
	for k in 1:d
		os = OpSum()
		os += 
		Coefficients[:,k] .= inner(psiUniform, MPO(os,Sites) ,psi)
	end
	
	return Coefficients
end

function RunTest()
	L = 3
	d = 1
	Sites = siteinds("Boson", L, dim=d)
	H = GetHamiltonian(Sites)

	UniformState = ["0" for k in 1:L]
	psiUniform = MPS(Sites, UniformState)

	# StartingState = [isodd(k) ? "1" : "2" for k in 1:L]
	# psi0 = MPS(Sites, StartingState)
	# E, psi = dmrg(H, psi0; nsweeps=2, maxdim=2, cutoff=1E-8)

	Populations = GetLocalPopulation(psiUniform, d)
	display(Populations)
end
