#!/usr/bin/julia

using Base.Threads
using DelimitedFiles  # For writedlm

# ------------------------------ Boundaries sweep ------------------------------

@doc raw"""
function BoundariesSweep(
		L::Int64,
		VV::Array{Float64},
		DMRGParameters::Vector{Any},                         
		FilePathOut::String; 
		μ0=0.0
	)
	
Returns: none (inline print).

# TODO To be completely understood..

This function computes the estimated positions for the gapped/gapless phase
transitions boundaries. Keeping the fermionic parity constant, we add/remove 2
femions from the chain and compute the respective ground states keeping each
time fixed the number of fermions. The difference in energy, apart from simple
chemical potential shifts, must be all due to the hopping-interaction part only
parametrized by V/t.
"""
function BoundariesSweep(
		L::Int64,
		VV::Array{Float64},
		DMRGParameters::Vector{Any},                         
		FilePathOut::String; 
		μ0=0.0
	)
    
    # ModelParameters: we set V=1.0, η=0.0
    DataFile = open(FilePathOut, "a")
    println("Extracting boundaries...")
    
    l = length(VV)
    for (j, V) in enumerate(VV)
        println("Running DMRG for V=$(round(V, digits=3)), μ=$μ0", 
        	" (simulation $j/$l for L=$L)")
        
        # Use @sync to wait for all tasks to complete
        @sync begin
            # Create tasks for each DMRG call
            task1 = @spawn RunDMRGAlgorithm(
            	[L, L/2-1, 1.0, V, μ0, 0.0],
				DMRGParameters,
				"Fast", # UserMode
				true;	# FixedN
				verbose=false
			)
            task2 = @spawn RunDMRGAlgorithm(
            	[L, L/2, 1.0, V, μ0, 0.0], 
				DMRGParameters,
				"Fast", # UserMode
				true;	# FixedN
				verbose=false
			)
            task3 = @spawn RunDMRGAlgorithm(
				[L, L/2+1, 1.0, V, μ0, 0.0],
				DMRGParameters,
				"Fast", # UserMode
				true;	# FixedN
				verbose=false
			)

            # Wait for all tasks to complete and collect results
            EDown = fetch(task1)
            E = fetch(task2)
            EUp = fetch(task3)

            # Calculate chemical potentials
            ΔUp = EUp - E
            ΔDown = E - EDown
            μUp = ΔUp/2 + μ0
            μDown = - ΔDown - 2*μ0 # Sign problem otherwise

            # Write results to the file
            write(DataFile, "$L; $V; $E; $μUp; $μDown\n")
        end
    end
    
    close(DataFile)
    println("Data for L=$L saved on file!")
end

# ------------------------------ Horizontal sweep ------------------------------

@doc raw"""
function HorizontalSweep(
		L::Int64,
		VV::Array{Float64},
		μ0::Float64,
		DMRGParameters::Vector{Any},                         
		FilePathOut::String
	)
	
Returns: none (inline print).

Run horizontal sweep (fixed μ0) to extract observables and correlator at
increasing V/t.
"""
function HorizontalSweep(
		L::Int64,
		VV::Array{Float64},
		μ0::Float64,
		DMRGParameters::Vector{Any},                         
		FilePathOut::String
	)
    
    DataFile = open(FilePathOut, "a")
    println("Performing horizontal sweep at μ=$μ0...")
    
    l = length(VV)
    for (j, V) in enumerate(VV)
    	println("Running DMRG for V=$(round(t, digits=3)), μ=$μ0", 
        	" (simulation $j/$l for L=$L)")
        
        # Use @sync to wait for all tasks to complete
        @sync begin
            E, Γ, eΓ, O, eO =  RunDMRGAlgorithm(
            	[L, L/2, 1.0, V, μ0, 0.0], 
				DMRGParameters,
				"Correlators",
				false;	# FixedN
				verbose=false
			)

            # Write results to the file
            write(DataFile, "$L; $V; $E; $Γ; $eΓ; $O, $eO\n")
        end
    end
    
    close(DataFile)
    println("Data for L=$L saved on file!")
end


# ----------------------------- Rectangular sweep ------------------------------

@doc raw"""

"""
function RectangularSweepBoundaries(
		L::Int64,
		N::Int64,
		VV::Vector{Float64},
		μμ::Vector{Float64},
		DMRGParametersXY::Vector{Any},
		DMRGParametersIF::Vector{Any},
		FilePathIn::String,
		FilePathOut::String
	)
    
    @warn "Function under construction."
#    
#    DataFile = open(FilePathOut,"a")

#	# Take data from fitting data of horizontal sweeps to separate MI from SF 
#	# ( use ΔE^+(∞) and ΔE^-(∞) )
#	BoundariesData = readdlm(FilePathIn, ',', Float64, '\n'; comments=true)
#    VVFitted = BoundariesData[:,1]

#    for (j,V) in enumerate(VV)
#		# We take the best approximating V ∈ VV, to assess whether we are in 
#		# the XY or IF phase.

#		Index = argmin(abs.(VVFitted .- V)) # this always works, gives the best approximation
#		μUp = BoundariesData[Index,2][1]
#		μDown = -BoundariesData[Index,3][1]
#    	
#		println("\nV=$V, phase boundaries: μ^+=$(μUp), μ^-=$(μDown)")

#		CachedRho = 1
#		
#        for (m,μ) in enumerate(μμ)
#        
#            ModelParameters = [L, N, 1.0, V, μ, 0.0]
#			inIsingRegion = false
#			
#			if (μ>=(μDown) && μ<=(μUp))
#				inIsingRegion = true
#			end

#			println("Running DMRG for L=$L, V=$(round(V, digits=3)), ",
#				"μ=$(round(μ, digits=3)) (simulation $m/$(length(μμ)) in μ, ",
#				"$j/$(length(VV)) in V, IF=$inIsingRegion)")

#			# Go on from here...

#			if inIsingRegion
#		        Results = RunDMRGAlgorithm(
#                    ModelParameters,
#                    DMRGParametersIF,
#                    "OrderParameters";
#                    FixedN=false,
#                    RandomPsi0=false
#                )
#				E, NTotAvg, nVariance, aAvg = Results
#				
#				if m>1
#					Rho = NTotAvg/L
#					K = (Rho - CachedRho) / (μ - μμ[m-1])	# Compressibility
#					K /= (Rho^2)
#					CachedRho = Rho							# Next step
#		        elseif m==1
#		        	K = NaN									# Avoid segmentation fault
#		        end
#		        
#		        write(DataFile,"$J; $μ; $E; $nVariance; $aAvg; $K # MI\n")
#		    else
#		    	Results = RunDMRGAlgorithm(ModelParameters,
#				                           DMRGParametersSF,
#				                           "OrderParameters";
#	  		                               FixedN=false,
#										   RandomPsi0=true)
#				E, NTotAvg, nVariance, aAvg = Results
#				
#				if m>1
#					Rho = NTotAvg/L
#					K = (Rho - CachedRho) / (μ - μμ[m-1])	# Compressibility
#					K /= (Rho^2)
#					CachedRho = Rho							# Next step
#		        elseif m==1
#		        	K = NaN									# Avoid segmentation fault
#		        end
#				
#		        write(DataFile,"$J; $μ; $E; $nVariance; $aAvg; $K # SF\n")
#		    end
#        end
#    end

#    close(DataFile)
end

@doc raw"""

"""
function RectangularSweep(
		L::Int64,
		N::Int64,
		VV::Vector{Float64},
		μμ::Vector{Float64},
		DMRGParameters::Vector{Any},
		FilePathOut::String
	)
    
    global ε
    DataFile = open(FilePathOut,"a")

    for (j,V) in enumerate(VV), (m,μ) in enumerate(μμ)
        
        ModelParameters = [L, N, 1.0, V, μ, 0.0]

		println("Running DMRG for L=$L, V=$(round(V, digits=3)),", 
			"μ=$(round(μ,digits=3)) (simulation $m/$(length(μμ)) in μ, ",
			"$j/$(length(VV)) in V)")

		E = 0 # Initalize energy to save variable
		k = 0 # Initalize compressibility to save variable
		D = 0 # Initalize stiffness to save variable
        @sync begin
            # Create tasks for each DMRG call
            task1 = @spawn RunDMRGAlgorithm(
            	[L, L/2, 1.0, V, μ, 0.0],		# Central
				DMRGParameters,
				"Fast", # UserMode
				true;	# FixedN
				verbose=false
			)
            task2 = @spawn RunDMRGAlgorithm(
            	[L, L/2+2, 1.0, V, μ, 0.0], 	# Add two particles
				DMRGParameters,
				"Fast", # UserMode
				true;	# FixedN
				verbose=false
			)
            task3 = @spawn RunDMRGAlgorithm(
				[L, L/2-2, 1.0, V, μ, 0.0],		# Remove two particles
				DMRGParameters,
				"Fast", # UserMode
				true;	# FixedN
				verbose=false
			)
			task4 = @spawn RunDMRGAlgorithm(
            	[L, L/2, 1.0, V, μ, ε], 		# Rotate clockwise
				DMRGParameters,
				"Fast", # UserMode
				false;	# FixedN
				verbose=false
			)
			task5 = @spawn RunDMRGAlgorithm(
            	[L, L/2, 1.0, V, μ, -ε], 		# Rotate counter-clockwise
				DMRGParameters,
				"Fast", # UserMode
				false;	# FixedN
				verbose=false
			)

            # Wait for all tasks to complete and collect results
            E = fetch(task1)
            EAdd = fetch(task2)
            ERemove = fetch(task3)
            EClock = fetch(task4)
            ECounterClock = fetch(task5)
            
            k = 4/( L * (EAdd+ERemove-2*E) )
            D = pi*L * (EClock+ECounterClock-2*E)/(4ε^2)
        end
        
        write(DataFile,"$V; $μ; $E; $k; $D\n")
	end

    close(DataFile)
end

# -------------------------------- Path sweep ----------------------------------

# UNDER CONSTRUCTION!

# function PathSweep(Path::NamedTuple(Name::String, Values::Vector{Tuple{Int64, Int64}}),
# 				   L::Int64,
# 				   N::Int64,
# 				   nmax::Int64,
# 				   DMRGParametersMI::Vector{Any},
# 				   DMRGParametersSF::Vector{Any},
# 				   FilePathIn::String,				# To evaluate if a given point is MI or SF
# 				   FilePathOut::String)
				   
# 	Sizes, Couplings, Energies, UpBoundaries, DownBoundaries, _, _ = readdlm(FilePathIn, ';', Float64, '\n'; comments=true)
# 	IndicesList = findall(==(L), Sizes)
				   
# 	for Point in Path.Values
# 		J, μ = Point
				                      
# 		# Possible improvement: we know how many simulations have been performed for each L
#     	Index = IndicesList[findall(==(J), Couplings[IndicesList])]
#     	μUp = UpBoundaries[Index]
#     	μDown = DownBoundaries[Index]
        
#         ModelParameters = [L, N, nmax, J, μ]
# 		inMottLobe = false
		
# 		if (μ>=μDown && μ<=μUp)
# 			inMottLobe=true
# 		end
		
# 		if inMottLobe
#             _, _, _, C = RunDMRGAlgorithm(ModelParameters,
#                                           DMRGParametersMI;
#                                           ComputeC = true,
#                                           FixedN = false)
#         else
#         	_, _, _, C = RunDMRGAlgorithm(ModelParameters,
#                                           DMRGParametersSF;
# 	                                      ComputeC = true,
# 	                                      FixedN = false)
#         end                
        
#         # Perform DFT
#         D = 0
#         q = 2*pi/L
#         for r in 1:length(C)
#         	D += (-1)^(2*r/L) * C[r] # Numerically smarter than using the imaginary unit
#         end
#         K = 1/(L*D)
        
#         DataFile = open(FilePathOut, "a")
#         write(DataFile, "$L, $J, $μ, $D, $K")
		
# 	end
# 	close(DataFile)
#     println("Data for L=$L saved on file!")
# end

