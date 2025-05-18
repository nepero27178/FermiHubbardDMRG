#!/usr/bin/julia

using Base.Threads
using DelimitedFiles  # For writedlm

# TODO Restart here after having finished checking for the XXZ mapping.

# ------------------------------ Boundaries sweep ------------------------------

@doc raw"""
function BoundariesSweep(
		L::Int64,
		nmax::Int64,
		JJ::Array{Float64},
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
parametrized by t/V.
"""
function BoundariesSweep(
		L::Int64,
		tt::Array{Float64},
		DMRGParameters::Vector{Any},                         
		FilePathOut::String; 
		μ0=0.0
	)
    
    @warn "This function is under construction!"
    
    # ModelParameters: we set V=1.0, η=0.0
    DataFile = open(FilePathOut, "a")
    println("Extracting boundaries...")
    
    l = length(tt)
    for (j, t) in enumerate(tt)
        println("Running DMRG for t=$(round(t, digits=3)), μ=$μ0". 
        	" (simulation $j/$l for L=$L)")
        
        # Use @sync to wait for all tasks to complete
        @sync begin
            # Create tasks for each DMRG call
            task1 = @spawn RunDMRGAlgorithm(
            	[L, L/2-2, t, 1.0, μ0, 0.0],
				DMRGParameters,
				"Fast", # UserMode
				true;	# FixedN
				verbose=false
			)
            task2 = @spawn RunDMRGAlgorithm(
            	[L, L/2, t, 1.0, μ0, 0.0], 
				DMRGParameters,
				"Fast", # UserMode
				true;	# FixedN
				verbose=false
			)
            task3 = @spawn RunDMRGAlgorithm(
				[L, L/2+2, t, 1.0, μ0, 0.0],
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
            μDown = - ΔDown/2 - μ0 # Sign problem otherwise

            # Write results to the file
            write(DataFile, "$L; $t; $E; $μUp; $μDown\n")
        end
    end
    
    close(DataFile)
    println("Data for L=$L saved on file!")
end

# ------------------------------ Horizontal sweep ------------------------------

@doc raw"""
function HorizontalSweep(
		L::Int64,
		tt::Array{Float64},
		μ0::Float64,
		DMRGParameters::Vector{Any},                         
		FilePathOut::String
	)
	
Returns: none (inline print).

Run horizontal sweep (fixed μ0) to extract observables and correlator at
increasing t.
"""
function HorizontalSweep(
		L::Int64,
		tt::Array{Float64},
		μ0::Float64,
		DMRGParameters::Vector{Any},                         
		FilePathOut::String
	)
    
    DataFile = open(FilePathOut, "a")
    println("Performing horizontal sweep at μ=$μ0...")
    
    l = length(tt)
    for (j, t) in enumerate(tt)
    	println("Running DMRG for t=$(round(t, digits=3)), μ=$μ0". 
        	" (simulation $j/$l for L=$L)")
        
        # Use @sync to wait for all tasks to complete
        @sync begin
            E, Γ, eΓ, O, eO =  RunDMRGAlgorithm(
            	[L, L/2, t, V, μ0, η], 
				DMRGParameters,
				"Correlators",
				false;	# FixedN
				verbose=false
			)

            # Write results to the file
            write(DataFile, "$L; $t; $E; $Γ; $eΓ; $O, $eO\n")
        end
    end
    
    close(DataFile)
    println("Data for L=$L saved on file!")
end


# ----------------------------- Rectangular sweep ------------------------------

# Go on from here...
# First task by now: run some BoundariesSweeps to see if we delimit a sensed phase domain

function RectangularSweepBoundaries(
		L::Int64,
		N::Int64,
		tt::Vector{Float64},
		μμ::Vector{Float64},
		DMRGParametersMI::Vector{Any},
		DMRGParametersSF::Vector{Any},
		FilePathIn::String,
		FilePathOut::String
	)
    """
    Calculate the variance of the number of particles on site i and ⟨b_i⟩ for a
    range of hopping J and chemical potential μ values. Results are saved to a file.
    """
    
    DataFile = open(FilePathOut,"a")

	# Take data from fitting data of horizontal sweeps to separate MI from SF 
	# ( use ΔE^+(∞) and ΔE^-(∞) )
	BoundariesData = readdlm(FilePathIn, ',', Float64, '\n'; comments=true)
    ttFitted = BoundariesData[:,1]

    for (j,t) in enumerate(tt)
		# We take the best approximating t ∈ tt, to assess whether we are in the MI or SF phase

		# Index = findall(==(t), BoundariesData[:,1]) # this would work if there is the EXACT t in the fit results
		Index = argmin(abs.(ttFitted .- t)) # this always works, gives the best approximation
		μUp = BoundariesData[Index,2][1]
		μDown = -BoundariesData[Index,3][1]
    	
		println("\nt=$t, phase boundaries: μ^+=$μUp, μ^-=$μDown")

		CachedRho = 1
		
        for (m,μ) in enumerate(μμ)
        
            ModelParameters = [L, N, t, V, μ, η]
			inMottLobe = false
			
			if (μ>=μDown && μ<=μUp)
				inMottLobe = true
			end

			println("Running DMRG for L=$L, t=$(round(t, digits=3)), ",
			"μ=$(round(μ, digits=3)) (simulation $m/$(length(μμ)) in μ, ",
			"$j/$(length(tt)) in t, MI=$inMottLobe)")

			if inMottLobe
		        Results = RunDMRGAlgorithm(
                    ModelParameters,
                    DMRGParametersMI,
                    "OrderParameters";
                    FixedN=false,
                    RandomPsi0=false
                )
				E, NTotAvg, nVariance, aAvg = Results
				
				if m>1
					Rho = NTotAvg/L
					K = (Rho - CachedRho) / (μ - μμ[m-1])	# Compressibility
					K /= (Rho^2)
					CachedRho = Rho							# Next step
		        elseif m==1
		        	K = NaN									# Avoid segmentation fault
		        end
		        
		        write(DataFile,"$J; $μ; $E; $nVariance; $aAvg; $K # MI\n")
		    else
		    	Results = RunDMRGAlgorithm(ModelParameters,
				                           DMRGParametersSF,
				                           "OrderParameters";
	  		                               FixedN=false,
										   RandomPsi0=true)
				E, NTotAvg, nVariance, aAvg = Results
				
				if m>1
					Rho = NTotAvg/L
					K = (Rho - CachedRho) / (μ - μμ[m-1])	# Compressibility
					K /= (Rho^2)
					CachedRho = Rho							# Next step
		        elseif m==1
		        	K = NaN									# Avoid segmentation fault
		        end
				
		        write(DataFile,"$J; $μ; $E; $nVariance; $aAvg; $K # SF\n")
		    end
        end
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

