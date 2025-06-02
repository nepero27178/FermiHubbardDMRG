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

@doc raw"""
function main()

Returns: none.

Note: the main() functions at this level are the only functions in the whole
project referring to externally defined global variables. All simulations
variables have been defined in `./setup/simulations-setup.jl` and are here
referred to with their names. Check the setup file to modify the definitions.
"""
function main()
    
	ModeErrorMsg = "Input error: use option --horizontal or --rectangular"
	
	if length(ARGS) != 1
		# If user does not specify the user mode
		error(ModeErrorMsg)
		exit()
	else
		
		UserMode = ARGS[1]
        
        # -------------------------- Horizontal sweep --------------------------
        
        if UserMode=="--horizontal"
		
			# Horizontal sweeps
			LL = HorizontalLL				# Imported from setup
	    	VV = HorizontalVV				# Imported from setup
	    	μ0μ0 = Horizontalμμ				# Imported from setup

			Waiting=true
			print("Choose your method: (Boundaries/StateAnalysis/Full) ")
			UserSubMode = readline()
			
			BoundariesDirPathOut = ""
			StateAnalysisDirPathOut = ""
			
			while Waiting
			
				if UserSubMode=="Boundaries" || UserSubMode=="StateAnalysis" || UserSubMode=="Full"
					
					Waiting=false
					BoundariesDirPathOut = PROJECT_ROOT * "/simulations/horizontal-sweep/boundaries"
					mkpath(BoundariesDirPathOut) 
					
					StateAnalysisDirPathOut = PROJECT_ROOT * "/simulations/horizontal-sweep/state-analysis"
					mkpath(StateAnalysisDirPathOut) 
					
				else
					print("Invalid input. Choose your method: (Boundaries/StateAnalysis/Full) ")
					UserSubMode = readline()
				end
			
			end
	    	
	    	for μ0 in μ0μ0
	    	
	    		# Prepare output files paths
	    		BoundariesFilePathOut = ""
	    		StateAnalysisFilePathOut = ""
	    		
	    		if UserSubMode=="Boundaries" || UserSubMode=="Full"
		    		BoundariesFilePathOut = BoundariesDirPathOut * "/μ0=$(μ0)_L=$(LL).txt"
		    		
					# Prepare file headers
					BoundariesFile = open(BoundariesFilePathOut,"w")
						write(BoundariesFile,"# SFH Charge gaps. μ0=$μ0, nsweeps=$nsweeps, cutoff=$cutoff\n")
						write(BoundariesFile,"# L; V; hΔ-1; hΔ+1; hΔ-2; hΔ+2; uΔ-1; uΔ-2 [calculated $(now())]\n")
					close(BoundariesFile)
		    	
		    	end	
		    		
				if UserSubMode=="StateAnalysis"	|| UserSubMode=="Full"
	    			StateAnalysisFilePathOut = StateAnalysisDirPathOut * "/μ0=$(μ0)_L=$(LL).txt"
	    			
	    			# Prepare file headers
					StateAnalysisFile = open(StateAnalysisFilePathOut,"w")
						write(StateAnalysisFile,"# SFH state properties. μ0=$μ0, nsweeps=$nsweeps, cutoff=$cutoff\n")
						write(StateAnalysisFile,"# L; V; E; δn_M^2; S [calculated $(now())]\n")
					close(StateAnalysisFile)
				end
	    	
	    		# Run horizontal sweep by selected UserSubMode
		    	for L in LL
    	    		println("Starting calculation of observables for L=$L...")
					HorizontalSweep(
						UserSubMode,	# Boundaries/StateAnalysis/Full
						L, 
						VV,
						μ0,
						DMRGParametersXY,
						BoundariesFilePathOut,
						StateAnalysisFilePathOut
					)
				end
				
				if UserSubMode=="Boundaries" || UserSubMode=="Full"
				
					BoundariesFile = open(BoundariesFilePathOut,"a")
						write(BoundariesFile,"# [finished at $(now())]\n")
					close(BoundariesFile)
		    		
				end
				
				if UserSubMode=="StateAnalysis"	|| UserSubMode=="Full"
					
					StateAnalysisFile = open(StateAnalysisFilePathOut,"a")
						write(StateAnalysisFile,"# [finished at $(now())]\n")
					close(StateAnalysisFile)
					
				end
			end
					
			printstyled("\rDone!\n", color=:green)

        # -------------------------- Rectangular sweep -------------------------
        
		elseif UserMode=="--rectangular"

			# Rectangular sweeps
	    	VV = RectangularVV				# Imported from setup
	    	LL = RectangularLL				# Imported from setup
	    	μμ = Rectangularμμ				# Imported from setup
			NN = Int64.(LL./2)				# Half filling

			DirPathOut = PROJECT_ROOT * "/simulations/rectangular-sweep"
    		mkpath(DirPathOut)

            DensityDirPathOut = ""
            FullDirPathOut = ""

            Waiting=true
			print("Choose your method: (Density/Full) ")
			UserSubMode = readline()

            while Waiting
			
				if UserSubMode=="Density"
					
					Waiting=false
					DensityDirPathOut = PROJECT_ROOT * "/simulations/rectangular-sweep/density"
					mkpath(DensityDirPathOut)
					
                elseif UserSubMode=="Full"

                    Waiting=false
					FullDirPathOut = PROJECT_ROOT * "/simulations/rectangular-sweep/full"
					mkpath(FullDirPathOut)

				else
					print("Invalid input. Choose your method: (Density/Full) ")
					UserSubMode = readline()
				end
			
			end

			for (j,L) in enumerate(LL)
				N = NN[j]
				
                if UserSubMode=="Density"
				    FilePathOut = DensityDirPathOut * "/L=$L.txt"
				    DataFile = open(FilePathOut, "w")
					    write(DataFile,"# Fermi-Hubbard model DMRG. L=$L, N=$N\n")
					    write(DataFile,"# V; μ; E; ρ [calculated $(now())]\n")
				    close(DataFile)
                elseif UserSubMode=="Full"
				    FilePathOut = FullDirPathOut * "/L=$L.txt"
				    DataFile = open(FilePathOut, "w")
					    write(DataFile,"# Fermi-Hubbard model DMRG. L=$L, N=$N\n")
					    write(DataFile,"# V; μ; E; ρ; k; D [calculated $(now())]\n")
				    close(DataFile)
                end
				
				println("Starting calculation of observables for L=$L...")
				RectangularSweep(
                    UserSubMode,
					L, 
					N,
					VV,
					μμ,
					DMRGParametersXY,
					FilePathOut
				)
				
				DataFile = open(FilePathOut,"a")
					write(DataFile,"# [finished at $(now())]\n")
				close(DataFile)
			end
			
			printstyled("\rDone!\n", color=:green)
			
		else
			error(ModeErrorMsg)
			exit()
		end
		
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
