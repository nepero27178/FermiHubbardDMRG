#!/usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../BoseHubbardDMRG/src

# Include setup
include(PROJECT_ROOT * "/setup/graphic_setup.jl")
include(PROJECT_ROOT * "/setup/simulations_setup.jl")

# Include modules
include(PROJECT_ROOT * "/modules/dmrg.jl")
include(PROJECT_ROOT * "/modules/plots.jl")
include(PROJECT_ROOT * "/modules/sweeps.jl")
include(PROJECT_ROOT * "/modules/subdomain_selection.jl")

# Move back
PROJECT_ROOT *= "/.."	# Absloute path up to .../BoseHubbardDMRG/

function main()

	"""
	Note: the main() functions at this level are the only functions in the whole
	project referring to externally defined global variables. All simulations
	variables have been defined in ./setup/simulations_setup.jl and are here
	referred to with their names. Check the setup file to modify the
	definitions.
	"""
    
	ModeErrorMsg = "Input error: use option --boundaries, --horizontal, " *
				   "--rectangular or --rectangular-selection"
	
	if length(ARGS) != 1
		# If user does not specify the user mode
		error(ModeErrorMsg)
		exit()
	else
		
		UserMode = ARGS[1]

        # -------------------------- Boundaries sweep --------------------------
        
        if UserMode=="--boundaries"

			# Horizontal sweeps
			LL = HorizontalLL				# Imported from setup
	    	JJ = HorizontalJJ				# Imported from setup
	    	μ0 = 0.0						# Imported from setup
	    	
	    	DirPathOut = PROJECT_ROOT * "/simulations/boundaries_sweep"
    		mkpath(DirPathOut)
	    		
    		FilePathOut = DirPathOut * "/μ0=$(μ0)_L=$(LL).txt"
    	
			DataFile = open(FilePathOut,"w")
			write(DataFile,"# Hubbard model DMRG. This file contains many sizes. nmax=$nmax, μ0=$μ0, nsweeps=$nsweeps, cutoff=$cutoff\n")
			write(DataFile,"# L; J; E; μ+; μ- [calculated $(now())]\n")
			close(DataFile)
    	
	    	for L in LL
	    		println("Starting calculation of observables for L=$L...")
	    		
	    		# Note: here we use superfluid DMRG parameters everywhere
				BoundariesSweep(L, nmax, JJ, DMRGParametersSF, FilePathOut)
			end
			
			DataFile = open(FilePathOut,"a")
			write(DataFile,"# [finished $(now())]\n")
			close(DataFile)
					
			println("Done!")
        
        # -------------------------- Horizontal sweep --------------------------
        
        elseif UserMode=="--horizontal"

			# Horizontal sweeps
			LL = HorizontalLL				# Imported from setup
	    	JJ = HorizontalJJ				# Imported from setup
	    	μ0μ0 = Horizontalμμ				# Imported from setup
	    	
	    	DirPathOut = PROJECT_ROOT * "/simulations/horizontal_sweep"
    		mkpath(DirPathOut)
	    	
	    	for μ0 in μ0μ0
	    		
	    		FilePathOut = DirPathOut * "/μ0=$(μ0)_L=$(LL).txt"
	    	
				DataFile = open(FilePathOut,"w")
				write(DataFile,"# Hubbard model DMRG. This file contains many sizes. nmax=$nmax, μ0=$μ0, nsweeps=$nsweeps, cutoff=$cutoff\n")
				write(DataFile,"# L; J; E; Γ; eΓ; [calculated $(now())]\n")
				close(DataFile)
	    	
		    	for L in LL
    	    		println("Starting calculation of observables for L=$L...")
    	    		
    	    		# Note: here we use superfluid DMRG parameters everywhere
					HorizontalSweep(L, nmax, JJ, μ0, DMRGParametersSF, FilePathOut)
				end
				
				DataFile = open(FilePathOut,"a")
				write(DataFile,"# [finished at $(now())]\n")
				close(DataFile)
				
			end
					
			println("Done!")

        # -------------------------- Rectangular sweep -------------------------
        
		elseif UserMode=="--rectangular"
		
			# Rectangular sweeps
	    	JJ = RectangularJJ				# Imported from setup
	    	LL = RectangularLL				# Imported from setup
	    	μμ = Rectangularμμ				# Imported from setup

			μ0 = 0.0 						# μ0 from which to take phase boundaries

			NN = LL							# Unitary filling

			DirPathOut = PROJECT_ROOT * "/simulations/rectangular_sweep"
    		mkpath(DirPathOut)

			# Uncomment here to use L-wise phase boundaries
			# FilePathIn =  PROJECT_ROOT * "/simulations/horizontal_sweep/L=$LL.txt"
			FilePathIn =  PROJECT_ROOT * "/analysis/phase_boundaries/μ0=$μ0/fitted_phase_boundaries_μ0=$μ0.txt"

			for (j,L) in enumerate(LL)
				N = NN[j]
				
				FilePathOut = DirPathOut * "/L=$L.txt"
				DataFile = open(FilePathOut, "w")
				write(DataFile,"# Hubbard model DMRG. L=$L, N=$N, nmax=$nmax\n")
				write(DataFile,"# NOTE: Different DMRG settings have been used for MI and SF. Calculation performed on central site.\n")
				write(DataFile,"# J; μ; E; n_variance; <a>; k[calculated $(now())]\n")
				close(DataFile)
				
				println("Starting calculation of observables for L=$L...")
				RectangularSweep(L, N, nmax, JJ, μμ, DMRGParametersMI, DMRGParametersSF, FilePathIn, FilePathOut)
				
				DataFile = open(FilePathOut,"a")
				write(DataFile,"# [finished at $(now())]\n")
				close(DataFile)
			end
			
			println("Done!")

        # --------------------- Rectangular selection sweep --------------------
        
        elseif UserMode=="--rectangular-selection"
		
			# Rectangular selection sweep
			LL = RectangularSelectionLL		# Imported from setup
			
			NN = LL							# Unitary filling
			
			μ0 = 0.0 						# μ0 from which to take phase boundaries
			
			println("Starting tip selection...")
			FilePathIn = PROJECT_ROOT * "/analysis/phase_boundaries/μ0=$μ0/fitted_phase_boundaries.txt"
			Selections = FindMottTip(FilePathIn; verbose=false)
			PlotSelection(FilePathIn, Selections)
			
			print("Should I run a rectangular simulation inside the selected region? (y/n) ")
			PerformTipSweep = readline()
			while true
				if PerformTipSweep == "y"
					println("Selection accepted. Starting simulations around Mott lobe tip...")
					break
				elseif PerformTipSweep == "n"
					println("Selection rejected. Change the setting of this program to improve selection.")
					exit()
				else
					print("Invalid input. Please use valid input. (y/n) ")
					PerformTipSweep = readline()
				end
			end
			
			if size(Selections, 1) > 1
				print("There are $(size(Selections,1)) intersection points. Which one do you select? ")
				UserTipSelection = parse(Int64, readline())
				if UserTipSelection>0 & UserTipSelection < size(Selections,1)
					Left, Right, Up, Down = Selections[UserTipSelection,:]
				else
					print("Invalid input. Please use valid input. (1, ..., $(size(Selections,1)))")
				end
			else
				Left, Right, Up, Down = Selections
			end

			JJ = collect(range(start=Left, stop=Right, length=Resolution["Medium"]))
			μμ = collect(range(start=Down, stop=Up, length=Resolution["Medium"]))

			DirPathOut = PROJECT_ROOT * "/simulations/rectangular_selection_sweep"
    		mkpath(DirPathOut)

			# Uncomment here to use L-wise phase boundaries
			# FilePathIn =  PROJECT_ROOT * "/simulations/horizontal_sweep/L=$LL.txt"
			FilePathIn =  PROJECT_ROOT * "/analysis/phase_boundaries/μ0=$μ0/fitted_phase_boundaries.txt"

			for (j,L) in enumerate(LL)
				N = NN[j]
				
				FilePathOut = DirPathOut * "/L=$L.txt"
				DataFile = open(FilePathOut, "w")
				write(DataFile,"# Hubbard model DMRG. L=$L, N=$N, nmax=$nmax\n")
				write(DataFile,"# NOTE: Different DMRG settings have been used for MI and SF. Calculation performed on central site.\n")
				write(DataFile,"# J; μ; E; n_variance; <a> [calculated $(now())]\n")
				close(DataFile)
				
				println("Starting calculation of observables for L=$L...")
				RectangularSweep(L, N, nmax, JJ, μμ, DMRGParametersMI, DMRGParametersSF, FilePathIn, FilePathOut)
				
				DataFile = open(FilePathOut,"a")
				write(DataFile,"# [finished at $(now())]\n")
				close(DataFile)
			end
			
			println("Done!")
			
		else
			error(ModeErrorMsg)
			exit()
		end
		
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
