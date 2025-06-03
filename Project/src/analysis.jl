#!/usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../FermiHubbardDMRG/src

# Include setup
include(PROJECT_ROOT * "/setup/graphic-setup.jl")
include(PROJECT_ROOT * "/setup/simulations-setup.jl")

# Include modules
include(PROJECT_ROOT * "/modules/methods-plotting.jl")

using DelimitedFiles
using Dates

function main()    

    ModeErrorMsg = "Input error: use option --heatmap, --boundaries or" * 
        "--chain-plots"
	
	if length(ARGS) != 1
		# If user does not specify the user mode
		error(ModeErrorMsg)
		exit()
	else
        UserMode = ARGS[1]

        # ---------------------------- Plot Heatmap ----------------------------

        if UserMode=="--heatmap"

            Full=false
            Waiting=true
			print("Choose your method: (Density/Full) ")
			UserSubMode = readline()
    
            while Waiting
			
				if UserSubMode=="Density"
					
					Waiting=false
					
                elseif UserSubMode=="Full"

                    Waiting=false
                    Full=true

				else
					print("Invalid input. Choose your method: (Density/Full) ")
					UserSubMode = readline()
				end
			
			end

            global HorizontalLL, RectangularLL # Imported from setup
            
			PhaseBoundariesLL = HorizontalLL
            LL = RectangularLL   
            μ0 = 0.0
            Cutoff = 0.05

            PhaseBoundariesFilePath = "" # PROJECT_ROOT * "/../simulations/horizontal_sweep/μ0=$(μ0)_L=$PhaseBoundariesLL.txt"
            HeatmapDir = PROJECT_ROOT * "/../analysis/heatmap/"
            mkpath(HeatmapDir)


			for L in LL
                if !Full
                    FilePathIn = PROJECT_ROOT * "/../simulations/rectangular-sweep/density/L=$(L).txt"#_High.txt"
                elseif Full
                    FilePathIn = PROJECT_ROOT * "/../simulations/rectangular-sweep/full/L=$(L).txt"#_High.txt"
                end

				EFilePathOut = HeatmapDir * "Energy_L=$(L).pdf" 			# Ground-state energy plot
                ρFilePathOut = HeatmapDir * "Density_L=$(L).pdf"       	    # Charge density plot
                
                if !Full
    	            PlotHeatmap(L, FilePathIn; PhaseBoundariesFilePath, EFilePathOut, ρFilePathOut)
                elseif Full
                    kFilePathOut = HeatmapDir * "Compressibility_L=$(L).pdf" 	# Compressibility plot
        	        DFilePathOut = HeatmapDir * "Stiffness_L=$(L).pdf"       	# Charge stiffness plot
        	        uPFilePathOut = HeatmapDir * "uMI_Projection_L=$(L).pdf"	# Unitary filling plot
                	hPFilePathOut = HeatmapDir * "hMI_Projection_L=$(L).pdf"	# Unitary filling plot
    	            PlotHeatmap(
    	            	L,
    	            	FilePathIn;
    	            	PhaseBoundariesFilePath,
    	            	EFilePathOut, 
    	            	kFilePathOut,
    	            	DFilePathOut,
    	            	ρFilePathOut,
    	            	uPFilePathOut,
    	            	hPFilePathOut
    	            )
                end

   			end

        # --------------------- Boundary between XY and MI ---------------------

        elseif UserMode=="--boundaries"
        	
			μ0 = Horizontalμμ[1] # Change!

            FilePathIn = PROJECT_ROOT * "/../simulations/horizontal-sweep/boundaries/μ0=$(μ0)_L=$HorizontalLL.txt"
            PhaseBoundariesDir = PROJECT_ROOT * "/../analysis/phase-boundaries/μ0=$(μ0)/"
            mkpath(PhaseBoundariesDir)

            PlotPhaseBoundaries(
            	FilePathIn;
                FilePathOut = PhaseBoundariesDir * "phase-boundaries_μ0=$(μ0).pdf",
                double=true
            )
            
            # PlotCompressibility..

        elseif UserMode=="--chain-plots"
            
            global StatePropertiesL
            L = StatePropertiesL # Imported from setup
            global XYPoint, IAFPoint, IAFPoint
 
            DirPathIn = PROJECT_ROOT * "/../simulations/states-properties/"
            DirPathOut = PROJECT_ROOT * "/../analysis/states-properties/"

            ChainPlots(
                DirPathIn,
	            DirPathOut,
	            L,
	            XYPoint,
	            IFPoint,
                IAFPoint
            )
		else
			error(ModeErrorMsg)
			exit()
        end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
