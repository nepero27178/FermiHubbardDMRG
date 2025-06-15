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

			Density=false
			Complementary=false
            Full=false
            
            Waiting=true
			print("Choose your method: (Density/Complementary/Full) ")
			UserSubMode = readline()
    
            while Waiting
			
				if UserSubMode=="Density"
					
					Waiting=false
					Density=true

                elseif UserSubMode=="Complementary"

                    Waiting=false
                    Complementary=true
					
                elseif UserSubMode=="Full"

                    Waiting=false
                    Full=true

				else
					print("Invalid input. Choose your method: (Density/Complementary/Full) ")
					UserSubMode = readline()
				end
			
			end

            global HorizontalLL, RectangularLL # Imported from setup
            
			PhaseBoundariesLL = HorizontalLL
            LL = RectangularLL   
            μ0 = 0.0
            Cutoff = 0.05

            PhaseBoundariesFilePath = "" # Add if needed
            HeatmapDir = PROJECT_ROOT * "/../analysis/heatmap/"
            mkpath(HeatmapDir)

			for L in LL
			
                if Density
                    FilePathIn = PROJECT_ROOT * "/../simulations/rectangular-sweep/density/L=$(L).txt"#_High.txt"
                elseif Complementary
                    FilePathIn = PROJECT_ROOT * "/../simulations/rectangular-sweep/complementary/L=$(L).txt"#_High.txt"
                elseif Full
                    FilePathIn = PROJECT_ROOT * "/../simulations/rectangular-sweep/full/L=$(L).txt"#_High.txt"
                end

				EFilePathOut = HeatmapDir * "energy_L=$(L).pdf" 				# Ground-state energy plot
                ρFilePathOut = HeatmapDir * "density_L=$(L).pdf"       	    	# Charge density plot
                δFilePathOut = HeatmapDir * "block-density-variance_L=$(L).pdf" # Block density variance plot
    	        uPFilePathOut = HeatmapDir * "uMI_projection_L=$(L).pdf"		# Unitary filling plot
            	hPFilePathOut = HeatmapDir * "hMI_projection_L=$(L).pdf"		# Unitary filling plot
            	kFilePathOut = HeatmapDir * "compressibility_L=$(L).pdf" 		# Compressibility plot
    	        DFilePathOut = HeatmapDir * "stiffness_L=$(L).pdf"       		# Charge stiffness plot
    	        KFilePathOut = HeatmapDir * "K-Luttinger_L=$(L).pdf"       		# K Luttinger parameter plot
    	        uFilePathOut = HeatmapDir * "u-Luttinger_L=$(L).pdf"       		# u Luttinger parameter plot
                
                if Density
                
    	            PlotHeatmap(
    	            	L,
    	            	FilePathIn;
    	            	PhaseBoundariesFilePath,
    	            	EFilePathOut,
    	            	ρFilePathOut,
    	            	δFilePathOut
    	            )
    	            
    	        elseif Complementary

    	            PlotHeatmap(
    	            	L,
    	            	FilePathIn;
    	            	PhaseBoundariesFilePath,
    	            	EFilePathOut, 
    	            	uPFilePathOut,
    	            	hPFilePathOut,
    	            	kFilePathOut,
    	            	DFilePathOut,
    	            	KFilePathOut,
    	            	uFilePathOut
    	            )
    	            
                elseif Full

    	            PlotHeatmap(
    	            	L,
    	            	FilePathIn;
    	            	PhaseBoundariesFilePath,
    	            	EFilePathOut,
    	            	ρFilePathOut,
    	            	δFilePathOut,
    	            	uPFilePathOut,
    	            	hPFilePathOut,
    	            	kFilePathOut,
    	            	DFilePathOut,
    	            	KFilePathOut,
    	            	uFilePathOut
    	            )
    	            
                end

   			end

        # --------------------- Boundary between XY and MI ---------------------

        elseif UserMode=="--boundaries"
        	
			μ0 = Horizontalμμ[1] # Change!

            FilePathIn = PROJECT_ROOT * "/../simulations/horizontal-sweep/boundaries/μ0=$(μ0)_L=$HorizontalLL.txt"
            PhaseBoundariesDir = PROJECT_ROOT * "/../analysis/phase-boundaries/"
            mkpath(PhaseBoundariesDir)

            PlotPhaseBoundaries(
            	FilePathIn;
                FilePathOut = PhaseBoundariesDir * "phase-boundaries_μ0=$(μ0)_L=$HorizontalLL.pdf",
                double=true
            )
            
            # PlotCompressibility..

        elseif UserMode=="--chain-plots"
            
            global StatePropertiesL
            L = StatePropertiesL # Imported from setup
            global XYPoint1, XYPoint2, IAFPoint, IAFPoint
 
            DirPathIn = PROJECT_ROOT * "/../simulations/states-properties/"
            DirPathOut = PROJECT_ROOT * "/../analysis/states-properties/"

			for CompactPlot in [false] #, true]
		        ChainPlots(
		            DirPathIn,
			        DirPathOut,
			        L,
			        XYPoint1,
			        XYPoint2,
			        IFPoint,
		            IAFPoint;
		            CompactPlot
		        )
		    end
		else
			error(ModeErrorMsg)
			exit()
        end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
