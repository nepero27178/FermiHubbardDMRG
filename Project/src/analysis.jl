#!/usr/bin/julia

PROJECT_ROOT = @__DIR__ # Absloute path up to .../FermiHubbardDMRG/src

# Include setup
include(PROJECT_ROOT * "/setup/graphic-setup.jl")
include(PROJECT_ROOT * "/setup/simulations-setup.jl")

# Include modules
include(PROJECT_ROOT * "/modules/methods-fitting.jl")
include(PROJECT_ROOT * "/modules/methods-plotting.jl")
include(PROJECT_ROOT * "/modules/subdomain-selector.jl")

using DelimitedFiles
using Dates

function main()    

    ModeErrorMsg = "Input error: use option --heatmap, " * 
    "--boundaries, --zero-field"
	
	if length(ARGS) != 1
		# If user does not specify the user mode
		error(ModeErrorMsg)
		exit()
	else
        UserMode = ARGS[1]

        # ----------------------------------------------------------------------
        # ---------------------------- Plot Heatmap ----------------------------
        # ----------------------------------------------------------------------

        if UserMode=="--heatmap"

            global HorizontalLL, RectangularLL # Imported from setup
            
			PhaseBoundariesLL = HorizontalLL
            LL = RectangularLL   
            μ0 = 0.0
            Cutoff = 0.05

            PhaseBoundariesFilePath = "" # PROJECT_ROOT * "/../simulations/horizontal_sweep/μ0=$(μ0)_L=$PhaseBoundariesLL.txt"
            HeatmapDir = PROJECT_ROOT * "/../analysis/heatmap/"
            mkpath(HeatmapDir)


			for L in LL
                FilePathIn = PROJECT_ROOT * "/../simulations/rectangular-sweep/L=$(L).txt"#_High.txt"

				EFilePathOut = HeatmapDir * "Energy_L=$(L).pdf" 			# Ground-state energy plot
	            kFilePathOut = HeatmapDir * "Compressibility_L=$(L).pdf" 	# Compressibility plot
    	        DFilePathOut = HeatmapDir * "Stiffness_L=$(L).pdf"       	# Charge stiffness plot
		        ρFilePathOut = ZeroFieldDir * "Density_L=$(L).pdf"       	# Charge density plot
				
	            PlotHeatmap(L, FilePathIn; PhaseBoundariesFilePath, EFilePathOut, kFilePathOut, DFilePathOut)
   			end
            
#            FilePathIn = PROJECT_ROOT * "/../simulations/rectangular-sweep/L=$(L)_High.txt"

#            EstimateFractionOrderParameter(L, FilePathIn, PhaseBoundariesFilePath, Cutoff)

		# ----------------------------------------------------------------------
        # ------------------------- Conserve particles -------------------------
        # ----------------------------------------------------------------------

        elseif UserMode=="--zero-field"

            global HorizontalLL, RectangularLL # Imported from setup
            
			PhaseBoundariesLL = HorizontalLL
            LL = RectangularLL          
            μ0 = 0.0
            Cutoff = 0.05

            PhaseBoundariesFilePath = "" # PROJECT_ROOT * "/../simulations/horizontal_sweep/μ0=$(μ0)_L=$PhaseBoundariesLL.txt"
            ZeroFieldDir = PROJECT_ROOT * "/../analysis/zero-field/"
            mkpath(HeatmapDir)

            FilePathIn = PROJECT_ROOT * "/../simulations/zero-field-sweep/μ0=$(μ0)_L=$(LL).txt"#_High.txt"

			@warn "Go on from: analysis.jl, Line 80"

			EFilePathOut = ZeroFieldDir * "Energy_L=$(L).pdf" 			# Ground-state energy plot
            kFilePathOut = ZeroFieldDir * "Compressibility_L=$(L).pdf" 	# Compressibility plot
	        DFilePathOut = ZeroFieldDir * "Stiffness_L=$(L).pdf"       	# Charge stiffness plot
			
            PlotZeroFieldSweep(L, FilePathIn; EFilePathOut, kFilePathOut, DFilePathOut)

        # ----------------------------------------------------------------------
        # --------------------- Boundary between XY and MI ---------------------
        # ----------------------------------------------------------------------

        elseif UserMode=="--boundaries"
        	
			μ0 = 0.0 #Horizontalμμ[3] # CHANGE!

            FilePathIn = PROJECT_ROOT * "/../simulations/boundaries-sweep/μ0=$(μ0)_L=$HorizontalLL.txt"
            PhaseBoundariesDir = PROJECT_ROOT * "/../analysis/phase-boundaries/μ0=$(μ0)/"
            mkpath(PhaseBoundariesDir)

            FilePathPlot = PhaseBoundariesDir * "phaseboundaries_μ0=$(μ0).pdf"
            FilePathFit = PhaseBoundariesDir * "fitted-phase-boundaries_μ0=$(μ0).txt"
            
            FilePathPlotOut = PhaseBoundariesDir * "phaseboundaries-fit_μ0=$(μ0).pdf"
            FilePathSinglePlotOut = PhaseBoundariesDir * "phaseboundaries-fit-single_μ0=$(μ0).pdf"

            # Uncomment the needed analysis.
            #PlotPhaseBoundaries(FilePathIn; gap=false, FilePathOut=FilePathPlot, μ0)
            # FitPhaseBoundaries(FilePathIn, FilePathFit; FilePathPlotOut, FilePathSinglePlotOut, μ0)
            PlotPhaseBoundaries(
            	FilePathIn;
                FilePathOut=PhaseBoundariesDir*"colored_PB.pdf",
                HideDataPoints=false,
                DrawMottLobe=false,
                MottLobeFilePath=PROJECT_ROOT * "/../analysis/phase-boundaries/μ0=$μ0/fitted-phase-boundaries_μ0=$μ0.txt",
            )

            # FindMottTip(FilePathFit, verbose=true)

		else
			error(ModeErrorMsg)
			exit()
        end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
