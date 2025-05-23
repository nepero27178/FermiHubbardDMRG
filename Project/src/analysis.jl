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

    ModeErrorMsg = "Input error: use option --heatmap, --boundaries, --gamma, gamma-K
    or --gamma-MISF"
	
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
            L = 10             
            μ0 = 0.0
            Cutoff = 0.05

            PhaseBoundariesFilePath = "" # PROJECT_ROOT * "/../simulations/horizontal_sweep/μ0=$(μ0)_L=$PhaseBoundariesLL.txt"
            HeatmapDir = PROJECT_ROOT * "/../analysis/heatmap/"
            mkpath(HeatmapDir)


			for L in LL
                FilePathIn = PROJECT_ROOT * "/../simulations/rectangular-sweep/L=$(L).txt"#_High.txt"

				EFilePathOut = HeatmapDir * "Energy_L=$(L)_High-Resolution.pdf" 			# Ground-state energy plot
	            kFilePathOut = HeatmapDir * "Compressibility_L=$(L)_High-Resolution.pdf" 	# Compressibility plot
    	        DFilePathOut = HeatmapDir * "Stiffness_L=$(L)_High-Resolution.pdf"       	# Charge stiffness plot
				
	            PlotHeatmap(L, FilePathIn; PhaseBoundariesFilePath, EFilePathOut, kFilePathOut, DFilePathOut)
   			end
            
#            FilePathIn = PROJECT_ROOT * "/../simulations/rectangular-sweep/L=$(L)_High.txt"

#            EstimateFractionOrderParameter(L, FilePathIn, PhaseBoundariesFilePath, Cutoff)


        # ----------------------------------------------------------------------
        # --------------------- Boundary between SF and MI ---------------------
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

        # ----------------------------------------------------------------------
        # ----------------------- Correlation function Γ -----------------------
        # ----------------------------------------------------------------------

        elseif UserMode=="--gamma"
    
            @warn "Mode under construction."
        
#            global HorizontalLL, RectangularLL # Imported from setup
#
#            μ0 = 0.08 # CHANGE!
#            rrMin = [2, 4, 6, 8]
#            rrMax = [12, 14, 16, 18]
#            JMin = 0.20
#            LMin = 20
#
#            GammaDir = PROJECT_ROOT * "/../analysis/gamma/μ0=$(μ0)/"
#            mkpath(GammaDir)
#
#            FilePathIn = PROJECT_ROOT * "/../simulations/horizontal_sweep/μ0=$(μ0)_L=$HorizontalLL.txt"
#            FilePathFit = GammaDir * "fitted_Luttinger_parameter_μ0=$μ0.txt"
#
#            # (Step 1) Perform all the fits for all J>J_min
#            FitRoutineGamma(FilePathIn, FilePathFit; rrMin, rrMax, JMin, LMin)
#
#            # (Step 2) Plot Γ(r) vs r for one  given (J, μ0) ( j ∈ [1, 50] )
#            # j = 43 corresponds to J=0.30
#            # ! NOTE: already plotted for all the μ0 we have. !
#            # for j in 12:2:50
#            #     mkpath(GammaDir*"data_plot/")
#            #     FileGammaPlot = GammaDir * "data_plot/data_gamma_j=$(j)_μ0=$μ0.pdf"
#            #     PlotPowerLawGamma(FilePathIn, μ0, j; FilePathOut=FileGammaPlot, overwrite=false)
#            # end
#
#            # (Step 3) Plot K_∞ vs J, reading data from Step 1
#            FilePathInK = GammaDir * "fitted_Luttinger_parameter_μ0=$μ0.txt"
#            for rMin in [2,4,6,8]
#                FilePathOutKInfty = GammaDir * "fitted_Luttinger_plot_μ0=$(μ0)_rMin=$rMin.pdf"
#                PlotFitResultsK(FilePathInK, μ0, rMin; FilePathOut=FilePathOutKInfty)
#            end
#        
#        elseif UserMode=="--gamma-MISF"
#            global HorizontalLL
#            μ0 = 0.6
#
#            FilePathIn = PROJECT_ROOT * "/../simulations/horizontal_sweep/μ0=$(μ0)_L=$HorizontalLL.txt"
#            PhaseBoundariesFilePath = PROJECT_ROOT * "/../simulations/boundaries_sweep/μ0=0.0_L=$HorizontalLL.txt"
#
#            GammaDir = PROJECT_ROOT * "/../analysis/gamma/illustrative_plots/"
#            mkpath(GammaDir)
#
#            FilePathOut1 = GammaDir * "gamma_MI_vs_SF.pdf"
#            FilePathOut2 = GammaDir * "gamma_MI_vs_SF_phaseboundaries.pdf"
#            PlotCorrelationFunctionsMIvsSF(FilePathIn; μ0, JMin=0.10, JMax=0.15, L=70, PhaseBoundariesFilePath,
#                FilePathOut1, FilePathOut2)
        
        elseif UserMode=="--gamma-K"

            @warn "Mode under construction."

#            # Make an illustrative plot of K_∞ vs J, for one
#            # rMin and one rMax, chosen with the definitions below.
#
#            global HorizontalLL
#            μ0 = 0.0
#            rMin = 8 
#            rMax = 18
#
#            GammaDir = PROJECT_ROOT * "/../analysis/gamma/μ0=$(μ0)/"
#            FilePathInK = GammaDir * "fitted_Luttinger_parameter_μ0=$μ0.txt"
#            FilePathOutKInfty = GammaDir * "illustrative_K_μ0=$μ0.pdf"
#            PlotIllustrativeResultK(FilePathInK, μ0, rMin, rMax; FilePathOut=FilePathOutKInfty)
#
		else
			error(ModeErrorMsg)
			exit()
        end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
