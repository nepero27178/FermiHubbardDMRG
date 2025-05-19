#!/usr/bin/julia

# ---------------------------------- Heatmap -----------------------------------

"""
Plot the variance of the number of particles from data saved
in `FilePathIn`.
"""
function PlotHeatmap(
        L::Int64,
        FilePathIn::String;
        PhaseBoundariesFilePath="",
        VarianceFilePathOut="",		    # Variance heatmap
        DFilePathOut="",				# Charge stiffness
        kFilePathOut=""                 # Compressibility heatmap
    )

    @warn "Mode under construction."
    
#    # Extract data coming from rectangular_sweep
#    Data = readdlm(FilePathIn, ';', '\n'; comments=true)
#    # TODO not uniform ; vs , everywhere
#
#    # display(Data)
#
#    JJ = Data[:,1]
#    μμ = Data[:,2]
#    EE = Data[:,3]
#    varvar = Data[:,4]
#    aa = Data[:,5]
#    kk = Data[:,6]
#
#    NumJ = length(unique(JJ))
#    Numμ = length(unique(μμ))
#
#    # Store variance and order parameter <a_i>
#    Variances = zeros(Numμ, NumJ)
#    OrderParameters = zeros(Numμ, NumJ)
#    FourierTransforms = zeros(Numμ, NumJ)
#    Compressibilities = zeros(Numμ, NumJ)
#
#    for jj in 1:NumJ
#        Variances[:,jj] = varvar[ Numμ*(jj-1)+1 : Numμ*jj ]
#        OrderParameters[:,jj] = aa[ Numμ*(jj-1)+1 : Numμ*jj ]
#        Compressibilities[:,jj] = kk[ Numμ*(jj-1)+1 : Numμ*jj ]	# First row is NaN!
#    end
#
#    i = (ceil(Int64, L/2)) # site index
#
#    if VarianceFilePathOut != ""
#        # Plot variance
#        heatmap(unique(JJ), unique(μμ), Variances, 
#                xlabel=L"J",
#                ylabel=L"$\mu$",
#                title=L"Variance $\delta n_i^2$ ($L=%$L, i=%$i$)")
#        ylabel!(L"μ")
#
#        if PhaseBoundariesFilePath != ""
#            HeatmapAddPhaseBoundaries(PhaseBoundariesFilePath, L, JJ, μμ)
#        end
#
#        savefig(VarianceFilePathOut)
#        println("Variance plot for L=$L saved on file!")
#    end
#
#    if AFilePathOut != ""
#        # Plot order parameter <a_i>
#        heatmap(unique(JJ), unique(μμ), OrderParameters, 
#                xlabel=L"J",
#                ylabel=L"$\mu$",
#                title=L"$\langle \hat b_i \rangle$ ($L=%$L, i=%$i$)")
#
#        # Add phase boundaries
#        if PhaseBoundariesFilePath != ""
#            HeatmapAddPhaseBoundaries(PhaseBoundariesFilePath, L, JJ, μμ)
#        end
#
#        savefig(AFilePathOut)
#        println("Order parameter plot for L=$L saved on file!")
#    end
#	
#    if KFilePathOut != ""
#		# Plot compressibility
#        heatmap(unique(JJ), unique(μμ)[3:end], Compressibilities[3:end,:], 
#                xlabel=L"J",
#                ylabel=L"$\mu$",
#                title=L"Compressibility $\kappa$ ($L=%$L, i=%$i$)",
#                clim=(0,5))
#                
#        # Add phase boundaries
#        if PhaseBoundariesFilePath != ""
#            HeatmapAddPhaseBoundaries(PhaseBoundariesFilePath, L, JJ, μμ)
#        end
#        
#        savefig(KFilePathOut)
#        println("Compressibility plot for L=$L saved on file!")        
#    end
#end
#
#"""
#Add the phase boundaries to the heatmap, for the closest size L available,
#and adjusting the xlimits and ylimits appropriately.
#"""
#function HeatmapAddPhaseBoundaries(PhaseBoundariesFilePath::String,
#                                   L::Int64,
#                                   JJ::Array{Float64},
#                                   μμ::Array{Float64};
#                                   μ0=0.0)
#                                   
#    BoundariesData = readdlm(PhaseBoundariesFilePath, ';', '\n'; comments=true)
#    LL = unique(BoundariesData[:, 1])
#
#    L_PB_index = argmin(abs.(LL .- L)) # best approximation of L
#    L_PhaseBoundaries = LL[L_PB_index]
#
#    if L_PhaseBoundaries !== L
#        println("L=$L is not among horizontal data. Using the closest available size, L=$L_PhaseBoundaries. (μ0=$μ0)")
#    end
#
#    # Filter data corresponding to L_PhaseBoundaries
#    indices = (BoundariesData[:, 1] .== L_PhaseBoundaries)
#    JJ_PB = BoundariesData[indices, 2]
#    μUp = BoundariesData[indices, 4]
#    μDown = BoundariesData[indices, 5]
#    
#    plot!(JJ_PB, -μDown .+ 2 * μ0, 
#        label=nothing,#"$L=%$L_PhaseBoundaries$",
#        seriestype=:scatter,
#        markersize=1.5,
#        color="white",
#        xlimits=(minimum(JJ), maximum(JJ)),
#        ylimits=(minimum(μμ), maximum(μμ)))
#    plot!(JJ_PB, μ0 .+ μUp, seriestype=:scatter,
#        label=nothing,
#        markersize=1.5,
#        color="white",
#        xlimits=(minimum(JJ), maximum(JJ)),
#        ylimits=(minimum(μμ), maximum(μμ)))      
end

# ----------------------------- Phase boundaries -------------------------------

"""
Plot the phase boundaries between the Mott Insulator (MI) and Superfluid (SF) 
phases, calculated from RectangularSweep. If `gap=true`, plot the charge gap
instead of the phase boundaries. If `overwrite=true`, clear previous plots. 
If `CustomLL` specified, plot only those sizes.
"""
function PlotPhaseBoundaries(
        FilePathIn::String; 
        FilePathOut="",
        gap=false,
        overwrite=true, 
        CustomLL=[],
        μ0=0.0,
        HideDataPoints=false,
        DrawMottLobe=false,
        MottLobeFilePath=PROJECT_ROOT * "/analysis/phase_boundaries/fitted-phase-boundaries.txt",
    )
    
    global HorizontalLL # Imported from setup
    BoundariesData = readdlm(FilePathIn, ';', '\n'; comments=true)
    
    # gr(legend_background_color=RGBA{Float64}(1, 1, 1, 0.8))

    if overwrite
        plot()
    end

    # Extract unique L values
    if CustomLL==[]
        LL = HorizontalLL # unique(BoundariesData[:, 1])
    else
        LL = CustomLL
    end

	tt = 0	# Initialize tt to use it outside the for loop
	if !HideDataPoints
		for (l, L) in enumerate(LL)
		    indices = (BoundariesData[:, 1] .== L)
		    tt = BoundariesData[indices,2]
		    μUp = BoundariesData[indices,4]
		    μDown = BoundariesData[indices,5]
		    
		    if gap
		        plot!(
                    tt,                 μUp - μDown, 
                    xlabel=L"$J$",      ylabel=L"$\Delta E_{\mathrm{gap}}$", 
                    title=L"Charge gap  as a function of $t$ ($\mu_0=%$μ0$)",
                    seriestype=:scatter,
                    markersize=1.5,
                    label=L"$L=%$L$",
                    color=MyColors[l % length(MyColors)]
                )
		    else
		        plot!(
                    tt,                 -μDown .+ 2 * μ0, 
                    xlabel=L"$J$",      ylabel=L"$\mu$",
                    label=L"$L=%$L$",
                    title=L"Extrapolation of $\mu_c^\pm(t)$ ($\mu_0=%$μ0$)",
                    seriestype=:scatter,
                    markersize=1.5,
                    color=MyColors[l % length(MyColors)]
                )

		        plot!(
                    tt,                 μUp, 
                    seriestype=:scatter,
                    label="",
                    markersize=1.5,
                    color=MyColors[l % length(MyColors)]
                )
		    end
		end
	elseif HideDataPoints
		plot()
    end
    
    if DrawMottLobe

        @warn "Mode under construction."

#    	MottLobeData = readdlm(MottLobeFilePath, ',', '\n'; comments=true)
#    	JJ = MottLobeData[:,1]
#    	ΔEp = MottLobeData[:,2]
#    	ΔEm = MottLobeData[:,3]
#    	
#    	plot!(JJ, [zeros(length(JJ)) ΔEp],
#			  linewidth=0,
#			  fillrange=[-ΔEm ones(length(JJ))],
#			  fillcolor="red",
#			  fillalpha=0.1,
#			  label=nothing,
#              xlimits=(minimum(JJ), maximum(JJ)),
#              ylimits=(minimum(-ΔEm), maximum(ΔEp)),
#              xlabel=L"$J$",
#              ylabel=L"$\mu$",
#              title="Phases of the model")
#    	
#    	plot!(JJ, -ΔEm,
#			  linewidth=0,
#			  fillrange=ΔEp,
#			  fillcolor="blue",
#			  fillalpha=0.1,
#			  label=nothing)
#			  
#        # MI / SF text
#		annotate!(0.07, 0.35, text(L"MI ($\rho=1$)", 10))
#        annotate!(0.18, 0.62, text("SF", 10))
#        annotate!(0.18, 0.09, text("SF", 10))
#
#        # Phase boundaries
#        plot!(JJ, [ΔEp, -ΔEm],
#              label=[L"\mu_c^+ \, (L \rightarrow \infty)" L"\mu_c^- \, (L \rightarrow \infty)"],
#              color=["black" "black"],
#              linestyle=[:dash :dashdot])
#
#        # Intersection point (Mott tip)
#        scatter!([0.308], [0.080], markersize=1.5, xerr=[0.004], yerr=[0.008],
#            label="Intersection", color="blue", msw=0.5, markerstrokecolor="blue")

    end
    
    if !=(FilePathOut,"")
        savefig(FilePathOut)
        if gap
            println("\nGap for L=$(Int.(LL)) plotted to ", FilePathOut)
        else
            println("\nPhase boundaries for L=$(Int.(LL)) plotted to ", FilePathOut)
        end
    else
    	gui()
    end
end			 

"""
Plot the phase boundaries and a specific point.
"""
function PlotPointAndPhaseBoundaries(MottLobeFilePath::String)

    @warn "Mode under construction."

#    MottLobeData = readdlm(MottLobeFilePath, ',', '\n'; comments=true)
#    JJ = MottLobeData[:,1]
#    ΔEp = MottLobeData[:,2]
#    ΔEm = MottLobeData[:,3]
#
#    plot(size=(352, 256))
#
#    # Phase boundaries
#    plot!(JJ, [ΔEp, -ΔEm],
#            #label=[L"\mu_c^+ \, (L \rightarrow \infty)" L"\mu_c^- \, (L \rightarrow \infty)"],
#            label=nothing,
#            color=["black" "black"],
#            xlabel=L"$J$",
#            ylabel=L"$\mu$",
#            background_color = :transparent,
#            minorticks=false)
#            #linestyle=[:dash :dashdot])
#
#    scatter!([0.06], [0.4], markersize=4, marker=:diamond,
#        label=nothing, color=MyColors[4])
#
#    scatter!([0.33], [0.8], markersize=4, marker=:diamond,
#       label=nothing, color=MyColors[4])
#    
#    # MI / SF text
#    annotate!(0.06, 0.48, text("MI", 9))
#    annotate!(0.31, 0.8, text("SF", 9))
#
#    gui()
end

# ------------------------------------------------------------------------------
# --------------------------- Correlation functions ----------------------------
# ------------------------------------------------------------------------------

""" Plot the correlation function Γ(r) for the chosen J (the j-th) and μ0."""
function PlotPowerLawGamma(FilePathIn::String,
                           μ0::Float64,
                           j::Int64;
                           FilePathOut="",
                           overwrite=true)

    @warn "Mode under construction."

#    plot()
#
#    # Read the input data
#    data = readdlm(FilePathIn, ';', '\n'; comments=true)
#
#    # Extract unique J, L values
#    LL = unique(data[:, 1])
#    JJ = unique(data[:, 2])
#
#    println("\nPlotting correlation function.")
#    println("From input file, there are $(length(JJ)) possible values of J.")
#    
#    # Mastruzzo to extract array of Γ
#    function parse_array(str)
#        return parse.(Float64, split(strip(str, ['[', ']', ' ']), ','))
#    end
#    Γall = [parse_array(row[4]) for row in eachrow(data)]
#
#    if overwrite
#        plot()
#    end
#
#    J = JJ[j] # choose the j-th J
#
#    println("Chosen value of J: $J, chosen value of μ0: $μ0")
#
#    for L in LL
#        # Filter data for the current J and L
#        filter = (data[:, 2] .== J) .& (data[:, 1] .== L)
#
#        # Extract the correlators {Γ(r,J,L)}_r for the selected J,L
#        Γeven = Γall[filter][1] # this is an array, Γ(r even)
#        r = range(start=2, step=2, length=length(Γeven)) # r even
#
#        # Check if Γeven has at least two elements
#        if length(Γeven) < 2
#            println("Warning: Not enough data points for L = $L. Skipping.")
#            continue
#        end
#
#        J_round = round(J, digits=3)
#
#        scatter!(r, Γeven,
#            xlabel=L"$r$",
#            ylabel=L"$\Gamma(r)$",
#            title=L"Correlation function ($J=%$J_round, \mu_0 = %$μ0$)",
#            label=L"L=%$L",
#            markersize=2,
#            xscale=:log10,
#            yscale=:log10,  
#            xticks=[1,10,100],
#            xlimits=(1,100),
#            legend=:topright)
#    end
#
#    if !=(FilePathOut,"")
#        savefig(FilePathOut)
#        println("Correlator vs r plotted to ", FilePathOut)
#    end
end

"""
Make a qualitative plot to put in the report: Gamma(r) vs r for different values
of J, showing the MI-SF crossing when from exponential Gamma becomes power-law
"""
function PlotCorrelationFunctionsMIvsSF(FilePathIn::String;
                                        FilePathOut1="",
                                        FilePathOut2="",
                                        PhaseBoundariesFilePath="",
                                        L=70,
                                        μ0=0.6,
                                        JMin = 0.1,
                                        JMax = 0.35)

    @warn "Mode under construction."

#    plot()
#
#    # Read the input data
#    data = readdlm(FilePathIn, ';', '\n'; comments=true)
#
#    # Extract unique J, L values
#    LL = unique(data[:, 1])
#    JJ = unique(data[:, 2])
#
#    # Select the couplings
#    JJ = JJ[ (JJ .>= JMin) .& (JJ .<= JMax) ]
#
#    # Mastruzzo to extract array of Γ
#    function parse_array(str)
#        return parse.(Float64, split(strip(str, ['[', ']', ' ']), ','))
#    end
#    Γall = [parse_array(row[4]) for row in eachrow(data)]
#
#    for (j, J) in enumerate(JJ)
#        # Filter data for the current J and L
#        filter = (data[:, 2] .== J) .& (data[:, 1] .== L)
#
#        # Extract the correlators {Γ(r,J,L)}_r for the selected J,L
#        Γeven = Γall[filter][1] # this is an array, Γ(r even)
#        r = range(start=2, step=2, length=length(Γeven)) # r even
#
#        # Check if Γeven has at least two elements
#        if length(Γeven) < 2
#            println("Warning: Not enough data points for L = $L. Skipping.")
#            continue
#        end
#
#        J_round = round(J, digits=3)
#
#        scatter!(r, Γeven,
#            xlabel=L"$r$",
#            ylabel=L"$\Gamma(r)$",
#            title=L"Correlation function ($L=%$L, \mu_0 = %$μ0$)",
#            label=L"J=%$J_round",
#            markersize=2,
#            xscale=:log10,
#            yscale=:log10,  
#            xticks=[1,10,100],
#            xlimits=(1,100),
#            legend=:bottomleft)
#    end
#
#    if FilePathOut1 != ""
#        savefig(FilePathOut1)
#        println("First plot saved on file.")
#    else
#        gui()
#    end
#
#    if PhaseBoundariesFilePath != ""
#        # NOTE: the code works if the μ0 = 0.0 boundaries data are used.
#        # gr()
#        plot()
#
#        BoundariesData = readdlm(PhaseBoundariesFilePath, ';', '\n'; comments=true)
#        LL = unique(BoundariesData[:, 1])
#    
#        L_PB_index = argmin(abs.(LL .- L)) # best approximation of L
#        L_PhaseBoundaries = LL[L_PB_index]
#    
#        if L_PhaseBoundaries !== L
#            println("L=$L is not among horizontal data. Using the closest available size, L=$L_PhaseBoundaries. (μ0=$μ0)")
#        end
#    
#        # Filter data corresponding to L_PhaseBoundaries
#        indices = (BoundariesData[:, 1] .== L_PhaseBoundaries)
#        JJ_PB = BoundariesData[indices, 2]
#        μUp = BoundariesData[indices, 4]
#        μDown = BoundariesData[indices, 5]
#        
#        plot!(JJ_PB, -μDown, 
#            label=L"$\mu_c^\pm$",
#            seriestype=:scatter,
#            markersize=1.5,
#            color="black",
#            xlabel=L"$J$",
#            ylabel=L"$\mu$",
#            title=L"Phase boundaries ($L=%$L, \mu_0=%$μ0$)"
#            )
#        plot!(JJ_PB, μUp, seriestype=:scatter,
#            label="",
#            markersize=1.5,
#            color="black",
#            )
#
#        vspan!([minimum(JJ), maximum(JJ)]; alpha=0.2, color=MyColors[4], label=nothing)
#        plot!([minimum(JJ), maximum(JJ)], [0.6, 0.6], color=MyColors[4], label=L"$(J,\mu_0)$")
#
#        if FilePathOut2 != ""
#            savefig(FilePathOut2)
#            println("Second plot saved on file.")
#        else
#            gui()
#        end
#    end
end

"""
Read the fit results from txt.
Plot the results of the fits, i.e., K_∞ vs J, for all the rMin, rMax available.
"""
function PlotFitResultsK(FilePathIn::String,
                         μ0::Float64,
                         rMin::Int64;
                         FilePathOut="")
    
    @warn "Mode under construction."

#    # Read the input data
#    FittedData = readdlm(FilePathIn, ',', '\n'; comments=true)
#
#    # Extract unique J, rMin values
#    JJ = unique(FittedData[:, 1])
#    rrMax = unique(FittedData[:, 3])
#    rrMin = unique(FittedData[:, 2])
#
#    println("Chosen rMin = $rMin")
#
#    # if length(rrMin) != 1
#    #     print("Error! More than one rMin. Which one should I put in the title?")
#    #     return
#    # end
#    
#    AllK_∞ = FittedData[:,4]
#    e_AllK_∞ = FittedData[:,4]
#
#    plot()
#
#    for rMax in rrMax
#        # Filter data corresponding to the current rMin
#        Filter = (FittedData[:, 3] .== rMax) .& (FittedData[:, 2] .== rMin)
#        K_∞ = AllK_∞[Filter]
#        e_K_∞ = e_AllK_∞[Filter]
#
#        scatter!(JJ, K_∞,
#                 xlabel=L"$J$",
#                 ylabel=L"$K_\infty$",
#                 title=L"Luttinger parameter vs $J$ ($\mu_0 = %$μ0$, $r_\mathrm{min}=%$rMin$)",
#                 label=L"$r_\mathrm{max}=%$(Int64(rMax))$",
#                 markersize=2,
#                 legend=:topright)
#
#    end
#
#    hline!([0.5], label=L"$K_c=1/2$", color="black", linestyle=:dash)
#
#    if FilePathOut != ""
#        savefig(FilePathOut)
#        println("\nPlot of K_∞ vs J plotted to ", FilePathOut)
#    else
#        gui()  
#    end
end

function PlotIllustrativeResultK(FilePathIn::String,
                                 μ0::Float64,
                                 rMin::Int64, 
                                 rMax::Int64;
                                 FilePathOut="")

    @warn "Mode under construction."

#    # Read the input data
#    FittedData = readdlm(FilePathIn, ',', '\n'; comments=true)
#
#    # Extract unique J values
#    JJ = unique(FittedData[:, 1])
#
#    # Extract K_∞ for chosen rMin and rMax
#    Filter = (FittedData[:, 3] .== rMax) .& (FittedData[:, 2] .== rMin)
#    K_∞ = FittedData[Filter,4]
#    e_K_∞ = FittedData[Filter,5]
#
#    scatter!(JJ, K_∞,
#    yerr=e_K_∞,
#    xlabel=L"$J$",
#    ylabel=L"$K_\infty$",
#    title=L"Luttinger parameter vs $J$ ($\mu_0 = %$μ0$, $%$rMin \le r \le %$rMax$)",
#    label=nothing,
#    markersize=2,
#    color="black",
#    legend=:topright)
#
#    if FilePathOut != ""
#        savefig(FilePathOut)
#        println("Illustrative plot saved on file.")
#    else
#        gui()
#    end
end

# ------------------------ Selection plot for check ----------------------------

function PlotSelection(FilePathIn::String,
					   Selections::Matrix{Float64})
	
    @warn "Mode under construction."

#	BoundariesData = readdlm(FilePathIn, ',', '\n'; comments=true)
#	JJ = BoundariesData[:,1]
#	ΔEp = BoundariesData[:,2]
#	ΔEm = BoundariesData[:,3]
#	
#	plot()
#    plot!(JJ,
#          [ΔEp, -ΔEm],
#          label=[L"\mu_c^+ \, (L \rightarrow \infty)" L"\mu_c^- \, (L \rightarrow \infty)"], 
#          xlabel=L"J", ylabel=L"$\mu$", 
#          title=L"Fitted phase boundaries",
#          alpha=1.0)
#     
#    rectangle(l, r, u, d) = Shape(l .+ [0,r-l,r-l,0], d .+ [0,0,u-d,u-d])
#    for i in 1:size(Selections,1)
#    	Left, Right, Up, Down = Selections[i,:] 
#    	plot!(rectangle(Left,Right,Up,Down),
#    		  label="Selection",
#    		  linewidth=0,
#    		  opacity=0.2)
#	end
#	gui()		# .pdf file saved on /tmp, erased on boot
end

# --------------------------- State properties plot ----------------------------

function PlotPopulations(DirPathOut::String,
						 P::Matrix{Float64},
						 ModelParameters::Vector{Float64},
						 SF::Bool,				# Redundant?
						 ConserveNumber::Bool)
	
    @warn "Mode under construction."

#	L, N, nmax = Int64.(ModelParameters[1:3])
#	J, μ = ModelParameters[4:5]
#	
#	plot(grid=false,
#		 xlabel=L"$i$ (site)",
#		 ylabel=L"$|n\rangle$ (state)",
#		 size = (440, 220),
#		 minorticks = false)
#
#	heatmap!(1:size(P,1), 0:(size(P,2)-1), P', clim=(0,1),
#			 label=L"$\mathrm{Tr}\lbrace | i;n \rangle \langle i;n | \rho \rbrace$",
#			 legend=true,
#			 color=:coolwarm)
#
#	# Draw by hand: grid
#	m, n = size(P')
#	vline!(0.5:(n+0.5), c=:black, alpha=0.2, linewidth=0.2, label=false)
#	hline!(-0.5:(m-0.5), c=:black, alpha=0.2, linewidth=0.2, label=false)
#
#	if SF
#		if ConserveNumber
#			title!(L"SF state population (fixed $N=%$N$, $J=%$J$, $\mu=%$μ$)")
#			savefig(DirPathOut * "/SF_J=$(J)_μ=$(μ)_fixedN_populations.pdf")
#		elseif !ConserveNumber
#			title!(L"SF state population (optimal $N$, $J=%$J$, $\mu=%$μ$)")
#			savefig(DirPathOut * "/SF_J=$(J)_μ=$(μ)_free_populations.pdf")
#		end
#	elseif !SF
#		if ConserveNumber
#			title!(L"MI state population (fixed $N=%$N$, $J=%$J$, $\mu=%$μ$)")
#			savefig(DirPathOut * "/MI_J=$(J)_μ=$(μ)_fixedN_populations.pdf")
#		elseif !ConserveNumber
#			title!(L"MI state population (optimal $N$, $J=%$J$, $\mu=%$μ$)")
#			savefig(DirPathOut * "/MI_J=$(J)_μ=$(μ)_free_populations.pdf")
#		end
#	end

end

function PlotBipartiteEntropy(DirPathOut::String,
						 	  S::Vector{Float64},
							  ModelParameters::Vector{Float64},
							  SF::Bool,
							  ConserveNumber::Bool)
	
    @warn "Mode under construction."

#	L, N, nmax = Int64.(ModelParameters[1:3])
#	J, μ = ModelParameters[4:5]
#				  
#	scatter([l for l in 1:length(S)], S,
#		 	 xlabel=L"$\ell$ (link index)",
#			 ylabel=L"$S$ (bipartite entropy)",
#			 legend=false)
#
#	if SF
#		if ConserveNumber
#			title!(L"SF bipartite entropy (fixed $N=%$N$, $J=%$J$, $\mu=%$μ$)")
#			savefig(DirPathOut * "/SF_J=$(J)_μ=$(μ)_fixedN_entropy.pdf")
#		elseif !ConserveNumber
#			title!(L"SF bipartite entropy (optimal $N$, $J=%$J$, $\mu=%$μ$)")
#			savefig(DirPathOut * "/SF_J=$(J)_μ=$(μ)_free_entropy.pdf")
#		end
#	elseif !SF
#		if ConserveNumber
#			title!(L"MI bipartite entropy (fixed $N=%$N$, $J=%$J$, $\mu=%$μ$)")
#			savefig(DirPathOut * "/MI_J=$(J)_μ=$(μ)_fixedN_entropy.pdf")
#		elseif !ConserveNumber
#			title!(L"MI bipartite entropy (optimal $N$, $J=%$J$, $\mu=%$μ$)")
#			savefig(DirPathOut * "/MI_J=$(J)_μ=$(μ)_free_entropy.pdf")
#		end
#	end

end

# Support function outside workflow for compared plots

function PlotBipartiteEntropyCompared(DirPathIn::String,
									  DirPathOut::String,
									  L::Int64,
									  N::Int64,
									  nmax::Int64,
									  SFPoint::Vector{Float64},
									  MIPoint::Vector{Float64})
				  
    @warn "Mode under construction."

#	SFFilePathIn = DirPathIn * "SF_J=$(SFPoint[1])_μ=$(SFPoint[2]).txt"
#	MIFilePathIn = DirPathIn * "MI_J=$(MIPoint[1])_μ=$(MIPoint[2]).txt"
#	
#	SFData = readdlm(SFFilePathIn, ';', '\n'; comments=true)
#	SFBool = SFData[:,1]
#	
#	# Mastruzzo to extract array of Γ
#    function ParseArray(str)
#        return parse.(Float64, split(strip(str, ['[', ']', ' ']), ','))
#    end
#	SFEntropy = [ParseArray(row[end]) for row in eachrow(SFData)]
#	
#	MIData = readdlm(MIFilePathIn, ';', '\n'; comments=true)
#	MIBool = MIData[:,1]
#	MIEntropy = [ParseArray(row[end]) for row in eachrow(MIData)]
#	
#	for (j,ConserveNumber) in enumerate(SFBool)
#		
#		scatter(legend=:right)
#		if ConserveNumber
#			scatter!([l for l in 1:length(SFEntropy[j])], SFEntropy[j],
#					 xlabel=L"$\ell$ (left partition length)",
#					 ylabel=L"$S$ (bipartite entropy)",
#					 markersize=2,
#					 label=L"$J=%$(SFPoint[1])$, $\mu=%$(SFPoint[2])$ (SF)")
#					
#			scatter!([l for l in 1:length(MIEntropy[j])], MIEntropy[j],
#			         markersize=2,
#					 label=L"$J=%$(MIPoint[1])$, $\mu=%$(MIPoint[2])$ (MI)")
#		
#			title!(L"SF bipartite entropy (fixed $N=%$N$, %$L sites)")
#			savefig(DirPathOut * "/fixedN_compared_entropy.pdf")
#		elseif !ConserveNumber
#			scatter!([l for l in 1:length(SFEntropy[2])], SFEntropy[2],
#					 xlabel=L"$\ell$ (left partition length)",
#					 ylabel=L"$S$ (bipartite entropy)",
#					 markersize=2,
#					 label=L"$J=%$(SFPoint[1])$, $\mu=%$(SFPoint[2])$ (SF)")
#					
#			scatter!([l for l in 1:length(MIEntropy[j])], MIEntropy[j],
#					 markersize=2,
#					 label=L"$J=%$(MIPoint[1])$, $\mu=%$(MIPoint[2])$ (MI)")
#		
#			title!(L"SF bipartite entropy (optimal $N$, %$L sites)")
#			savefig(DirPathOut * "/free_compared_entropy.pdf")
#		end
#	end
end
