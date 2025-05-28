#!/usr/bin/julia

# ---------------------------------- Heatmap -----------------------------------

@doc raw"""
function PlotHeatmap(
        L::Int64,
        FilePathIn::String;
        PhaseBoundariesFilePath="",
        EFilePathOut="",
        kFilePathOut="",
        DFilePathOut=""
    )
    
Returns: none (plots in folder).    
    
Plot heatmaps of ground-state energy, compressibility and charge stiffness of
the system from data saved in `FilePathIn`.
"""
function PlotHeatmap(
        L::Int64,
        FilePathIn::String;
        PhaseBoundariesFilePath="",
        EFilePathOut="",		    # Energy heatmap
        kFilePathOut="",			# Charge stiffness
        DFilePathOut=""             # Compressibility heatmap
    )
    
    # Extract data coming from rectangular-sweep
    Data = readdlm(FilePathIn, ';', '\n'; comments=true)
	@info "Data" Data

    VV = Data[:,1]
    μμ = Data[:,2]
    EE = Data[:,3]
    kk = Data[:,4]
    DD = Data[:,5]

    NumV = length(unique(VV))
    Numμ = length(unique(μμ))
    
    Energies = zeros(Numμ, NumV)
    Compressibilities = zeros(Numμ, NumV)
    Stiffnesses = zeros(Numμ, NumV)

    for jj in 1:NumV
        Energies[:,jj] = EE[ Numμ*(jj-1)+1 : Numμ*jj ]
        Compressibilities[:,jj] = kk[ Numμ*(jj-1)+1 : Numμ*jj ]
        Stiffnesses[:,jj] = DD[ Numμ*(jj-1)+1 : Numμ*jj ]
    end

    if EFilePathOut != ""
        hE = heatmap(
        	unique(VV), unique(μμ), Energies, 
			xlabel=L"V/t",
			ylabel=L"$\mu/t$",
			title=L"Ground-state energy ($L=%$L$)"
		)

#        if PhaseBoundariesFilePath != ""
#            HeatmapAddPhaseBoundaries(PhaseBoundariesFilePath, L, JJ, μμ)
#        end

        savefig(hE, EFilePathOut)
        println("Energy plot for L=$L saved on file!")
    end

    if kFilePathOut != ""
        hk = heatmap(
        	unique(VV), unique(μμ)[3:end], Compressibilities[3:end,:], 
			xlabel=L"V/t",
			ylabel=L"$\mu/t$",
			title=L"Compressibility ($L=%$L$)"
		)

#        if PhaseBoundariesFilePath != ""
#            HeatmapAddPhaseBoundaries(PhaseBoundariesFilePath, L, JJ, μμ)
#        end

        savefig(hk, kFilePathOut)
        println("Compressibility plot for L=$L saved on file!")
    end
    
    if DFilePathOut != ""
        hD = heatmap(
        	unique(VV), unique(μμ), Stiffnesses, 
			xlabel=L"V/t",
			ylabel=L"$\mu/t$",
			title=L"Charge stiffness ($L=%$L$)"
		)

#        if PhaseBoundariesFilePath != ""
#            HeatmapAddPhaseBoundaries(PhaseBoundariesFilePath, L, JJ, μμ)
#        end

        savefig(hD, DFilePathOut)
        println("Charge stiffness plot for L=$L saved on file!")
    end
end

"""
Add the phase boundaries to the heatmap, for the closest size L available,
and adjusting the xlimits and ylimits appropriately.
"""
function HeatmapAddPhaseBoundaries(PhaseBoundariesFilePath::String,
                                   L::Int64,
                                   JJ::Array{Float64},
                                   μμ::Array{Float64};
                                   μ0=0.0)
	
	@warn "Function under construction..."	
	                                   
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
        MottLobeFilePath=PROJECT_ROOT * "/analysis/phase-boundaries/fitted-phase-boundaries.txt",
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
                    xlabel=L"$t$",      ylabel=L"$\Delta E_{\mathrm{gap}}$", 
                    title=L"Charge gap as a function of $t$ ($\mu_0=%$μ0$)",
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

@doc raw"""
function PlotPopulations(
		DirPathOut::String,
		P::Matrix{Float64},
		ModelParameters::Vector{Float64},
		XY::Bool,
		ConserveNumber::Bool
	)
	
Returns: none (plots saved).

This function plots, for a given parametrization `ModelParameters`, the state
populations stored in `P`.
"""
function PlotPopulations(
		DirPathOut::String,
		P::Matrix{Float64},
		ModelParameters::Vector{Float64},
		XY::Bool,
		ConserveNumber::Bool
	)
	
	# If I succeed in predicting the phase boundaries.

	L, N = Int64.(ModelParameters[1:2])
	t, V, μ, η = ModelParameters[3:6]
	
	plot(
		grid=false,
		xlabel=L"$i$ (site)",
		ylabel=L"$|n\rangle$ (state)",
		size = (440, 220),
		minorticks = false
	)

	heatmap!(
		1:size(P,1), 0:(size(P,2)-1), P', clim=(0,1),
		label=L"$\mathrm{Tr}\lbrace | i;n \rangle \langle i;n | \rho \rbrace$",
		legend=true,
		color=:coolwarm
	)

	# Draw by hand: grid
	m, n = size(P')
	vline!(0.5:(n+0.5), c=:black, alpha=0.2, linewidth=0.2, label=false)
	hline!(-0.5:(m-0.5), c=:black, alpha=0.2, linewidth=0.2, label=false)

	if XY
		if ConserveNumber
			title!(L"XY state population (fixed $N=%$N$, $V=%$V$, $\mu=%$μ$)")
			savefig(DirPathOut * "/XY_V=$(V)_μ=$(μ)_fixedN-populations.pdf")
		elseif !ConserveNumber
			title!(L"XY state population (optimal $N$, $V=%$V$, $\mu=%$μ$)")
			savefig(DirPathOut * "/XY_V=$(V)_μ=$(μ)_free-populations.pdf")
		end
	elseif !XY
		if ConserveNumber
			title!(L"IF state population (fixed $N=%$N$, $V=%$V$, $\mu=%$μ$)")
			savefig(DirPathOut * "/IF_V=$(V)_μ=$(μ)_fixedN-populations.pdf")
		elseif !ConserveNumber
			title!(L"IF state population (optimal $N$, $V=%$V$, $\mu=%$μ$)")
			savefig(DirPathOut * "/IF_V=$(V)_μ=$(μ)_free-populations.pdf")
		end
	end

end

@doc raw"""
function PlotBipartiteEntropy(
		DirPathOut::String,
		S::Vector{Float64},
		ModelParameters::Vector{Float64},
		XY::Bool,
		ConserveNumber::Bool
	)
	
Returns: none (plots saved).

This function plots, for a given parametrization `ModelParameters`, the state
entropy stored in `S`.
"""
function PlotBipartiteEntropy(
		DirPathOut::String,
		S::Vector{Float64},
		ModelParameters::Vector{Float64},
		XY::Bool,
		ConserveNumber::Bool
	)

	L, N = Int64.(ModelParameters[1:2])
	t, V, μ, η = ModelParameters[3:6]
				  
	scatter(
		[l for l in 1:length(S)], S,
		xlabel=L"$\ell$ (link index)",
		ylabel=L"$S$ (bipartite entropy)",
		legend=false
	)

	if XY
		if ConserveNumber
			title!(L"XY bipartite entropy (fixed $N=%$N$, $V=%$V$, $\mu=%$μ$)")
			savefig(DirPathOut * "/XY_V=$(V)_μ=$(μ)_fixedN-entropy.pdf")
		elseif !ConserveNumber
			title!(L"SF bipartite entropy (optimal $N$, $V=%$V$, $\mu=%$μ$)")
			savefig(DirPathOut * "/XY_V=$(V)_μ=$(μ)_free-entropy.pdf")
		end
	elseif !XY
		if ConserveNumber
			title!(L"IF bipartite entropy (fixed $N=%$N$, $V=%$V$, $\mu=%$μ$)")
			savefig(DirPathOut * "/IF_V=$(V)_μ=$(μ)_fixedN-entropy.pdf")
		elseif !ConserveNumber
			title!(L"IF bipartite entropy (optimal $N$, $V=%$V$, $\mu=%$μ$)")
			savefig(DirPathOut * "/IF_V=$(V)_μ=$(μ)_free-entropy.pdf")
		end
	end

end

@doc raw"""
function PlotBipartiteEntropyCompared(
		DirPathIn::String,
		DirPathOut::String,
		L::Int64,
		N::Int64,
		XYPoint::Vector{Float64},
		IFPoint::Vector{Float64}
	)

Returns: none (plots saved).

This function plots, for given parametrizations (points ,`XYPoint` and `IFPoint`
on the `(V,μ)` plane, the state entropy recovering it from the appropriate
FilePathIn. `PlotBipartiteEntropyCompared` is a support function outside the 
workflow to produce superimposed compared plots.
"""
function PlotBipartiteEntropyCompared(
		DirPathIn::String,
		DirPathOut::String,
		L::Int64,
		N::Int64,
		XYPoint::Vector{Float64},
		IFPoint::Vector{Float64}
	)

	XYFilePathIn = DirPathIn * "XY_V=$(XYPoint[1])_μ=$(XYPoint[2]).txt"
	IFFilePathIn = DirPathIn * "IF_V=$(IFPoint[1])_μ=$(IFPoint[2]).txt"
	
	XYData = readdlm(SFFilePathIn, ';', '\n'; comments=true)
	XYBool = SFData[:,1]
	
	# Mastruzzo to extract array of Γ
    function ParseArray(str)
        return parse.(Float64, split(strip(str, ['[', ']', ' ']), ','))
    end
	XYEntropy = [ParseArray(row[end]) for row in eachrow(XYData)]
	
	IFData = readdlm(IFFilePathIn, ';', '\n'; comments=true)
	IFBool = IFData[:,1]
	IFEntropy = [ParseArray(row[end]) for row in eachrow(IFData)]
	
	for (j,ConserveNumber) in enumerate(XYBool)
		
		scatter(legend=:right)
		if ConserveNumber
			scatter!([l for l in 1:length(XYEntropy[j])], XYEntropy[j],
					 xlabel=L"$\ell$ (left partition length)",
					 ylabel=L"$S$ (bipartite entropy)",
					 markersize=2,
					 label=L"$V=%$(XYPoint[1])$, $\mu=%$(XYPoint[2])$ (XY)")
					
			scatter!([l for l in 1:length(MIEntropy[j])], MIEntropy[j],
			         markersize=2,
					 label=L"$V=%$(IFPoint[1])$, $\mu=%$(IFPoint[2])$ (IF)")
		
			title!(L"XY bipartite entropy (fixed $N=%$N$, %$L sites)")
			savefig(DirPathOut * "/fixedN-compared-entropy.pdf")
		elseif !ConserveNumber
			scatter!([l for l in 1:length(XYEntropy[2])], XYEntropy[2],
					 xlabel=L"$\ell$ (left partition length)",
					 ylabel=L"$S$ (bipartite entropy)",
					 markersize=2,
					 label=L"$V=%$(XYPoint[1])$, $\mu=%$(XYPoint[2])$ (XY)")
					
			scatter!([l for l in 1:length(IFEntropy[j])], IFEntropy[j],
					 markersize=2,
					 label=L"$V=%$(IFPoint[1])$, $\mu=%$(IFPoint[2])$ (IF)")
		
			title!(L"XY bipartite entropy (optimal $N$, %$L sites)")
			savefig(DirPathOut * "/free-compared-entropy.pdf")
		end
	end
end
