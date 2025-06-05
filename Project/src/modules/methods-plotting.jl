#!/usr/bin/julia

# ---------------------------------- Heatmap -----------------------------------

@doc raw"""
function PlotHeatmap(
        L::Int64,
        FilePathIn::String;
        PhaseBoundariesFilePath="",
        EFilePathOut="",
        ρFilePathOut="",
        δFilePathOut="",
        uPFilePathOut="",
    	hPFilePathOut="",
    	kFilePathOut="",
        DFilePathOut="",
        verbose=false
    )
    
Returns: none (plots in folder).    
    
Plot heatmaps of ground-state energy, compressibility, charge stiffness and
particle density of the system from data saved in `FilePathIn`.
"""
function PlotHeatmap(
        L::Int64,
        FilePathIn::String;
        PhaseBoundariesFilePath="",
        EFilePathOut="",		    # Energy heatmap
        ρFilePathOut="",			# Density heatmap
        δFilePathOut="",			# Block density variance heatmap
        uPFilePathOut="",			# Unitary filling MI projection heatmap
    	hPFilePathOut="",			# Half filling MI projection heatmap
    	kFilePathOut="",			# Charge stiffness
        DFilePathOut="",            # Compressibility heatmap
        verbose=false
    )
    
    # Extract data coming from rectangular-sweep
    Data = readdlm(FilePathIn, ';', '\n'; comments=true)
    if verbose
		@info "Data" Data
	end
	
	# Look for phase boundaries
	PhaseBoundariesDirPath = PROJECT_ROOT * "/../simulations/horizontal-sweep/boundaries/"
	F = filter(x -> occursin("$L",x), readdir(PhaseBoundariesDirPath))
	
	PlotBoundaries = false
	if length(F)==0
	
		@warn "No boundaries data found, skipping phase boundaries plot."
		PhaseBoundariesFilePathIn=""
	
	elseif length(F)==1
	
		PlotBoundaries = true
		@info "Boundaries data found, plotting phase boundaries."
		PhaseBoundariesFilePathIn = PhaseBoundariesDirPath * F[1]
	
	elseif length(F)>1
		
		PlotBoundaries = true
		@warn "Multiple boundaries data found, using the first one. Check the code methods-plotting.jl:64 to change this setting."
		jj = 1 # Change
		PhaseBoundariesFilePathIn = PhaseBoundariesDirPath * F[jj]
	
	end

    Full=false
    if kFilePathOut!=="" || DFilePathOut!==""
        Full=true
    end

    VV = Data[:,1]
    μμ = Data[:,2]
    EE = Data[:,3]
    ρρ = Data[:,4]
    δδ = Data[:,5]
    
    if Full
    	uPuP = Data[:,6]
    	hPhP = Data[:,7]
        kk = Data[:,8]
        DD = Data[:,9]
    end

    NumV = length(unique(VV))
    Numμ = length(unique(μμ))
    
    Energies = zeros(Numμ, NumV)
    Densities = zeros(Numμ, NumV)
    VarVar = zeros(Numμ, NumV)

    if Full
        Compressibilities = zeros(Numμ, NumV)
        Stiffnesses = zeros(Numμ, NumV)
        UnitaryProjections = zeros(Numμ, NumV)
        HalfProjections = zeros(Numμ, NumV)
    end

    for jj in 1:NumV
        Energies[:,jj] = EE[ Numμ*(jj-1)+1 : Numμ*jj ]
        Densities[:,jj] = ρρ[ Numμ*(jj-1)+1 : Numμ*jj ]
        VarVar[:,jj] = δδ[ Numμ*(jj-1)+1 : Numμ*jj ]
        
        if Full
            Compressibilities[:,jj] = kk[ Numμ*(jj-1)+1 : Numμ*jj ]
            Stiffnesses[:,jj] = DD[ Numμ*(jj-1)+1 : Numμ*jj ]
            UnitaryProjections[:,jj] = uPuP[ Numμ*(jj-1)+1 : Numμ*jj ]
            HalfProjections[:,jj] = hPhP[ Numμ*(jj-1)+1 : Numμ*jj ]
        end
    end

    if EFilePathOut != ""
    	
    	printstyled("Plotting energy...", color=:yellow)
    
        hE = heatmap(
        	unique(VV), unique(μμ), Energies, 
			xlabel=L"V/t",
			ylabel=L"$\mu/t$",
			label=L"E",
			title=L"Ground-state energy ($L=%$L$)"
		)

    	if PlotBoundaries                                      
			PlotPhaseBoundaries(
				PhaseBoundariesFilePathIn; 
				FilePathOut="",
				gap=false,
				double=true,
				CustomLL=[L],
				μ0=0.0,
				StandAlone=false
			)
		end

        savefig(hE, EFilePathOut)
        printstyled("\e[2K\e[1GEnergy plot for L=$L saved on file!\n", color=:green)
    end
    
    if ρFilePathOut != ""
    
    	printstyled("Plotting density...", color=:yellow)
    
        hρ = heatmap(
        	unique(VV), unique(μμ), Densities, 
			xlabel=L"V/t",
			ylabel=L"$\mu/t$",
			label=L"$\rho$",
			title=L"Charge density ($L=%$L$)"
		)

        if PlotBoundaries                                      
			PlotPhaseBoundaries(
				PhaseBoundariesFilePathIn; 
				FilePathOut="",
				gap=false,
				double=true,
				CustomLL=[L],
				μ0=0.0,
				StandAlone=false
			)
		end


        savefig(hρ, ρFilePathOut)
        printstyled("\e[2K\e[1GCharge density plot for L=$L saved on file!\n", color=:green)
    end
    
    if δFilePathOut != ""
    
    	printstyled("Plotting density block variance...", color=:yellow)
    
        hδ = heatmap(
        	unique(VV), unique(μμ), VarVar, 
			xlabel=L"V/t",
			ylabel=L"$\mu/t$",
			label=L"\delta_{L/4}^2",
			title=L"Charge density block variance ($L=%$L$, $M=\lfloor L/4 \rfloor$)"
		)

		if PlotBoundaries                                      
			PlotPhaseBoundaries(
				PhaseBoundariesFilePathIn; 
				FilePathOut="",
				gap=false,
				double=true,
				CustomLL=[L],
				μ0=0.0,
				StandAlone=false
			)
		end


        savefig(hδ, δFilePathOut)
        printstyled("\e[2K\e[1GCharge density block variance plot for L=$L saved on file!\n", color=:green)
    end

    
    if uPFilePathOut != ""
    
    	printstyled("Plotting unitary filling projection...", color=:yellow)
    
        huP = heatmap(
        	unique(VV), unique(μμ), UnitaryProjections, 
			xlabel=L"V/t",
			ylabel=L"$\mu/t$",
			label=L"$\langle \Psi | \hat P_{\mathrm{MI}_1} | \Psi \rangle$",
			title=L"$\mathrm{MI}_1$ projection ($L=%$L$)"
		)

		if PlotBoundaries                                      
			PlotPhaseBoundaries(
				PhaseBoundariesFilePathIn; 
				FilePathOut="",
				gap=false,
				double=true,
				CustomLL=[L],
				μ0=0.0,
				StandAlone=false
			)
		end


        savefig(huP, uPFilePathOut)
        printstyled("\e[2K\e[1GUnitary-filling projection plot for L=$L saved on file!\n", color=:green)
    end
    
    if hPFilePathOut != ""
    
    	printstyled("Plotting half filling projection...", color=:yellow)
    
        hhP = heatmap(
        	unique(VV), unique(μμ), HalfProjections, 
			xlabel=L"V/t",
			ylabel=L"$\mu/t$",
			label=L"$\langle \Psi | \hat P_{\mathrm{MI}_{1/2}} | \Psi \rangle$",
			title=L"$\mathrm{MI}_{1/2}$ projection ($L=%$L$)"
		)

		if PlotBoundaries                                      
			PlotPhaseBoundaries(
				PhaseBoundariesFilePathIn; 
				FilePathOut="",
				gap=false,
				double=true,
				CustomLL=[L],
				μ0=0.0,
				StandAlone=false
			)
		end

        savefig(hhP, hPFilePathOut)
        printstyled("\e[2K\e[1GHalf-filling projection plot for L=$L saved on file!\n", color=:green)
    end

    if kFilePathOut != ""
    
    	printstyled("Plotting charge compressibility...", color=:yellow)
    	
        hk = heatmap(
        	unique(VV), unique(μμ)[3:end], Compressibilities[3:end,:], 
			xlabel=L"V/t",
			ylabel=L"$\mu/t$",
			label=L"\kappa",
			title=L"Compressibility ($L=%$L$)",
			#clim=(0,1)
		)

		if PlotBoundaries                                      
			PlotPhaseBoundaries(
				PhaseBoundariesFilePathIn; 
				FilePathOut="",
				gap=false,
				double=true,
				CustomLL=[L],
				μ0=0.0,
				StandAlone=false
			)
		end


        savefig(hk, kFilePathOut)
        printstyled("\e[2K\e[1GCharge compressibility plot for L=$L saved on file!\n", color=:green)
    end
    
    if DFilePathOut != ""
    
    	printstyled("Plotting charge stiffness...", color=:yellow)
    	
        hD = heatmap(
        	unique(VV), unique(μμ), Stiffnesses, 
			xlabel=L"V/t",
			ylabel=L"$\mu/t$",
			label=L"\mathcal{D}",
			title=L"Charge stiffness ($L=%$L$)",
			clim=(-100,100)
		)

		if PlotBoundaries                                      
			PlotPhaseBoundaries(
				PhaseBoundariesFilePathIn; 
				FilePathOut="",
				gap=false,
				double=true,
				CustomLL=[L],
				μ0=0.0,
				StandAlone=false
			)
		end


        savefig(hD, DFilePathOut)
        printstyled("\e[2K\e[1GCharge stiffness plot for L=$L saved on file!\n", color=:green)
    end

end

# ----------------------------- Phase boundaries -------------------------------

@doc raw"""

"""
function PlotPhaseBoundaries(
        FilePathIn::String; 
        FilePathOut="",
        gap=false,
        double=false,
        CustomLL=[],
        μ0=0.0,
        StandAlone=true
    )
    
    BoundariesData = readdlm(FilePathIn, ';', '\n'; comments=true)

    # Extract unique L values
    if CustomLL==[]
        LL = unique(BoundariesData[:, 1])
    else
        LL = CustomLL
    end

	if StandAlone
		P = plot()
	end

	VV = 0	# Initialize VV to use it outside the for loop
	for (l, L) in enumerate(LL)
	    jj =  (BoundariesData[:, 1] .== L)
	    VV =   BoundariesData[jj,2]
	    uΔm1 = BoundariesData[jj,7]
		uΔm2 = BoundariesData[jj,8]
		
		# Half-Δ boundary
		hΔm2 = BoundariesData[jj,5]
		hΔm1 = BoundariesData[jj,3]
		hΔp1 = BoundariesData[jj,4]
		hΔp2 = BoundariesData[jj,6]
		
		if double
			uμm = uΔm2/2
			hμp = hΔp2/2
			hμm = hΔm2/2
		elseif !double
			uμm = uΔm1
			hμp = hΔp1
			hμm = hΔm1
		end
	    
	    if gap
	        plot!(
                VV, hμp - hμm, 
                xlabel=L"$V/t$", ylabel=L"$\Delta E_{\mathrm{gap}}$", 
                title=L"Charge gap as a function of $V/t$ ($\mu_0=%$μ0$)",
                seriestype=:scatter,
                markersize=1.5,
                label=L"$L=%$L$",
                color=MyColors[l % length(MyColors)]
            )
   	    else
   	    
   	    	if StandAlone
   	    		plot!(
		            xlabel=L"$V/t$", ylabel=L"$\mu$",
		            title=L"Extrapolation of $\mu_c^\pm$ ($\mu_0=%$μ0$, \texttt{double}=%$(double))",
		            legend=:outertopright,
		        )
   	    	end

            plot!(
            	xlims=[minimum(RectangularVV), maximum(RectangularVV)],
	            ylims=[minimum(Rectangularμμ), maximum(Rectangularμμ)],
				VV, -hμm .+ μ0,
				label="",
                linewidth=0.5,
                color=MyColors[l % length(MyColors)]
            )
            
			# Separate points and lines to correct legend            
            plot!(
				VV, -hμm .+ μ0,
				label=L"$L=%$(Int64(L))$",
				seriestype=:scatter,
                markershape=:circle,
                markersize=1.5,
                linewidth=0.5,
                color=MyColors[l % length(MyColors)]
            )

	        plot!(
                VV, hμp .+ μ0, 
                markershape=:circle,
                label="",
                markersize=1.5,
                linewidth=0.5,
                color=MyColors[l % length(MyColors)]
            )
            
            plot!(
                VV, -uμm .+ μ0, 
                markershape=:circle,
                linestyle=:dash,
                label="",
                markersize=1.5,
                linewidth=0.5,
                color=MyColors[l % length(MyColors)]
            )
	    end
	    
	end
	
	if !StandAlone
    	plot!(
    		legend=false
    	)
	end
    
    if !=(FilePathOut,"")
        savefig(FilePathOut)
        printstyled("Done!\n", color=:green)
    end
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

# -------------------------- State populations plot ----------------------------

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

# -------------------------------- Entropy plot --------------------------------

function ParseRawString(
            RawString::SubString{String}
        )
    return parse.(Float64, split(strip(RawString, ['[', ']', ' ']), ','))
end

@doc raw"""
function ChainPlots(
		DirPathIn::String,
		DirPathOut::String,
		L::Int64,
		XYPoint::Vector{Float64},
		IFPoint::Vector{Float64},
        IAFPoint::Vector{Float64}
	)

Returns: none (plots saved).

This function plots, for given parametrizations (points ,`XYPoint`, `IFPoint` 
and `IAFPoint` on the `(V,μ)` plane, the state entropy recovering it from the
appropriate `FilePathIn` built thanks to the function variables.
"""
function ChainPlots(
		DirPathIn::String,
		DirPathOut::String,
		L::Int64,
		XYPoint::Vector{Float64},
		IFPoint::Vector{Float64},
        IAFPoint::Vector{Float64}
	)

    PointDict = Dict([
        ("XY", XYPoint),
        ("IF", IFPoint),
        ("IAF", IAFPoint)
    ])

    LocalEnergy = Dict([])
    Density = Dict([])
    BlockDensityVariance = Dict([])
    Entropy = Dict([])
	
    for Phase in ["XY", "IF", "IAF"]

        Point = PointDict[Phase]
        FilePathIn = DirPathIn * Phase * "_V=$(Point[1])_μ=$(Point[2]).txt"
    	Data = readdlm(FilePathIn, ';', '\n'; comments=true)
        ConserveNumber = Data[:,1]
        ConserveParity = Data[:,2]
        
        Index = findall(.!ConserveNumber .&& .!ConserveParity)  # Both false
        
        # [1] is necessary to read the SubString
        LocalEnergy[Phase] = ParseRawString(Data[Index,8][1])
        Density[Phase] = ParseRawString(Data[Index,9][1])        
        BlockDensityVariance[Phase] = ParseRawString(Data[Index,10][1])
        Entropy[Phase] = ParseRawString(Data[Index,11][1])
    end

    MainDict = Dict([
        ("Local energy", LocalEnergy),
        ("Density", Density),
        ("Block density variance", BlockDensityVariance),
        ("Bipartite entropy", Entropy)
    ])
    MarkerStyleDict = Dict([
        ("XY", [MyColors[1], 1.5, :transparent, 0]),
        ("IF", [MyColors[4], 1.5, :transparent, 0]),
        ("IAF", [MyColors[3], 1.5, :transparent, 0])
#        ("XY", [MyColors[1], 2, :transparent, 0]),
#        ("IF", [:transparent, 2.5, MyColors[4], 0.5]),
#        ("IAF", [MyColors[3], 1.5, :transparent, 0])
    ])
    xLabelsDict = Dict([
        ("Local energy", L"$j$"),
        ("Density", L"$j$"),
        ("Block density variance", L"$M$"),
        ("Bipartite entropy", L"$\ell$")
    ])
    yLabelsDict = Dict([
        ("Local energy", L"$e_j$"),
        ("Density", L"$n_j$"),
        ("Block density variance", L"$\delta n_M^2$"),
        ("Bipartite entropy", L"$S_\ell$")
    ])

    for Observable in [k for k in keys(MainDict)] # Needed to get Vector{String}
        plot(
            title = Observable * " ($L sites)",
            legend=:best
        )
    
        for Phase in ["XY", "IF", "IAF"]

            Point = PointDict[Phase]

            yy = MainDict[Observable][Phase]
            xx = [l for l in 1:length(yy)]
		    plot!(
                xx, yy,
                xlabel = xLabelsDict[Observable],
                ylabel = yLabelsDict[Observable],
                markershape = :circle,                
                markercolor = MarkerStyleDict[Phase][1],
                markersize = MarkerStyleDict[Phase][2],
                markerstrokecolor = MarkerStyleDict[Phase][3],
                markerstrokewidth = MarkerStyleDict[Phase][4],
                linewidth = 0.5,
                label = L"$V=%$(Point[1])$, $\mu=%$(Point[2])$ (%$(Phase))",
            )

        end
    
		savefig(DirPathOut * "/" * lowercase(replace(Observable, ' '=>'-')) * ".pdf")
        printstyled("$(Observable) plot done!\n", color=:green)

    end

end
