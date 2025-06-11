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
        KFilePathOut="",
        uFilePathOut="",
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
		#
        EFilePathOut="",		    # Energy heatmap
        ρFilePathOut="",			# Density heatmap
        δFilePathOut="",			# Block density variance heatmap
        #
        uPFilePathOut="",			# Unitary filling MI projection heatmap
    	hPFilePathOut="",			# Half filling MI projection heatmap
    	kFilePathOut="",			# Charge stiffness
        DFilePathOut="",            # Compressibility heatmap
        KFilePathOut="",			# K Luttinger heatmap
    	uFilePathOut="",				# u Luttinger heatmap
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

    Full=false	# Simpler conditioning
    if kFilePathOut!=="" || DFilePathOut!==""
        Full=true
    end
    
    # Import data

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
    	KK = sqrt.(abs.(pi .* kk .* DD))
    	uu = sqrt.(abs.(DD ./ (pi .* kk)))
    end
    
    KK[1:2,:] .= 0	# Avoid errors later
    uu[1:2,:] .= 0	# Avoid errors later

	# Filter data

    NumV = length(unique(VV))
    Numμ = length(unique(μμ))
    
    Energies = zeros(Numμ, NumV)
    Densities = zeros(Numμ, NumV)
    VarVar = zeros(Numμ, NumV)

    if Full
    
        UnitaryProjections = zeros(Numμ, NumV)
        HalfProjections = zeros(Numμ, NumV)    
        Compressibilities = zeros(Numμ, NumV)
        Stiffnesses = zeros(Numμ, NumV)    
	    Interactions = zeros(Numμ, NumV)
	    Velocities = zeros(Numμ, NumV)
    end

    for jj in 1:NumV
        Energies[:,jj] = EE[ Numμ*(jj-1)+1 : Numμ*jj ]
        Densities[:,jj] = ρρ[ Numμ*(jj-1)+1 : Numμ*jj ]
        VarVar[:,jj] = δδ[ Numμ*(jj-1)+1 : Numμ*jj ]
        
        if Full
        
        	UnitaryProjections[:,jj] = uPuP[ Numμ*(jj-1)+1 : Numμ*jj ]
            HalfProjections[:,jj] = hPhP[ Numμ*(jj-1)+1 : Numμ*jj ]
            Compressibilities[:,jj] = kk[ Numμ*(jj-1)+1 : Numμ*jj ]
            Stiffnesses[:,jj] = DD[ Numμ*(jj-1)+1 : Numμ*jj ]
	    	Interactions[:,jj] = KK[ Numμ*(jj-1)+1 : Numμ*jj ]
        	Velocities[:,jj] = uu[ Numμ*(jj-1)+1 : Numμ*jj ]
        	
        end
    end

    if EFilePathOut != ""
    	
    	printstyled("Plotting energy...", color=:yellow)
    
        hE = heatmap(
        	unique(VV), unique(μμ), Energies, 
			xlabel=L"V/t",
			ylabel=L"$\mu/t$",
			label=L"$E$",
			title=L"$E$ ($L=%$L$)"
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
			title=L"$\rho$ ($L=%$L$)",
			clim=(0,1)
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
			label=L"$\delta n_{L/4}^2$",
			title=L"$\delta n_{L/4}^2$ ($L=%$L$, $M=\lfloor L/4 \rfloor$)",
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
			title=L"$\langle \Psi | \hat P_{\mathrm{MI}_1} | \Psi \rangle$ ($L=%$L$)",
			clim=(0,1)
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
			title=L"$\langle \Psi | \hat P_{\mathrm{MI}_{1/2}} | \Psi \rangle$ ($L=%$L$)",
			clim=(0,0.1)
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
			title=L"$\kappa$ ($L=%$L$)",
			clim=(0,1)
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
			label=L"$\mathcal{D}$",
			title=L"$\mathcal{D}$ ($L=%$L$)",
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
    
    if KFilePathOut != "" && DFilePathOut != "" && kFilePathOut != ""
    
    	printstyled("Plotting K Luttinger parameter...", color=:yellow)
    	
        hK = heatmap(
        	unique(VV), unique(μμ)[3:end], Interactions[3:end,:], 
			xlabel=L"V/t",
			ylabel=L"$\mu/t$",
			label=L"$K$",
			title=L"$K$ ($L=%$L$)",
			clim=(0,2),
			#ylim=[Rectangularμμ[2], Rectangularμμ[end]]
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

        savefig(hK, KFilePathOut)
        printstyled("\e[2K\e[1GK Luttinger parameter plot for L=$L saved on file!\n", color=:green)
    end

    if uFilePathOut != "" && DFilePathOut != "" && kFilePathOut != ""
    
    	printstyled("Plotting u Luttinger parameter...", color=:yellow)
    	
        hu = heatmap(
        	unique(VV), unique(μμ)[3:end], Velocities[3:end,:], 
			xlabel=L"V/t",
			ylabel=L"$\mu/t$",
			label=L"$u$",
			title=L"$u$ ($L=%$L$)",
			clim=(0,100),
			#ylim=[Rectangularμμ[2], Rectangularμμ[end]]
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

        savefig(hu, uFilePathOut)
        printstyled("\e[2K\e[1Gu Luttinger parameter plot for L=$L saved on file!\n", color=:green)
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
		        UserColor = MyColors[l % length(MyColors)]
		    elseif !StandAlone
		    	UserColor = :yellow
   	    	end

            plot!(
            	xlims=[RectangularVV[1], RectangularVV[end]],
	            ylims=[Rectangularμμ[1], Rectangularμμ[end]],
				VV, -hμm .+ μ0,
				label="",
                linewidth=0.5,
                color=UserColor
            )
            
            if StandAlone
	            SizeLabel = L"$L=%$(Int64(L))$"
            elseif !StandAlone
            	SizeLabel = ""
            end
            
			# Separate points and lines to correct legend            
            plot!(
				VV, -hμm .+ μ0,
				label=SizeLabel,
				seriestype=:scatter,
                markershape=:circle,
                markersize=1.5,
                linewidth=0.5,
                color=UserColor
            )

	        plot!(
                VV, hμp .+ μ0, 
                markershape=:circle,
                label="",
                markersize=1.5,
                linewidth=0.5,
                color=UserColor
            )
            
            plot!(
                VV, -uμm .+ μ0, 
                markershape=:circle,
                linestyle=:dash,
                label="",
                markersize=1.5,
                linewidth=0.5,
                color=UserColor
            )
	    end
	    
	end
    
    if !=(FilePathOut,"")
        savefig(FilePathOut)
        printstyled("Done!\n", color=:green)
    end
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

# --------------------------- State-properties plots ---------------------------

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
		XYPointUp::Vector{Float64},
		XYPointDown::Vector{Float64},
		IFPoint::Vector{Float64},
        IAFPoint::Vector{Float64};
        CompactPlot=false
	)

Returns: none (plots saved).

This function plots, for given parametrizations (points ,`XYPoint`, `IFPoint` 
and `IAFPoint` on the `(V,μ)` plane, all the \"chain variables\": microscopic
energy, local density, bipartite entropy, block density variance and the CDW
and SU correlators. If `CompactPlot=true`, only the last four are plotted in a
grid scheme.
"""
function ChainPlots(
		DirPathIn::String,
		DirPathOut::String,
		L::Int64,
		XYPointUp::Vector{Float64},
		XYPointDown::Vector{Float64},
		IFPoint::Vector{Float64},
        IAFPoint::Vector{Float64};
        CompactPlot=false
	)

    PointDict = Dict([
        ("XY-Up", XYPointUp),
        ("XY-Down", XYPointDown),
        ("IF", IFPoint),
        ("IAF", IAFPoint)
    ])

    LocalEnergy = Dict([])
    Density = Dict([])
    BlockDensityVariance = Dict([])
    Entropy = Dict([])
    CDWCorrelator = Dict([])
	SUCorrelator = Dict([])
	
    for Phase in ["XY-Up", "XY-Down", "IF", "IAF"]

        Point = PointDict[Phase]
        FilePathIn = DirPathIn * Phase * "_V=$(Point[1])_μ=$(Point[2])_L=$(L).txt"
    	Data = readdlm(FilePathIn, ';', '\n'; comments=true)
        ConserveNumber = Data[:,1]
        ConserveParity = Data[:,2]
        
        Index = findall(.!ConserveNumber .&& .!ConserveParity)  # Both false
        
        # [1] is necessary to read the SubString
        LocalEnergy[Phase] = ParseRawString(Data[Index,8][1])
        Density[Phase] = ParseRawString(Data[Index,9][1])
        Entropy[Phase] = ParseRawString(Data[Index,10][1])        
        BlockDensityVariance[Phase] = ParseRawString(Data[Index,11][1])
        CDWCorrelator[Phase] = ParseRawString(Data[Index,12][1])
        SUCorrelator[Phase] = ParseRawString(Data[Index,13][1])
    end

    MainDict = Dict([
        ("Local energy", LocalEnergy),
        ("Density", Density),
        ("Block density variance", BlockDensityVariance),
        ("Bipartite entropy", Entropy),
        ("CDW Correlator", CDWCorrelator),
        ("SU Correlator", SUCorrelator)
    ])
    StyleDict = Dict([
      	# Name: color, size, border color, border size, line style, label
        ("XY-Up", [MyColors[1], 1.5, :transparent, 0, :solid, L"$\mathrm{XY}_1$"]),
        ("XY-Down", [MyColors[1], 1.5, :transparent, 0, :dash, L"$\mathrm{XY}_2$"]),
        ("IF", [MyColors[4], 1.5, :transparent, 0, :solid, L"$\mathrm{IF}$"]),
        ("IAF", [MyColors[3], 1.5, :transparent, 0, :solid, L"$\mathrm{IAF}$"])
    ])
    xLabelsDict = Dict([
        ("Local energy", L"$j$"),
        ("Density", L"$j$"),
        ("Block density variance", L"$M$"),
        ("Bipartite entropy", L"$\ell$"),
        ("CDW Correlator", L"$r$"),
        ("SU Correlator", L"$r$")
    ])
    yLabelsDict = Dict([
        ("Local energy", L"$e_j$"),
        ("Density", L"$n_j$"),
        ("Block density variance", L"$\delta n_M^2$"),
        ("Bipartite entropy", L"$S_\ell$"),
        ("CDW Correlator", L"$\mathcal{C}_\mathrm{CDW}(r)$"),
        ("SU Correlator", L"$\mathcal{C}_\mathrm{SU}(r)$")
    ])

	if !CompactPlot

		for (o,Observable) in enumerate([k for k in keys(MainDict)]) # Needed to get Vector{String}
		    plot(
		        title = Observable * " ($L sites)",
		        legend=:best
		    )
		    
		    if o > 6	# Correlators are the last two to be plotted
		    	plot!(
		    		xaxis=:log,
		    		yaxis=:log,
		    	)
		    end
		
		    for Phase in ["XY-Up", "XY-Down"]#, "IF", "IAF"]

		        Point = PointDict[Phase]

		        yy = MainDict[Observable][Phase]
		        xx = [l for l in 1:length(yy)]
				plot!(
		            xx, yy,
		            xlabel = xLabelsDict[Observable],
		            ylabel = yLabelsDict[Observable],
		            markershape = :circle,
		            markercolor = StyleDict[Phase][1],
		            markersize = StyleDict[Phase][2],
		            markerstrokecolor = StyleDict[Phase][3],
		            markerstrokewidth = StyleDict[Phase][4],
		            linewidth = 0.5,
		            linestyle = StyleDict[Phase][5],
		            linecolor= StyleDict[Phase][1],
		            label = L"$V=%$(Point[1])$, $\mu=%$(Point[2])$ (%$(StyleDict[Phase][6]))",
		        )

		    end
		
			savefig(DirPathOut * "/" * lowercase(replace(Observable, ' '=>'-')) * "_L=$L.pdf")
		    printstyled("$(Observable) plot done!\n", color=:green)

		end

	elseif CompactPlot
		
		OrderedObservables = [
			"Bipartite entropy",
			"Block density variance",
			"CDW Correlator",
			"SU Correlator"
		]
		PlotsList = Any[]
		
		for (o,Observable) in enumerate(OrderedObservables) # Needed to get Vector{String}
		    P = plot(
		    	size = (800,600),
		        title = Observable * " ($L sites)",
		    )
		
		    for Phase in ["XY-Up", "XY-Down", "IF", "IAF"]

		        Point = PointDict[Phase]

				if o==2
		        	LabelString = L"$V=%$(Point[1])$, $\mu=%$(Point[2])$ (%$(StyleDict[Phase][6]))"
		        else
					LabelString = ""
		        end
					

		        yy = MainDict[Observable][Phase]
		        xx = [l for l in 1:length(yy)]
				plot!(
		            xx, yy,
		            xlabel = xLabelsDict[Observable],
		            ylabel = yLabelsDict[Observable],
		            markershape = :circle,
		            markercolor = StyleDict[Phase][1],
		            markersize = StyleDict[Phase][2]+0.5,
		            markerstrokecolor = StyleDict[Phase][3],
		            markerstrokewidth = StyleDict[Phase][4],
		            linewidth = 0.5,
		            linestyle = StyleDict[Phase][5],
		            linecolor= StyleDict[Phase][1],
		            label = LabelString,
		            minorticks = false
		        )

		    end
			push!(PlotsList,P)
		end
		
		plot(
			PlotsList[1], PlotsList[2],
			PlotsList[3], PlotsList[4],
	        layout=4,
	        margin=3*Plots.mm,
	        legend=:outertopright
	    )
		
		savefig(DirPathOut * "/compact-plot_L=$L.pdf")
	    printstyled("Compact plot done!\n", color=:green)
		
	end

end
