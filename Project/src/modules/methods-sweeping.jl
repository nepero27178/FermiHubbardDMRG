#!/usr/bin/julia

using Base.Threads
using LinearAlgebra
using DelimitedFiles  # For writedlm
using GeometricalPredicates

# ------------------------------ Horizontal sweep ------------------------------

@doc raw"""
function HorizontalSweep(
		UserSubMode::String,	# Boundaries/StateAnalysis/Full
		L::Int64, 
		VV::Array{Float64},
		μ0::Float64,
		DMRGParametersXY::Vector{Any},
		DirPathOut::String
	)
	
Returns: none (inline print and data saved on file).

This function runs an horizontal sweep (fixed `μ0`) to extract observables and 
correlator at increasing `V/t`. The string variable `UserSubMode` can take up
the values \"Boundaries\", \"StateAnalysis\" and \"Full\".
- \"Boundaries\": extract only the phase boundaries;
- \"StateAnalysis\": only perform state analysis;
- \"Full\": run both algorithms.
"""
function HorizontalSweep(
		UserSubMode::String,	# Boundaries/StateAnalysis/Full
		L::Int64, 
		VV::Array{Float64},
		μ0::Float64,
		DMRGParametersXY::Vector{Any},
		BoundariesFilePathOut,
		StateAnalysisFilePathOut
	)
	
	# DMRG metrhods dictionary for later selection
	DMRGMethod = Dict([
		("Boundaries", "Fast"),
		("StateAnalysis", "StateAnalyzer"),
		("Full", "StateAnalyzer"),
	])
    println("Performing horizontal sweep at μ=$μ0...")
    
    # Initialize matrices
    l = length(VV)
    
    for (j, V) in enumerate(VV)
    	printstyled("\rRunning DMRG for V=$(round(V, digits=3)), ",
            "μ=$μ0 (simulation $j/$l for L=$L)",
            color=:yellow)
        
        # First DMRG call: real gran-canonical ground-stae
        GroundState = RunDMRGAlgorithm(
        	[L, L/2, 1.0, V, μ0, 0.0],
			DMRGParametersXY,
			DMRGMethod[UserSubMode];	# Selected mode
			verbose=false,
			FixedN=false,				# Let the particles number vary
			FixedParity=false			# Let the particles number parity vary
		)
		
		if  UserSubMode=="StateAnalyzer" || UserSubMode=="Full" 
			
            #TODO Add errors to the chain observables
    
			E, e, _, n, Cnn, S = GroundState
			ρ = sum(n)/L
			δn2M = GetBlockVariance(n, Cnn)	# Imported from physics-definitions
			StateAnalysisFile = open(StateAnalysisFilePathOut, "a")
				write(StateAnalysisFile, "$L; $V; $E; $(δn2M); $S\n")
			close(StateAnalysisFile)
            print("\e[2K")
			
		elseif  UserSubMode=="Boundaries"
			E, _ = GroundState
		end
		
		if  UserSubMode=="Boundaries" ||  UserSubMode=="Full"
		
			# Additional DMRG calls, needed to extract phase boundaries and compressibility
			EF = zeros(3) 	# Col 1: L, col 2: L-1, col 3: L-2
			EAF = zeros(5)	# Col 1: L/2+2, col 2: L/2+1, col 3: L/2, col 4: L/2-1, col 5: L/2-2
		
			# Three DMRG calls (HF boundary)
			for (k,n) in enumerate([0,1,2])
				TmpE, _ = RunDMRGAlgorithm(
			    	[L, L-n, 1.0, V, μ0, 0.0],
					DMRGParametersXY,
					"Fast"; # UserMode
					verbose=false,
					FixedN=true,
					FixedParity=true
				)
				EF[k] = TmpE
			end
			
			# Five DMRG calls (HAF boundary)
			for (k,n) in enumerate([-2,-1,0,1,2])
				TmpE, _ = RunDMRGAlgorithm(
			    	[L, L/2+n, 1.0, V, μ0, 0.0],
					DMRGParametersXY,
					"Fast"; # UserMode
					verbose=false,
					FixedN=true,
					FixedParity=true
				)
				EAF[k] = TmpE
			end
		
			uΔm1 = EF[2] - EF[1] 	# Unitary-fililng, minus 1 particle
			uΔm2 = EF[3] - EF[1]	# Unitary-fililng, minus 2 particles
			
			# Half-Δ boundary
			hΔm2 = EAF[1] - EAF[3] # Half-filling, minus 2 particles
			hΔm1 = EAF[2] - EAF[3] # Half-filling, minus 1 particle
			hΔp1 = EAF[4] - EAF[3] # Half-filling, plus 1 particle
			hΔp2 = EAF[5] - EAF[3] # Half-filling, plus 2 particles
		
			BoundariesFile = open(BoundariesFilePathOut, "a")
				write(BoundariesFile, "$L; $V; $(hΔm1); $(hΔp1); $(hΔm2); $(hΔp2); $(uΔm1); $(uΔm2)\n")
			close(BoundariesFile)
            print("\e[2K")
		
		end
    end
    printstyled("Data for L=$L saved on file!\n", color=:green)
end

# ----------------------------- Rectangular sweep ------------------------------

@doc raw"""
function RectangularSweep(
        UserSubMode::String,
		L::Int64,
		N::Int64,
		VV::Vector{Float64},
		μμ::Vector{Float64},
		DMRGParameters::Vector{Any},
		FilePathOut::String;
		FilePathIn=""
	)
	
Returns: none (data saved on file).

The first variable `UserSubMode` is chosen between \"Density\" (only the charge
density is extracted) and \"Full\" (also charge compressibility and charge
stiffness are). `Full` mode requires (much) more computational resources. The
function performs a rectangular sweep in the `[V/t,μ/t]` space for a closed 
chain with `L` site and initialized to `N` particles states. Data are saved at
`FilePathOut` after DMRG simulations with parameters `DMRGParameters`.

TODO: Edit documentation for this function.
"""
function RectangularSweep(
        UserSubMode::String,	# Density/Full
		L::Int64,
		N::Int64,
		VV::Vector{Float64},
		μμ::Vector{Float64},
		DMRGParametersArray::Vector{Vector{Any}},
		FilePathOut::String;
		FilePathIn="",			# Boundaries filepath (if absent, ignore)
		double=false			# Read +/- 1 particle boundaries
	)
    
    global ε
    DataFile = open(FilePathOut,"a")
    ρCache = 1/2	# Cache half-filling
    NTotAvg = 0
    UseBoundaries = false
    
    if FilePathIn!==""
    
    	UseBoundaries = true
    
    	# DMRGParameters::Vector{Vector{Any}}
    	DMRGParametersXY = DMRGParametersArray[1]
    	DMRGParametersIF = DMRGParametersArray[2]
    	DMRGParametersIAF = DMRGParametersArray[3]
    	
		BoundariesData = readdlm(FilePathIn, ';', '\n'; comments=true)
		uμm = 0	# Initialize
		hμp = 0 # Initialize
		hμm = 0 # Initialize
		
		jj =  (BoundariesData[:, 1] .== L)
		HotizontalVV =  BoundariesData[jj,2]
		uΔm1 = BoundariesData[jj,7]
		uΔm2 = BoundariesData[jj,8]
		
		# Half-Δ boundary
		hΔm2 = BoundariesData[jj,5]
		hΔm1 = BoundariesData[jj,3]
		hΔp1 = BoundariesData[jj,4]
		hΔp2 = BoundariesData[jj,6]
		
		if double
			Up = -uΔm2/2
			Intermediate = hΔp2/2
			Down = -hΔm2/2
		elseif !double
			Up = -uΔm1
			Intermediate = hΔp1
			Down = -hΔm1
		end
		
		# Correct by +/- 0.1 to include borders (arbitrary)
		DownLeft = Point( HotizontalVV[1]-0.1, Rectangularμμ[1]-0.1 )
		UpLeft = Point( HotizontalVV[1]-0.1, Rectangularμμ[end]+0.1 )
		DownRight = Point( HotizontalVV[end]+0.1, Rectangularμμ[1]-0.1 )
		UpRight = Point( HotizontalVV[end]+0.1, Rectangularμμ[end]+0.1 )
		
		IFPoints = vcat(
			[UpLeft, DownLeft], 
			[Point( HotizontalVV[jj], Up[jj] ) for jj in 1:length(Up)])
		IFPolygon = Polygon(IFPoints...)


		# Note: from Beta testing it appears the IAF polygon creates some 
		# difficulty in binary recognition of points. This is not so important 
		# however since the phase is complicated and relatively small.
		IAFPoints = vcat(
			[Point( VV[jj], Intermediate[jj] ) for jj in 1:length(Intermediate)],
			[Point( VV[end-jj+1], Down[end-jj+1] ) for jj in 1:length(Down)]
		)
		IAFPolygon = Polygon(IAFPoints...)
	
	end

    for (j,V) in enumerate(VV), (m,μ) in enumerate(μμ)
    
    	if UseBoundaries
    	
			P = Point(V,μ)
			
			# Evaluate expected phase
			isIF = inpolygon(IFPolygon, P)
			isIAF = inpolygon(IAFPolygon, P)
			if isIF && isIAF
				@warn "This point appears to be both IF and IAF... using IAF parameters." P
				isIF = false
			end
			
			print("\e[2K")
			if isIF && !isIAF
				DMRGParameters = DMRGParametersIF
				@info "Phase: IF"
			elseif isIAF && !isIF
				DMRGParameters = DMRGParametersIAF
				@info "Phase: IAF"
			elseif !isIF && !isIAF
				DMRGParameters = DMRGParametersXY
				@info "Phase: XY"
			end
			
		elseif !UseBoundaries
			
			DMRGParameters = DMRGParametersArray[1]
			
		end
            
        ModelParameters = [L, N, 1.0, V, μ, 0.0]
        
        if j < length(VV) || m < length(μμ)
			printstyled("\e[2KRunning DMRG for L=$L, V=$(round(V, digits=3)), ", 
				"μ=$(round(μ,digits=3)) (simulation $m/$(length(μμ)) in μ, ",
				"$j/$(length(VV)) in V)\e[1F",
		        color=:yellow)
		elseif j == length(VV) && m == length(μμ)	# Remove final escape code
			printstyled("\e[2KRunning DMRG for L=$L, V=$(round(V, digits=3)), ", 
				"μ=$(round(μ,digits=3)) (simulation $m/$(length(μμ)) in μ, ",
				"$j/$(length(VV)) in V)",
		        color=:yellow)
		end

		E = 0 # Initalize energy to save variableJulia comes with

        if UserSubMode=="Density"

            E, psi = RunDMRGAlgorithm(
            	[L, L/2, 1.0, V, μ, 0.0],		# Central
			    DMRGParameters,
			    "Fast"; # UserMode
			    verbose=false,
			    FixedN=false,
			    FixedParity=false
		    )
            
            n = expect(psi, "n")
            ρ = sum(n)/L
            
            Cnn = GetDensityCorrelator(psi)
            δ = GetBlockVariance(
            	n,
            	Cnn;
            	kVector=[floor(Int64, L/4)]
            )
            
            #TODO Eventually, print entire δ vector
            write(DataFile,"$V; $μ; $E; $ρ; $(δ[1])\n")

        elseif UserSubMode=="Complementary"

		    k = 0 # Initalize compressibility to save variable
		    D = 0 # Initalize stiffness to save variable
            @sync begin
                # Create tasks for each DMRG call
                task1 = @spawn RunDMRGAlgorithm(
                	[L, L/2, 1.0, V, μ, 0.0],		# Central
				    DMRGParameters,
				    "Fast"; # UserMode
				    verbose=false,
				    FixedN=false,
				    FixedParity=false
			    )
			    task2 = @spawn RunDMRGAlgorithm(
                	[L, L/2, 1.0, V, μ, ε], 		# Rotate clockwise
				    DMRGParameters,
				    "Fast"; # UserMode
				    verbose=false,
				    FixedN=false,
				    FixedParity=false
			    )
			    task3 = @spawn RunDMRGAlgorithm(
                	[L, L/2, 1.0, V, μ, -ε], 		# Rotate counter-clockwise
				    DMRGParameters,
				    "Fast"; # UserMode
				    verbose=false,
				    FixedN=false,
				    FixedParity=false
			    )

                # Wait for all tasks to complete and collect results
                E, psi = fetch(task1)
                EClock, _ = fetch(task2)
                ECounterClock, _ = fetch(task3)
                
                D = pi * L * (EClock+ECounterClock-2*E)/(4ε^2)
                NTotAvg = GetTotalFermionNumber(psi)
            end
            
            ρ = sum(n)/L
            if m>1
			    k = (ρ - ρCache) / (μ - μμ[m-1])	# Compressibility
			    ρCache = ρ							# Next step
            elseif m==1
            	k = NaN								# Avoid segmentation fault
            end
            
            # Unitary and half projection
            sites = siteinds(psi)
            UP = GetUnitaryMIProjector(sites)
            HP = GetHalfMIProjector(sites)
            
			uP = inner(psi', UP, psi)
			hP = inner(psi', HP, psi)            
            
            write(DataFile,"$V; $μ; $E; $(uP); $(hP); $k; $D\n")
        
		elseif UserSubMode=="Full"

		    k = 0 # Initalize compressibility to save variable
		    D = 0 # Initalize stiffness to save variable
            @sync begin
                # Create tasks for each DMRG call
                task1 = @spawn RunDMRGAlgorithm(
                	[L, L/2, 1.0, V, μ, 0.0],		# Central
				    DMRGParameters,
				    "Fast"; # UserMode
				    verbose=false,
				    FixedN=false,
				    FixedParity=false
			    )
			    task2 = @spawn RunDMRGAlgorithm(
                	[L, L/2, 1.0, V, μ, ε], 		# Rotate clockwise
				    DMRGParameters,
				    "Fast"; # UserMode
				    verbose=false,
				    FixedN=false,
				    FixedParity=false
			    )
			    task3 = @spawn RunDMRGAlgorithm(
                	[L, L/2, 1.0, V, μ, -ε], 		# Rotate counter-clockwise
				    DMRGParameters,
				    "Fast"; # UserMode
				    verbose=false,
				    FixedN=false,
				    FixedParity=false
			    )

                # Wait for all tasks to complete and collect results
                E, psi = fetch(task1)
                EClock, _ = fetch(task2)
                ECounterClock, _ = fetch(task3)
                
                D = pi * L * (EClock+ECounterClock-2*E)/(4ε^2)
                n = expect(psi, "n")
            end
            
            ρ = sum(n)/L
            if m>1
			    k = (ρ - ρCache) / (μ - μμ[m-1])	# Compressibility
			    ρCache = ρ							# Next step
            elseif m==1
            	k = NaN								# Avoid segmentation fault
            end
            
            Cnn = GetDensityCorrelator(psi)
            δ = GetBlockVariance(
            	n,
            	Cnn;
            	kVector=[floor(Int64, L/4)]
            )
            
            # Unitary and half projection
            sites = siteinds(psi)
            UP = GetUnitaryMIProjector(sites)
            HP = GetHalfMIProjector(sites)
            
			uP = inner(psi', UP, psi)
			hP = inner(psi', HP, psi)            
            
            write(DataFile,"$V; $μ; $E; $ρ; $(δ[1]); $(uP); $(hP); $k; $D\n")
        
        end

	end

    close(DataFile)
end

function RectangularSweepTest(
        UserSubMode::String,	# Density/Full
		L::Int64,
		N::Int64,
		VV::Vector{Float64},
		μμ::Vector{Float64},
		DMRGParametersArray::Vector{Vector{Any}},
		FilePathOut::String;
		FilePathIn="",			# Boundaries filepath (if absent, ignore)
		double=false			# Read +/- 1 particle boundaries
	)
    
    global ε
    DataFile = open(FilePathOut,"a")
    ρCache = 1/2	# Cache half-filling
    NTotAvg = 0
    UseBoundaries = false
    
    if FilePathIn!==""
    
    	UseBoundaries = true
    
    	# DMRGParameters::Vector{Vector{Any}}
    	DMRGParametersXY = DMRGParametersArray[1]
    	DMRGParametersIF = DMRGParametersArray[2]
    	DMRGParametersIAF = DMRGParametersArray[3]
    	
		BoundariesData = readdlm(FilePathIn, ';', '\n'; comments=true)
		uμm = 0	# Initialize
		hμp = 0 # Initialize
		hμm = 0 # Initialize
		
		jj =  (BoundariesData[:, 1] .== L)
		HotizontalVV =  BoundariesData[jj,2]
		uΔm1 = BoundariesData[jj,7]
		uΔm2 = BoundariesData[jj,8]
		
		# Half-Δ boundary
		hΔm2 = BoundariesData[jj,5]
		hΔm1 = BoundariesData[jj,3]
		hΔp1 = BoundariesData[jj,4]
		hΔp2 = BoundariesData[jj,6]
		
		if double
			Up = -uΔm2/2
			Intermediate = hΔp2/2
			Down = -hΔm2/2
		elseif !double
			Up = -uΔm1
			Intermediate = hΔp1
			Down = -hΔm1
		end
		
		# Correct by +/- 0.1 to include borders (arbitrary)
		DownLeft = Point( HotizontalVV[1]-0.1, Rectangularμμ[1]-0.1 )
		UpLeft = Point( HotizontalVV[1]-0.1, Rectangularμμ[end]+0.1 )
		DownRight = Point( HotizontalVV[end]+0.1, Rectangularμμ[1]-0.1 )
		UpRight = Point( HotizontalVV[end]+0.1, Rectangularμμ[end]+0.1 )
		
		IFPoints = vcat(
			[UpLeft, DownLeft], 
			[Point( HotizontalVV[jj], Up[jj] ) for jj in 1:length(Up)])
		IFPolygon = Polygon(IFPoints...)


		# Note: from Beta testing it appears the IAF polygon creates some 
		# difficulty in binary recognition of points. This is not so important 
		# however since the phase is complicated and relatively small.
		IAFPoints = vcat(
			[Point( VV[jj], Intermediate[jj] ) for jj in 1:length(Intermediate)],
			[Point( VV[end-jj+1], Down[end-jj+1] ) for jj in 1:length(Down)]
		)
		IAFPolygon = Polygon(IAFPoints...)
	
	end

    for (j,V) in enumerate(VV), (m,μ) in enumerate(μμ)
    
    	if UseBoundaries
    	
			P = Point(V,μ)
			
			# Evaluate expected phase
			isIF = inpolygon(IFPolygon, P)
			isIAF = inpolygon(IAFPolygon, P)
			if isIF && isIAF
				@warn "This point appears to be both IF and IAF... using IAF parameters." P
				isIF = false
			end
			
			print("\e[2K")
			if isIF && !isIAF
				DMRGParameters = DMRGParametersIF
				@info "Phase: IF"
			elseif isIAF && !isIF
				DMRGParameters = DMRGParametersIAF
				@info "Phase: IAF"
			elseif !isIF && !isIAF
				DMRGParameters = DMRGParametersXY
				@info "Phase: XY"
			end
			
		elseif !UseBoundaries
			
			DMRGParameters = DMRGParametersArray[1]
			
		end
            
        ModelParameters = [L, N, 1.0, V, μ, 0.0]
        
        if j < length(VV) || m < length(μμ)
			printstyled("\e[2KRunning DMRG for L=$L, V=$(round(V, digits=3)), ", 
				"μ=$(round(μ,digits=3)) (simulation $m/$(length(μμ)) in μ, ",
				"$j/$(length(VV)) in V)\e[1F",
		        color=:yellow)
		elseif j == length(VV) && m == length(μμ)	# Remove final escape code
			printstyled("\e[2KRunning DMRG for L=$L, V=$(round(V, digits=3)), ", 
				"μ=$(round(μ,digits=3)) (simulation $m/$(length(μμ)) in μ, ",
				"$j/$(length(VV)) in V)",
		        color=:yellow)
		end

		E = 0 # Initalize energy to save variableJulia comes with
       
		E, psi  = RunDMRGAlgorithm(
			[L, L/2, 1.0, V, μ, 0.0],		# Central
			DMRGParameters,
			"Fast"; # UserMode
			verbose=false,
			FixedN=false,
			FixedParity=false
		)
                
        # Unitary and half projection
        sites = siteinds(psi)
        UP = GetUnitaryMIProjector(sites)
        HP = GetHalfMIProjector(sites)
        
		uP = inner(psi', UP, psi)
		hP = inner(psi', HP, psi)            
        
        write(DataFile,"$V; $μ; $E; $(uP); $(hP);\n")
    
    end

    close(DataFile)
end
