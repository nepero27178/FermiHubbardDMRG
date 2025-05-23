#!/usr/bin/julia

LL = collect(2^l for l in 3:1:6)		# Sizes
NN = LL									# Unitary filling

# Ising Ferromagnet (analog) DMRG parameters
nsweeps = 10
maxlinkdim = [10,50,75,200,500]
cutoff = [1E-8]
DMRGParametersIF = [nsweeps, maxlinkdim, cutoff]

# XY (analog) DMRG paramters
nsweeps = 20
maxlinkdim = [10,50,75,200,500]
cutoff = [1E-8]
DMRGParametersXY = [nsweeps, maxlinkdim, cutoff]

Resolution = Dict([
				("Test", 3),			# 3*3 steps
			 	("Low", 10),			# 10*10 steps
			 	("Medium", 30),			# 30*30 steps
			 	("High", 50)			# 50*50 steps
			 ])
			 
# Horizontal sweeps
HorizontalLL = [10, 20, 30] 
HorizontalVV = collect(range(start=-3.0, stop=3.0, length=Resolution["Low"]))
Horizontalμμ = [0.6, 0.8]

# Rectangular sweeps

# RectangularLL = [10] # [10, 15, 20] 
# RectangularVV = collect(range(start=0.0, stop=0.35, length=Resolution["Medium"]))
# Rectangularμμ = collect(range(start=0.0, stop=1.0, length=Resolution["Medium"]))

RectangularLL = [16]
RectangularVV = collect(range(start=-3.0, stop=3.0, length=Resolution["Test"]))
Rectangularμμ = collect(range(start=0.0, stop=3.0, length=Resolution["Test"]))

# Rectangular selection sweeps
RectangularSelectionLL = [50] # [40, 50, 60]
# VV and μμ defined inside the script

ε = 0.01 # Arbitrary phase shift
