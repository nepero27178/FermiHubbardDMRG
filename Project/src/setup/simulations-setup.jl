#!/usr/bin/julia

Resolution = Dict([
	("Test", 3),			# 3*3 steps
	("Low", 10),			# 10*10 steps
	("Medium", 25),			# 25*25 steps
	("High", 50)			# 50*50 steps
])
ε = 0.01

# ------------------------------ State properties ------------------------------

StatePropertiesL = 30
IFPoint = [-2.5, 1.0]
IAFPoint = [5.0, 4.0]
XYPoint = [0.5, 1.0]

# ----------------  Ising Ferromagnet (analog) DMRG parameters -----------------

IFnSweeps = 5
IFMaxLinkDim = [10,20,30,40,50]
IFCutoff = [1E-8]
DMRGParametersIF = [IFnSweeps, IFMaxLinkDim, IFCutoff]

# -------------  Ising Anti-Ferromagnet (analog) DMRG parameters ---------------

IAFnSweeps = 5
IAFMaxLinkDim = [10,20,30,40,50]
IAFCutoff = [1E-8]
DMRGParametersIAF = [IAFnSweeps, IAFMaxLinkDim, IAFCutoff]

# ------------------------- XY (analog) DMRG paramters -------------------------

XYnSweeps = 20
XYMaxLinkDim = [10,50,75,200,500]
XYCutoff = [1E-8]
DMRGParametersXY = [XYnSweeps, XYMaxLinkDim, XYCutoff]
			 
# ----------------------------- Horizontal sweeps ------------------------------

HorizontalLL = [x for x in 14:8:38]
HorizontalVV = collect(range(start=-2.0, stop=4.0, length=Resolution["Medium"]))
Horizontalμμ = [0.0]

# ----------------------------- Rectangular sweeps -----------------------------

RectangularLL = [26]
RectangularVV = collect(range(start=-2.0, stop=4.0, length=Resolution["High"]))
Rectangularμμ = collect(range(start=0.0, stop=4.0, length=Resolution["High"]))
