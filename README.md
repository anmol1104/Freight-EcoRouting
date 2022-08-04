[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5140095.svg)](https://doi.org/10.5281/zenodo.5140095)

# Freight Eco-Routing
This repo contains pertinent files for the analyses performed in Pahwa, A., & Jaller, M. Can Truck Eco-Routing Bridge the Gap in the Transition to Zero-Emission?. Available at SSRN 4113893.

The analyses requires two tools,
## point-to-point routing tool 

To evaluate private impacts of eco-routing for the carrier hauling trucks

## multi-class traffic assignment tool 

Developed using paired alternative segements (mTAPAS) to evaluate network-wide effects of system-wide freight eco-routing

```julia
traffic_assignment(;network, assignment=:UE, tol=1e-5, maxiters=20, maxruntime=600, log=:on)
```

multi-class Traffic Assignment by Paired Alternative Segments (iTAPAS) algorithm

Returns output.csv file with arc flows and arc costs for each vehicle class

Returns report.csv file summarzing iteration-wise total flow, total cost, relative gap and run time
 
### Generalized link cost function: `cᵏᵢⱼ = fᵏ(vᵢⱼ)tᵢⱼ`
- `cᵏᵢⱼ` : generalized link cost for link 𝑖𝑗 , vehicle class 𝑚
- `tᵢⱼ`  : travel time on link 𝑖𝑗
- `vᵢⱼ`  : travel speed on link 𝑖𝑗
- `fᵏ(v)`: rate function (\$ per hr) for vehicle class k, `fᵏ(v) = ∑ₙ ηⁿvⁿ`

### Required properties of the generalized cost function
- Strictly positive
- Monotonically non-decreasing
- Continuously differentiable

### Arguments
- `network::String`             : network (availabe at: https://github.com/anmol1104/Freight-EcoRouting/src/TA/network)
- `assignment::Symbol=:UE`      : User Equilibrium (UE) or System Optimal (SO) assigment
- `tol::Float=1e-5`             : tolerance level for relative gap convergence
- `maxiters::Integer=20`        : maximum number of iterations
- `maxruntime::Integer=600`     : maximum wall clock run time (s)
- `log::Symbol=:on`             : shows results for every iteration if log is on

### DataFiles (available at: https://github.com/anmol1104/Freight-EcoRouting/src/TA/network)
- class   : Enlists coefficients of `fᵐ(v)` for each class
- network : Details the topology of the network
- demand  : Enlists OD pairs and corresponding demand for each class in passenger car equivalent (PCE)

### IO Units
- length  : miles
- time    : hour
- volume  : litre
- mass    : kg
- cost    : \$

The analysis and results from this work were presented at Transportation Research Board (TRB) Annual Meeting - 2022. 

To clone and run this project on your local machine refer to Julia documentation: http://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project 
