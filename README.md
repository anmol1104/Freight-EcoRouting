[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5140095.svg)](https://doi.org/10.5281/zenodo.5140095)

# Freight Eco-Routing
This repo contains pertinent files for the analyses performed in Pahwa, A., & Jaller, M. Can Truck Eco-Routing Bridge the Gap in the Transition to Zero-Emission?. Available at SSRN 4113893.

The analyses requires two tools,
## point-to-point routing tool 

To evaluate private impacts of eco-routing for the carrier hauling trucks

```julia
ssp(origin, destination; network, parameter=["TT"], paradigm="expected value", distribution=Weibull(), threshold=1.0, leastcount=1/1000, numsims=100, showpath=false)
```

For a given paradigm, engine modes to operate in and parameters for the cost function, ssp performs numsims simulations for a vehicle traveling between origin-destination and returns simulated travel statisitcs for travel distance, travel time, fuel consumed, and emissions.

### Arguments
- `origin::Integer`                                 : origin node
- `destination::Integer`                            : destination node
- `network::String`                                 : network
- `parameter::Array{String}=["TT"]`                 : distance (TD), time (TT), energy (FC), emissions (CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ)
- `paradigm::String="expected value"`               : determnistic, expected value, variance, reliability
- `distribution::UnivariateDistribution=Weibull()`  : link speed distribution function
- `threshold::Float64=1.0`                          : threshold cost for reliability analysis
- `leastcount::Float64=1/1000`                      : smallest value of discretized cost
- `numsims::Integer=100`                            : number of simulations
- `showpath::Bool=false`                            : if true shows every path simulated

### DataFiles (available at: https://github.com/anmol1104/Freight-EcoRouting/tree/master/src/PPR/Network)
- coef    : Enlists parameter coefficients
- network : Details the topology of the network
- goefence: Enlists arcs in the geofence (optional)

### IO Units
- distance  : miles
- energy    : litre of fuel
- emissions : kg

## multi-class traffic assignment tool 

Developed using paired alternative segements (mTAPAS) to evaluate network-wide effects of system-wide freight eco-routing

```julia
traffic_assignment(;network, assignment=:UE, tol=1e-5, maxiters=20, maxruntime=600, log=:on)
```

multi-class Traffic Assignment by Paired Alternative Segments (iTAPAS) algorithm

Returns output.csv file with arc flows and arc costs for each vehicle class and report.csv file summarzing iteration-wise total flow, total cost, relative gap and run time
 
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
- `network::String`             : network (availabe at: https://github.com/anmol1104/Freight-EcoRouting/tree/master/src/TA/Network)
- `assignment::Symbol=:UE`      : User Equilibrium (UE) or System Optimal (SO) assigment
- `tol::Float=1e-5`             : tolerance level for relative gap convergence
- `maxiters::Integer=20`        : maximum number of iterations
- `maxruntime::Integer=600`     : maximum wall clock run time (s)
- `log::Symbol=:on`             : shows results for every iteration if log is on

### DataFiles (available at: https://github.com/anmol1104/Freight-EcoRouting/tree/master/src/TA/Network)
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
