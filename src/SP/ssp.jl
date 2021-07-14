using Random
using CSV
using DataFrames
using Statistics
using StatsBase
using Distributions
using Printf
cd(@__DIR__)


"""
    ssp(origin, destination[; networkName, parameter=["TT"], paradigm="expected value", distribution=Weibull(), threshold=1, leastcount=1/1000, η=100, showPath=false])

For a given paradigm, engine modes to operate in and parameters for the cost function, ssp performs η simulations 
for a vehicle traveling between origin-destination.

# Arguments
- `origin`::Integer                                 : Origin node
- `destination`::Integer                            : Destination node
- `networkName`::String                             : network
- `parameter`::Array{String}                        : distance (TD), time (TT), energy (FC), emissions (CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ)
- `paradigm`::String='expected value'               : determnistic, expected value, variance, reliability
- `threshold`::Float                                : threshold cost for reliability analysis
- `leastcount`::Float64                             : smallest value of discretized cost
- `distribution`::UnivariateDistribution            : link speed distribution function
- `η`::Integer                                      : Number of simulations
- `showPath`::Bool                                  : if true shows every path simulated

# IO Units
- distance  : miles
- energy    : litre of fuel
- emissions : kg
"""
function ssp(origin, destination; networkName, parameter=["TT"],
    paradigm="expected value", distribution=Weibull(), threshold=1,
    leastcount=1/1000, numSims=100, showPath=false)

    printstyled("\n----- $paradigm - $(join(parameter, ", "))-----\n", color=:cyan)

    # Algorithm parameters
    r, s = origin, destination
    C̅, δ = threshold, leastcount
    N = Int64[]                             # Nodes
    A = Array{Int64,1}[]                    # Arcs as adjacency list
    M = Array{Int64,1}[]                    # Link class
    d = Array{Float64,1}[]                  # Link length
    α = Array{Float64,1}[]                  # Link scale parameter
    β = Array{Float64,1}[]                  # Link shape parameter
    μ = Array{Float64,1}[]                  # Average link cost
    σ² = Array{Float64,1}[]                 # Var/iance in link cost
    pr = Array{Array{Float64,1},1}[]        # Link cost probabilities
    kₗ = Array{Int64,1}[]                    # Lower-bounds for pr → 0
    kᵤ = Array{Int64,1}[]                   # Upper-bounds for pr → 0
    prune = Int64[]                         # Pruned nodes
    parameters = String[]                   # Parameters
    ℿ = Float64[]                           # Cost parameters
    η = Float64[]                           # Coefficients
    γ = Array{Int64,1}[]                    # Binary geofence indicator
    θ = Float64[]                           # Binary internal array

    function popall!(arr) while !isempty(arr) pop!(arr) end end

    # Network build
    # Fetches netowrk files and builds network related vectors
    function build()
        println("\nBuilding network...")
        # Coefficient file
        coefFile = CSV.File("Network\\$networkName\\coef.csv")
        df₁ = DataFrame(coefFile)
        for r in 1:nrow(df₁)
            p = df₁[r,1]::String
            push!(parameters, p)
            for c in 2:4
                append!(ℿ, df₁[r,5])
                append!(η, df₁[r,c])
                if p in parameter append!(θ, 1.0) else append!(θ, 0.0) end
            end
        end

        # Network file
        networkFile = CSV.File("Network\\$networkName\\network.csv")
        df₂ = DataFrame(networkFile)
        head = df₂[!,1]::Array{Int64,1}
        tail = df₂[!,2]::Array{Int64,1}
        linkClass = df₂[!,3]::Array{Int64,1}
        linkLength = df₂[!,4]::Array{Float64,1}
        distShape = df₂[!,5]::Array{Float64,1}
        distScale = df₂[!,6]::Array{Float64,1}
        dist = typeof(distribution)
        n = max(maximum(head), maximum(tail))
        for i in 1:n
            append!(N, i)
            push!(A, [])
            push!(M, [])
            push!(d, [])
            push!(α, [])
            push!(β, [])
            push!(μ, [])
            push!(σ², [])
            push!(pr, [])
            push!(kₗ, [])
            push!(kᵤ, [])
        end

        for i in 1:length(head)
            append!(A[head[i]], tail[i])
            append!(M[head[i]], linkClass[i])
            append!(d[head[i]], linkLength[i])
            append!(α[head[i]], distShape[i])
            append!(β[head[i]], distScale[i])
        end

        # Geofence file
        append!(γ, [[0 for j in A[i]] for i in N])
        if "geofence.csv" in readdir("Network\\$networkName\\")
            geofFile = "Network\\$networkName\\geofence.csv"
            csv₃ = CSV.File(geofFile)
            df₃ = DataFrame(csv₃)
            for r in 1:nrow(df₃)
                i = df₃[r,1]::Int64
                j = df₃[r,2]::Int64
                k = findfirst(x -> (x == j), A[i])
                γ[i][k] =  df₃[r,3]
            end
        end

        # Generating cost metrics - mean, variance and proability distribution
        X(v) = [v^y for x in 1:length(parameters) for y in 0:2]
        L = Int(round(C̅/δ))
        Zₘ = Array{Float64,1}[[] for class in sort(unique(linkClass))]
        # Generating random instances of link cost for a 1 mile long link of each link class
        for (m, class) in enumerate(sort(unique(linkClass)))
            k = findfirst(x -> (x == class), linkClass)::Int64
            i, j = head[k], tail[k]
            k = findfirst(x -> (x == j), A[i])::Int64
            Random.seed!(i * k)
            V = rand(dist(α[i][k], β[i][k]), 1500)::Array{Float64,1}
            Z = [sum(η .* X(v) .* 1/v .* ℿ .* θ) for v in V]
            Zₘ[m] = Z
        end
        # Using the random instances to generate link cost metrics for each link
        Threads.@threads for i in N
             for (k,j) in enumerate(A[i])
                 m = M[i][k]
                 C = Zₘ[m] .* d[i][k] .* (1 + γ[i][k]*1e6)
                 append!(μ[i], mean(C))
                 append!(σ²[i], var(C))
                 if paradigm == "reliability" 
                    # Generating probability distribution
                     Lᶜ = Int(round(maximum(C) / δ))
                     h = fit(Histogram, C, (0.5:1:(max(L, Lᶜ) + 0.5)) * δ).weights
                     if sum(h) == 0 push!(pr[i], h)
                     else push!(pr[i], h / sum(h)) end
                     maxm, argmaxm = findmax(pr[i][k])
                     # Approximating pr → 0
                     ixₗ = findfirst(x -> (maxm/x  < 1000), pr[i][k][1:argmaxm])
                     ixᵤ = findfirst(x -> (maxm/x < 1000), reverse(pr[i][k][argmaxm:L]))
                     if (isnothing(ixₗ)) ixₗ = 1 end
                     if (isnothing(ixᵤ)) ixᵤ = 1  end
                     append!(kₗ[i], ixₗ)
                     append!(kᵤ[i], L - ixᵤ + 1)
                end
            end
        end

        println("Built fin.")
    end

    # Djikstra's label setting algorithm
    # Returns predecessor and cost labels L,C for every node i for least cost path from/to node r given arc costs cₐ
    function djk(cₐ, r, type)
        𝐴 = Array{Int64,1}[[] for i in N]
        𝑐ₐ= Array{Float64,1}[[] for i in N]

        for i in N
            for (k,j) in enumerate(A[i])
                if type == "source"
                    push!(𝐴[i], A[i][k])
                    push!(𝑐ₐ[i], cₐ[i][k])
                else
                    push!(𝐴[j], i)
                    push!(𝑐ₐ[j], cₐ[i][k])
                end
            end
        end
 
        L = [if i == r r else -1 end for i in N]       # Predecessor label
        C = [if i == r 0.0 else Inf end for i in N]    # Cost label
        X = copy(N)                                    # Set of open nodes

        i = r
        deleteat!(X, i)
        while !isempty(X)
            for (k,j) in enumerate(𝐴[i])
                c = C[i] + 𝑐ₐ[i][k]
                if c < C[j]  && j in X L[j], C[j] = i, c end
            end
            index = argmin([C[j] for j in X])
            i = X[index]
            deleteat!(X, index)
        end
        return L, C
    end

    # Djikstra's shortest path
    # Returns Djikstra's shortest path from/to node r to/from node s using label L
    function djkpath(L, r, s)
        source, sink = false, false
        if L[r] == r source = true end
        if L[s] == s sink = true end

        if sink r, s = s, r end

        p = [s]
        i = s
        while i ≠ r
            i = Int(L[i])
            append!(p, i)
        end

        if source reverse!(p) end
        return p
    end

    # Increasing Order of Time Budget
    # Returns reliability and successor node index matrix (Policy tables)
    function f()
        L = Int(round(C̅/δ))
        N′ = filter(x -> !(x in prune), N)
        filter!(x -> (x ≠ s), N′)
        ρ = [[if i==s 1.0 else 0.0 end for _ in 0:L] for i in N]    # Reliability
        λ = [[if i==s s else -1 end for _ in 0:L] for i in N]       # Index of next node in 𝒜[i]
        for l in 1:L
            Threads.@threads for i in N′
                ρ[i][l + 1], λ[i][l + 1] = findmax([sum([pr[i][k][kₒ] * ρ[j][l - kₒ + 1] for kₒ in min(l, kₗ[i][k]):min(l, kᵤ[i][k])]) for (k,j) in enumerate(A[i])])
            end
        end
        return ρ, λ
    end

    # Path reliability
    # Returns reliability of a given path (p)
    function reliability(p)
        L = Int(round(C̅/δ))
        ρ = [[if i==s 1.0 else 0.0 end for _ in 0:L] for i in N]    # Reliability
        for l in 1:L
            for m in length(p):-1:2
                i, j = p[m-1], p[m]
                k = findfirst(x -> (x == j), A[i])::Int64
                ρ[i][l+1] = sum([pr[i][k][kₒ] * ρ[j][l - kₒ + 1] for kₒ in min(l, kₗ[i][k]):min(l,kᵤ[i][k])])
            end
        end
        return ρ
    end

    # Monte-Carlo Simulation
    # Return simulated parameters (distances, travel time and emissions), and simulated paths
    function mcs()
        println("\nSimulating...")
        dist = typeof(distribution)
        L = Int(round(C̅/δ))
        X(v) = [v^y for x in 1:length(parameters) for y in 0:2]

        paths = Array{Int64,1}[[] for _ in 1:numSims]
        Z = [[0.0 for _ in 1:numSims] for _ in parameters]
        C = [0.0 for _ in 1:numSims]

        if paradigm == "deterministic"
            for n in 1:numSims
                Zₗ = [[[-1.0 for _ in 1:length(A[i])] for i in N] for _ in 1:length(parameters)]
                Cₗ = [[-1.0 for _ in 1:length(A[i])] for i in N]
                for i in N
                    for (k,j) in enumerate(A[i])
                        Random.seed!(i * n)
                        v = rand(dist(α[i][k], β[i][k]))
                        z = η .* X(v) .* (d[i][k]/v)
                        for x in 1:length(parameters) Zₗ[x][i][k] = sum(z[(3x-2):3x]) end
                        Cₗ[i][k] = sum(z .* ℿ .* θ)
                    end
                end
                Lᵣ, _ = djk(Cₗ, r, "source")
                p = djkpath(Lᵣ, r, s)
                for m in 2:length(p)
                    i, j = p[m-1], p[m]
                    k = findfirst(x -> (x == j), A[i])::Int64
                    for x in 1:length(parameters) Z[x][n] += Zₗ[x][i][j] end
                    C[n] += Cₗ[i][k]
                end
                paths[n] = p
            end

        elseif paradigm == "expected value"
            println("Algorithm run time")
            @time Lᵣ, _ = djk(μ, r, "source")
            @time p = djkpath(Lᵣ, r, s)
            #ρ = reliability(p)
            for n in 1:numSims
                for m in 2:length(p)
                    i, j = p[m-1],  p[m]
                    k = findfirst(x -> (x == j), A[i])::Int64
                    Random.seed!(i * n)
                    v = rand(dist(α[i][k], β[i][k]))
                    z = η .* X(v) .* (d[i][k]/v)
                    for x in 1:length(parameters) Z[x][n] += sum(z[(3x-2):3x]) end
                    C[n] += sum(z .* ℿ .* θ)
                end
                paths[n] = p
            end

        elseif paradigm == "variance"
            println("Algorithm run time")
            @time Lᵣ, _ = djk(σ², r, "source")
            @time p = djkpath(Lᵣ, r, s)
            #ρ = reliability(p)
            for n in 1:numSims
                for m in 2:length(p)
                    i, j = p[m-1], p[m]
                    k = findfirst(x -> (x == j), A[i])::Int64
                    Random.seed!(i * n)
                    v = rand(dist(α[i][k], β[i][k]))
                    z = η .* X(v) .* (d[i][k]/v)
                    for x in 1:length(parameters) Z[x][n] += sum(z[(3x-2):3x]) end
                    C[n] += sum(z .* ℿ .* θ)
                end
                paths[n] = p
            end

        elseif paradigm == "reliability"
            println("Pre-processing...")
            # Pruning
            cₘᵢₙ = [[kₗ[i][k] * δ for k in 1:length(A[i])] for i in N]
            _, Cᵣ = djk(cₘᵢₙ, r, "source")
            _, Cₛ = djk(cₘᵢₙ, s, "sink")
            for i in N if Cᵣ[i] + Cₛ[i] > C̅ append!(prune, i) end end

            # Tie breaker
            _, Cₒ = djk(μ, s, "sink")
            for i in N
                z = [if (A[i][k] in prune) Inf else μ[i][k] + Cₒ[j] end for (k,j) in enumerate(A[i])]
                k = sortperm(z)
                A[i] = A[i][k]
                d[i] = d[i][k]
                α[i] = α[i][k]
                β[i] = β[i][k]
                μ[i] = μ[i][k]
                σ²[i] = σ²[i][k]
                pr[i] = pr[i][k]
                kₗ[i] = kₗ[i][k]
                kᵤ[i] = kᵤ[i][k]
            end
            println("Pre-processing fin.")

            println("Algorithm run time")
            @time Lₛ, _ = djk(μ, s, "sink")
            @time ρ, λ = f()
            for n in 1:numSims
                i = r
                append!(paths[n], i)
                while i ≠ s
                    if (C̅ - C[n] ≥ δ && (prune ≠ N))     # Second condtion protects against very rare case when all nodes are pruned
                        l = Int(round((C̅ - C[n])/δ))
                        k = λ[i][l + 1]
                        Random.seed!(i * n)
                        v = rand(dist(α[i][k], β[i][k]))
                        z = η .* X(v) .* (d[i][k]/v)
                        for x in 1:length(parameters) Z[x][n] += sum(z[(3x-2):3x]) end
                        C[n] += sum(z .* ℿ .* θ)
                        i = A[i][k]
                        append!(paths[n], i)
                    else
                        p = djkpath(Lₛ, i, s)             # Remaining route is completed as Least Expected Path
                        append!(paths[n], p[2:end])
                        for m in 2:length(p)
                            i, j = p[m-1], p[m]
                            k = findfirst(x -> (x == j), A[i])::Int64
                            Random.seed!(i * n)
                            v = rand(dist(α[i][k], β[i][k]))
                            z = η .* X(v) .* (d[i][k]/v)
                            for x in 1:length(parameters) Z[x][n] += sum(z[(3x-2):3x]) end
                            C[n] += sum(z .* ℿ .* θ)
                        end
                        i = s
                    end
                end
            end
        end

        parameterize = Dict{String, String}("TD" => "TD", "TT" => "TT", "FC" => "FC",
        "CH4" => "CH₄", "CO"  => "CO", "CO2" => "CO₂" , "N2O" => "N₂O", "NOx" => "NOₓ", 
        "PM"  => "PM", "ROG" => "ROG", "SOx" => "SOₓ")
        uniquePaths = unique(paths)
        pathCount = [length(findall(x -> (x == path), paths)) for path in uniquePaths]
        df = DataFrame(stat = ["mean", "var", "median", "min", "max"])
        for i in 1:length(parameters)
            df[!,Symbol(parameterize[parameters[i]])] = [@sprintf("%.2E",mean(Z[i])), @sprintf("%.2E", var(Z[i])),
            @sprintf("%.2E", median(Z[i])), @sprintf("%.2E", minimum(Z[i])), @sprintf("%.2E", maximum(Z[i]))]
        end
        df[!,Symbol("Cost")] = [@sprintf("%.2E", mean(C)), @sprintf("%.2E", var(C)), @sprintf("%.2E", median(C)), @sprintf("%.2E", minimum(C)), @sprintf("%.2E", maximum(C))]
        println("\nTravel Statistics:")
        println(df)

        #if paradigm in ("expected value", "variance", "reliability") println("\nactual reliability: $(round(ρ[r][L + 1], digits=5))") else println("\nactual reliability: ", "NA") end
        println("calculated reliability: $(round(length(findall(x -> (x ≤ C̅), C))/numSims, digits=5))")
        if showPath
            println("\nPaths:")
            for path in uniquePaths println("   Path => $path") end
            println("Path count: $pathCount")
        end
        println("\nSimulation fin.")

        return Z, C, uniquePaths
    end

    printstyled("\nOrigin-Destination: $r - $s\n"; color=:red)
    @time build()
    @time mcs()
end
# ────────────────────────────────────────────────────────────────────────────────
