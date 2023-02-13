using Random
using CSV
using DataFrames
using Statistics
using StatsBase
using Distributions
using Printf


"""
    ssp(origin, destination; network, parameter=["TT"], paradigm="expected value", distribution=Weibull(), threshold=1.0, leastcount=1/1000, numsims=100, showpath=false)

Stochastic Shortest Path
For a given paradigm, engine modes to operate in and parameters for the cost function, ssp performs numsims simulations for a vehicle traveling between origin-destination
and returns simulated travel statisitcs for travel distance, travel time, fuel consumed and emissions.

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

### IO Units
- distance  : miles
- energy    : litre of fuel
- emissions : kg
"""
function ssp(origin, destination; network, parameter=["TT"],
    paradigm="expected value", distribution=Weibull(), threshold=1.0,
    leastcount=1/1000, numsims=100, showpath=false)

    printstyled("\n----- $paradigm - $(join(parameter, ", "))-----\n", color=:cyan)
    dir = joinpath(@__DIR__, "Network\\$network")

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
    kₗ = Array{Int64,1}[]                   # Lower-bounds for pr → 0
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
        coefFile = joinpath(dir, "coef.csv")
        csv₁ = CSV.File(coefFile)
        df₁ = DataFrame(csv₁)
        for r ∈ 1:nrow(df₁)
            p = df₁[r,1]::String
            push!(parameters, p)
            for c ∈ 2:4
                append!(ℿ, df₁[r,5])
                append!(η, df₁[r,c])
                if p ∈ parameter append!(θ, 1.0) else append!(θ, 0.0) end
            end
        end

        # Network file
        ntwkFile = joinpath(dir, "network.csv")
        csv₂ = CSV.File(ntwkFile)
        df₂ = DataFrame(csv₂)
        head = df₂[!,1]::Array{Int64,1}
        tail = df₂[!,2]::Array{Int64,1}
        linkClass = df₂[!,3]::Array{Int64,1}
        linkLength = df₂[!,4]::Array{Float64,1}
        distShape = df₂[!,5]::Array{Float64,1}
        distScale = df₂[!,6]::Array{Float64,1}
        dist = typeof(distribution)
        n = max(maximum(head), maximum(tail))
        for i ∈ 1:n
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

        for i ∈ eachindex(head)
            append!(A[head[i]], tail[i])
            append!(M[head[i]], linkClass[i])
            append!(d[head[i]], linkLength[i])
            append!(α[head[i]], distShape[i])
            append!(β[head[i]], distScale[i])
        end

        # Geofence file
        append!(γ, [[0 for j ∈ A[i]] for i ∈ N])
        if "geofence.csv" ∈ readdir(dir)
            geofFile = joinpath(dir, "geofence.csv")
            csv₃ = CSV.File(geofFile)
            df₃ = DataFrame(csv₃)
            for r ∈ 1:nrow(df₃)
                i = df₃[r,1]::Int64
                j = df₃[r,2]::Int64
                k = findfirst(isequal(j), A[i])
                γ[i][k] =  df₃[r,3]
            end
        end

        # Generating cost metrics - mean, variance and proability distribution
        X(v) = [v^y for x ∈ eachindex(parameters) for y ∈ 0:2]
        L = Int(round(C̅/δ))
        Zₘ = Array{Float64,1}[[] for class ∈ sort(unique(linkClass))]
        # Generating random instances of link cost for a 1 mile long link of each link class
        for (m, class) ∈ enumerate(sort(unique(linkClass)))
            k = findfirst(isequal(class), linkClass)::Int64
            i = head[k]
            j = tail[k]
            k = findfirst(isequal(j), A[i])::Int64
            Random.seed!(i * k)
            V = rand(dist(α[i][k], β[i][k]), 1500)::Array{Float64,1}
            Z = [sum(η .* X(v) .* 1/v .* ℿ .* θ) for v ∈ V]
            Zₘ[m] = Z
        end
        # Using the random instances to generate link cost metrics for each link
        Threads.@threads for i ∈ N
             for (k,j) ∈ enumerate(A[i])
                 m = M[i][k]
                 C = Zₘ[m] .* d[i][k] .* (1 + γ[i][k]*1e6)
                 append!(μ[i], mean(C))
                 append!(σ²[i], var(C))
                 if isequal(paradigm, "reliability")
                    # Generating probability distribution
                     Lᶜ = Int(round(maximum(C) / δ))
                     h = fit(Histogram, C, (0.5:1:(max(L, Lᶜ) + 0.5)) * δ).weights
                     if iszero(sum(h)) push!(pr[i], h)
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
        𝐴 = Array{Int64,1}[[] for i ∈ N]
        𝑐ₐ= Array{Float64,1}[[] for i ∈ N]

        for i ∈ N
            for (k,j) ∈ enumerate(A[i])
                if isequal(type, "source")
                    push!(𝐴[i], A[i][k])
                    push!(𝑐ₐ[i], cₐ[i][k])
                else
                    push!(𝐴[j], i)
                    push!(𝑐ₐ[j], cₐ[i][k])
                end
            end
        end
 
        L = [if isequal(i, r) r else -1 end for i ∈ N]      # Predecessor label
        C = [if isequal(i, r) 0. else Inf end for i ∈ N]    # Cost label
        X = copy(N)                                         # Set of open nodes

        i = r
        deleteat!(X, i)
        while !isempty(X)
            for (k,j) ∈ enumerate(𝐴[i])
                c = C[i] + 𝑐ₐ[i][k]
                if c < C[j]  && j ∈ X L[j], C[j] = i, c end
            end
            index = argmin([C[j] for j ∈ X])
            i = X[index]
            deleteat!(X, index)
        end
        return L, C
    end

    # Djikstra's shortest path
    # Returns Djikstra's shortest path from/to node r to/from node s using label L
    function djkpath(L, r, s)
        source, sink = false, false
        if isequal(L[r], r) source = true end
        if isequal(L[s], s) sink = true end

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
        N′ = filter(x -> !(x ∈ prune), N)
        filter!(x -> (x ≠ s), N′)
        ρ = [[if isequal(i, s) 1.0 else 0.0 end for _ ∈ 0:L] for i ∈ N]    # Reliability
        λ = [[if isequal(i, s) s else -1 end for _ ∈ 0:L] for i ∈ N]       # Index of next node ∈ A[i]
        for l ∈ 1:L
            Threads.@threads for i ∈ N′
                ρ[i][l + 1], λ[i][l + 1] = findmax([sum([pr[i][k][kₒ] * ρ[j][l - kₒ + 1] for kₒ ∈ min(l, kₗ[i][k]):min(l, kᵤ[i][k])]) for (k,j) ∈ enumerate(A[i])])
            end
        end
        return ρ, λ
    end

    # Path reliability
    # Returns reliability of a given path (p)
    function reliability(p)
        L = Int(round(C̅/δ))
        ρ = [[if isequal(i, s) 1.0 else 0.0 end for _ ∈ 0:L] for i ∈ N]    # Reliability
        for l ∈ 1:L
            for m ∈ reverse(eachindex(p))
                if isone(m) continue end
                i, j = p[m-1], p[m]
                k = findfirst(isequal(j), A[i])::Int64
                ρ[i][l+1] = sum([pr[i][k][kₒ] * ρ[j][l - kₒ + 1] for kₒ ∈ min(l, kₗ[i][k]):min(l,kᵤ[i][k])])
            end
        end
        return ρ
    end

    # Monte-Carlo Simulation
    # Return simulated parameters (distances, travel time and emissions), and simulated paths
    function mcs()
        println("\nSimulating...")
        dist = typeof(distribution)
        L    = Int(round(C̅/δ))
        X(v) = [v^y for x ∈ eachindex(parameters) for y ∈ 0:2]

        paths = Array{Int64,1}[[] for _ ∈ 1:numsims]
        Z     = [[0.0 for _ ∈ 1:numsims] for _ ∈ eachindex(parameters)]
        C     = [0.0 for _ ∈ 1:numsims]

        if isequal(paradigm, "deterministic")
            for n ∈ 1:numsims
                Zₗ = [[[-1.0 for _ ∈ eachindex(A[i])] for i ∈ N] for _ ∈ eachindex(parameters)]
                Cₗ = [[-1.0 for _ ∈ eachindex(A[i])] for i ∈ N]
                for i ∈ N
                    for (k,j) ∈ enumerate(A[i])
                        Random.seed!(i * n)
                        v = rand(dist(α[i][k], β[i][k]))
                        z = η .* X(v) .* (d[i][k]/v)
                        for x ∈ eachindex(parameters) Zₗ[x][i][k] = sum(z[(3x-2):3x]) end
                        Cₗ[i][k] = sum(z .* ℿ .* θ)
                    end
                end
                Lᵣ, _ = djk(Cₗ, r, "source")
                p = djkpath(Lᵣ, r, s)
                for m ∈ eachindex(p)
                    if isone(m) continue end
                    i, j = p[m-1], p[m]
                    k = findfirst(isequal(j), A[i])::Int64
                    for x ∈ 1:length(parameters) Z[x][n] += Zₗ[x][i][j] end
                    C[n] += Cₗ[i][k]
                end
                paths[n] = p
            end

        elseif isequal(paradigm, "expected value")
            println("Algorithm run time")
            @time Lᵣ, _ = djk(μ, r, "source")
            @time p = djkpath(Lᵣ, r, s)
            #ρ = reliability(p)
            for n ∈ 1:numsims
                for m ∈ eachindex(p)
                    if isone(m) continue end
                    i, j = p[m-1],  p[m]
                    k = findfirst(isequal(j), A[i])::Int64
                    Random.seed!(i * n)
                    v = rand(dist(α[i][k], β[i][k]))
                    z = η .* X(v) .* (d[i][k]/v)
                    for x ∈ eachindex(parameters) Z[x][n] += sum(z[(3x-2):3x]) end
                    C[n] += sum(z .* ℿ .* θ)
                end
                paths[n] = p
            end

        elseif isequal(paradigm, "variance")
            println("Algorithm run time")
            @time Lᵣ, _ = djk(σ², r, "source")
            @time p = djkpath(Lᵣ, r, s)
            #ρ = reliability(p)
            for n ∈ 1:numsims
                for m ∈ eachindex(p)
                    if isone(m) continue end
                    i, j = p[m-1], p[m]
                    k = findfirst(isequal(j), A[i])::Int64
                    Random.seed!(i * n)
                    v = rand(dist(α[i][k], β[i][k]))
                    z = η .* X(v) .* (d[i][k]/v)
                    for x ∈ eachindex(parameters) Z[x][n] += sum(z[(3x-2):3x]) end
                    C[n] += sum(z .* ℿ .* θ)
                end
                paths[n] = p
            end

        elseif isequal(paradigm, "reliability")
            println("Pre-processing...")
            # Pruning
            cₘᵢₙ = [[kₗ[i][k] * δ for k ∈ eachindex(A[i])] for i ∈ N]
            _, Cᵣ = djk(cₘᵢₙ, r, "source")
            _, Cₛ = djk(cₘᵢₙ, s, "sink")
            for i ∈ N if Cᵣ[i] + Cₛ[i] > C̅ append!(prune, i) end end

            # Tie breaker
            _, Cₒ = djk(μ, s, "sink")
            for i ∈ N
                z = [if (A[i][k] ∈ prune) Inf else μ[i][k] + Cₒ[j] end for (k,j) ∈ enumerate(A[i])]
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
            for n ∈ 1:numsims
                i = r
                append!(paths[n], i)
                while i ≠ s
                    if (C̅ - C[n] ≥ δ && (prune ≠ N))     # Second condtion protects against very rare case when all nodes are pruned
                        l = Int(round((C̅ - C[n])/δ))
                        k = λ[i][l + 1]
                        Random.seed!(i * n)
                        v = rand(dist(α[i][k], β[i][k]))
                        z = η .* X(v) .* (d[i][k]/v)
                        for x ∈ eachindex(parameters) Z[x][n] += sum(z[(3x-2):3x]) end
                        C[n] += sum(z .* ℿ .* θ)
                        i = A[i][k]
                        append!(paths[n], i)
                    else
                        p = djkpath(Lₛ, i, s)             # Remaining route is completed as Least Expected Path
                        append!(paths[n], p[2:end])
                        for m ∈ eachindex(p)
                            if isone(m) continue end
                            i, j = p[m-1], p[m]
                            k = findfirst(isequal(j), A[i])::Int64
                            Random.seed!(i * n)
                            v = rand(dist(α[i][k], β[i][k]))
                            z = η .* X(v) .* (d[i][k]/v)
                            for x ∈ eachindex(parameters) Z[x][n] += sum(z[(3x-2):3x]) end
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
        pathCount = [length(findall(isequal(path), paths)) for path ∈ uniquePaths]
        df = DataFrame(stat = ["mean", "var", "median", "min", "max"])
        for i ∈ eachindex(parameters)
            df[!,Symbol(parameterize[parameters[i]])] = [@sprintf("%.2E",mean(Z[i])), @sprintf("%.2E", var(Z[i])),
            @sprintf("%.2E", median(Z[i])), @sprintf("%.2E", minimum(Z[i])), @sprintf("%.2E", maximum(Z[i]))]
        end
        df[!,Symbol("Cost")] = [@sprintf("%.2E", mean(C)), @sprintf("%.2E", var(C)), @sprintf("%.2E", median(C)), @sprintf("%.2E", minimum(C)), @sprintf("%.2E", maximum(C))]
        println("\nTravel Statistics:")
        println(df)

        #if paradigm ∈ ("expected value", "variance", "reliability") println("\nactual reliability: $(round(ρ[r][L + 1], digits=5))") else println("\nactual reliability: ", "NA") end
        println("calculated reliability: $(round(length(findall(x -> (x ≤ C̅), C))/numsims, digits=5))")
        if showpath
            println("\nPaths:")
            for path ∈ uniquePaths println("   Path => $path") end
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
