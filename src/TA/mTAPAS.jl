using DataFrames
using CSV
using Random
using Printf
using StatsBase
using Dates
Random.seed!(1104)

"""
    traffic_assignment(;network, assignment=:UE, tol=1e-5, maxiters=20, maxruntime=600, log=:on)

multi-class Traffic Assignment by Paired Alternative Segments (mTAPAS) algorithm
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
"""
function traffic_assignment(;network, assignment=:UE, tol=1e-5, maxiters=20, maxruntime=600, log=:on)
    println()
    printstyled("\nmTAPAS Algorithm", color=:blue)
    dir = joinpath(@__DIR__, "Network\\$network")

    # Algorithm parameters
    N   = Int64[]                               # Nodes
    A   = Array{Int64,1}[]                      # Arcs as adjacency list
    A′  = Array{Int64,1}[]                      # Arcs as adjacency list (head node as index and tail node in the vector)
    M   = Int64[]                               # Modes/classes
    Vᵢⱼ = Array{Float64,1}[]                    # Link volume capcity
    dᵢⱼ = Array{Float64,1}[]                    # Link length
    tᵢⱼ = Array{Float64,1}[]                    # Link free flow travel time
    αᵢⱼ = Array{Float64,1}[]                    # BPR parameters
    βᵢⱼ = Array{Float64,1}[]                    # BPR parameters
    R   = Int64[]                               # Origins
    Mᵣ  = Dict{Int64,Int64}()                   # Mode/class assosciated with every origin
    Sᵣ  = Dict{Int64,Array{Int64,1}}()          # Destinations for every origin
    qᵣ  = Dict{Tuple{Int64,Int64},Float64}()    # Demand between OD pairs
    η   = Array{Float64,1}[]                    # Coefficients for the polynomial fᵐ(v)

    ϕ   = Array{Int64,1}[]
    γ   = Array{Array{Int64,1},1}[]

    # Network build
    # Fetches betwork files and builds network related vectors
    function build()
        # class file
        clssFile = joinpath(dir,"class.csv")
        csv₁ = CSV.File(clssFile)
        df₁ = DataFrame(csv₁)
        for m ∈ 1:nrow(df₁)
            push!(M,m)
            push!(η, [df₁[m, c] for c ∈ 2:ncol(df₁)])
        end

        # network file
        ntwkFile = joinpath(dir,"network.csv")
        csv₂ = CSV.File(ntwkFile, types=[Int64, Int64, Float64, Float64, Float64, Float64, Float64])
        df₂ = DataFrame(csv₂)
        head = df₂[!, 1]::Array{Int64,1}
        tail = df₂[!, 2]::Array{Int64,1}
        linkcapacity = df₂[!, 3]::Array{Float64,1}
        linklength = df₂[!, 4]::Array{Float64,1}
        linkfft = df₂[!, 5]::Array{Float64,1}
        alpha = df₂[!, 6]::Array{Float64,1}
        beta = df₂[!, 7]::Array{Float64,1}
        n = max(maximum(head), maximum(tail))
        for i ∈ 1:n
            push!(N, i)
            push!(A, [])
            push!(A′,[])
            push!(Vᵢⱼ, [])
            push!(dᵢⱼ, [])
            push!(tᵢⱼ, [])
            push!(αᵢⱼ, [])
            push!(βᵢⱼ, [])
            push!(ϕ, [])
        end
        for i ∈ eachindex(head)
            push!(A[head[i]], tail[i])
            push!(A′[tail[i]], head[i])
            push!(Vᵢⱼ[head[i]], linkcapacity[i])
            push!(dᵢⱼ[head[i]], linklength[i])
            push!(tᵢⱼ[head[i]], linkfft[i])
            push!(αᵢⱼ[head[i]], alpha[i])
            push!(βᵢⱼ[head[i]], beta[i])
            push!(ϕ[head[i]], 1)
        end

        # geofencing file
        for m ∈ M push!(γ, [[0 for j ∈ A[i]] for i ∈ N]) end
        if "geofence.csv" ∈ readdir(dir)
            geofFile = joinpath(dir,"geofence.csv")
            csv₃ = CSV.File(geofFile)
            df₃ = DataFrame(csv₃)
            for m ∈ M
                for r ∈ 1:nrow(df₃)
                    i = df₃[r,1]::Int64
                    j = df₃[r,2]::Int64
                    k = findfirst(isequal(j), A[i])
                    γ[m][i][k] = df₃[r,m+2]
                end
            end
        end

        # demand file
        dmndFile = joinpath(dir,"demand.csv")
        csv₄ = CSV.File(dmndFile)
        df₄ = DataFrame(csv₄)
        origin = df₄[!, 1]::Array{Int64,1}
        destination = df₄[!, 2]::Array{Int64,1}
        flows = df₄[!, 3:ncol(df₄)]::DataFrame
        dict = Dict{Int64,Array{Int64,1}}(r => [r] for r ∈ unique(origin))
        # Making copies of origin nodes for every vehicle class
        for m ∈ eachindex(M)
            if isone(m) continue end
            for rₒ ∈ unique(origin)
                r = length(N) + 1
                push!(N, r)
                push!(A[rₒ], r), push!(A, [rₒ])
                push!(Vᵢⱼ[rₒ], 1.0), push!(Vᵢⱼ, [1.0])
                push!(dᵢⱼ[rₒ], 0.0), push!(dᵢⱼ, [0.0])
                push!(tᵢⱼ[rₒ], 0.001), push!(tᵢⱼ, [0.001])
                push!(αᵢⱼ[rₒ], 0.0), push!(αᵢⱼ, [0.0])
                push!(βᵢⱼ[rₒ], 0.0), push!(βᵢⱼ, [0.0])
                push!(ϕ[rₒ], 0), push!(ϕ, [0])
                push!(dict[rₒ], r)
                for k ∈ M push!(γ[k][rₒ], 0), push!(γ[k], [0]) end
            end
        end
        for i ∈ 1:nrow(df₄)
            rₒ = origin[i]
            for j ∈ 1:(ncol(df₄)-2)
                r, s, m = dict[rₒ][j], destination[i], j
                if r ∉ R Sᵣ[r] = [] end
                if r ∉ R push!(R, r) end
                push!(Sᵣ[r], s)
                Mᵣ[r] = m
                qᵣ[(r,s)] = flows[i,j]
            end
        end
    end

    # Arc cost
    # Returns cost of arc (i,j) for class m given arc flow x (k = A[i]⁻¹(j))
    function cᵢⱼ(i, k, m, x, method)
        #j = A[i][k]

        α = αᵢⱼ[i][k]
        β = βᵢⱼ[i][k]
        tₒ= tᵢⱼ[i][k]
        V = Vᵢⱼ[i][k]

        t = tₒ * (1 + α * (abs(x)/V) ^ β)
        d = dᵢⱼ[i][k]
        v = d/t
        if isinf(v) v = 1.0e6 end

        c = 0.0
        if isequal(method, :UE) 
            for k ∈ 0:(length(η[m])-1) c += η[m][k+1] * v^k * t end 
        end
        if isequal(method, :SO) 
            t′ = tₒ * α * β * (abs(x) ^ (β - 1))/(V ^ β)
            if iszero(β) t′ = 0.0 end
            if isinf(t′) t′ = 1.0e6 end
            for k ∈ 0:(length(η[m])-1) c += η[m][k+1] * v^k * (t + x * t′ * (1 - k)) end
        end
        c = c * (1 + γ[m][i][k])
        return c
    end

    # Arc cost derivative
    # Returns derivative of cost of arc (i,j) at arc flow x (k = A[i]⁻¹(j))
    function c′ᵢⱼ(i, k, m, x, method)
        #j = A[i][k]

        α = αᵢⱼ[i][k]
        β = βᵢⱼ[i][k]
        tₒ= tᵢⱼ[i][k]
        V = Vᵢⱼ[i][k]

        t = tₒ * (1 + α * (abs(x)/V) ^ β)
        d = dᵢⱼ[i][k]
        v = d/t
        if isinf(v) v = 1.0e6 end

        t′ = tₒ * α * β * (abs(x) ^ (β - 1))/(V ^ β)
        if iszero(β) t′ = 0.0 end
        if isinf(t′) t′ = 1.0e6 end

        c′ = 0.0
        if isequal(method, :UE) 
            for k ∈ 0:(length(η[m])-1) c′ += η[m][k+1] * v^k * (1 - k) * t′ end 
        end
        if isequal(method, :SO) 
            t′′ = tₒ * α * β * (β - 1) * (abs(x) ^ (β - 2))/(V ^ β)
            if iszero(β) || isone(β) t′′ = 0.0 end
            if isinf(t′′) t′′ = 1.0e6 end
            for k ∈ 0:(length(η[m])-1)
                c′ += η[m][k+1] * v^k * (1-k) * (2t′ + x*(t′′ - k*(t′^2)/t))
            end
        end
        c′ = c′ * (1 + γ[m][i][k])
        return c′
    end

    # Segment cost
    # Returns cost for segment e given arc flows xₐ and arc costs c
    function cₑ(e, cₐ)
        c = 0.0
        for (n,i) ∈ enumerate(e[1:end-1])
            j = e[n+1]
            k = findfirst(isequal(j), A[i])::Int64
            c += cₐ[i][k]
        end
        return c
    end

    # Segment flow
    # Returns flow on segment e given arc flows xₐ
    function fₑ(e, xₐ)
        f = zeros(length(e)-1)
        for (n,i) ∈ enumerate(e[1:end-1])
            j = e[n+1]
            k = findfirst(isequal(j), A[i])::Int64
            f[n] = xₐ[i][k]
        end
        return minimum(f)
    end

    # Djikstra's label setting algorithm
    # Returns predecessor label L for every node i for least cost path from node r given arc costs cₐ
    function djk(cₐ, r)
        L = [if isequal(i, r) r else -1 end for i ∈ N]       # Predecessor label
        C = [if isequal(i, r) 0.0 else Inf end for i ∈ N]    # Cost label
        X = copy(N)                                          # Set of open nodes
        i = r
        deleteat!(X, i)
        while !isempty(X)
            for (k,j) ∈ enumerate(A[i])
                c = C[i] + cₐ[i][k]
                if c < C[j] && j ∈ X L[j], C[j] = i, c end
            end
            index = argmin([C[i] for i ∈ X])
            i = X[index]
            deleteat!(X, index)
        end
        return L
    end

    # Tree
    # Returns tree rooted at r given predecessor label L
    function tree(L, r)
        T = Array{Int64,1}[[] for j ∈ N]
        for j ∈ N
            i = L[j]
            if i ≠ j && i ≠ -1 push!(T[i], j) end
        end
        return T
    end

    # Path
    # Returns path between node r and s using predecessor label L
    function path(L, r, s)
        p = Int64[]
        i = s
        push!(p, i)
        while i ≠ r
            i = Int(L[i])
            push!(p, i)
        end
        reverse!(p)
        return p
    end

    # multi-class Traffic Assignment with Paired Alterantive Segments (iTAPAS)
    # Returns excel file with arc flows and arc cost for each class, and a log of iterations
    function mTAPAS(ϵ, θ, writeout)
        report = Dict("TF" => Float64[], "TC" => Float64[], "RG" => Float64[], "WT" => Float64[])

        xʳₐ = Dict(r => [[0.0 for j ∈ A[i]] for i ∈ N] for r ∈ R)                                     # Stores origin-based arc flows
        xₐ  = [[sum([xʳₐ[r][i][k] for r ∈ R]) for k ∈ eachindex(A[i])] for i ∈ N]                      # Stores arc flows
        cₐ  = [[[cᵢⱼ(i, k, m, xₐ[i][k], assignment) for k ∈ eachindex(A[i])] for i ∈ N] for m ∈ M]     # Stores arc cost
        c′ₐ = [[[c′ᵢⱼ(i, k, m, xₐ[i][k], assignment) for k ∈ eachindex(A[i])] for i ∈ N] for m ∈ M]    # Stores derivative of arc cost
        πʳₐ = Dict(r => [[0.0 for j ∈ A[i]] for i ∈ N] for r ∈ R)                                     # Stores arc reduced cost
        P   = Tuple{Array{Int64,1},Array{Int64,1}}[]                                                  # Stores PAS
        Lᵣ  = Dict(r => [if isequal(i, r) r else -1 end for i ∈ N] for r ∈ R)                         # Stores origin-based least cost lables
        rₚ  = Int64[]                                                                                 # Stores origin for PAS p

        # Checks if arc a fails reduced cost optimal conditions for origin r
        function ispotential(a, r)
            i, j = a
            m = Mᵣ[r]
            k = findfirst(isequal(j), A[i])::Int64
            pᵣᵢ = path(Lᵣ[r], r, i)
            pᵣⱼ = path(Lᵣ[r], r, j)
            uʳᵢ = cₑ(pᵣᵢ, cₐ[m])
            uʳⱼ = cₑ(pᵣⱼ, cₐ[m])
            cʳᵢⱼ = cₐ[m][i][k]
            πʳᵢⱼ = uʳᵢ + cʳᵢⱼ - uʳⱼ
            if xʳₐ[r][i][k] > ϵ && πʳᵢⱼ > θ return (true, πʳᵢⱼ)
            else return (false, 0.0) end
        end

        # Checks if PAS p assosciated with origin rₒ can be eliminated
        function isbad(p, rₒ)
            e₁, e₂ = p
            m = Mᵣ[rₒ]
            c₁, c₂ = cₑ(e₁, cₐ[m]), cₑ(e₂, cₐ[m])
            f₁, f₂ = fₑ(e₁, xʳₐ[rₒ]), fₑ(e₂, xʳₐ[rₒ])
            if (f₁ < ϵ || f₂ < ϵ) && (c₁ ≠ c₂) return true
            else return false end
        end

        # Shifts flows from higher cost segment to lower cost segment of PAS p 
        # on its assosciated origin rₒ, given cost difference is greater than λ
        function shift(p, rₒ, λ)
            e₁, e₂ = p
            m = Mᵣ[rₒ]

            c₁, c₂ = cₑ(e₁, cₐ[m]), cₑ(e₂, cₐ[m])
            if abs(c₂ - c₁) < λ return end

            c′₁, c′₂ = cₑ(e₁, c′ₐ[m]), cₑ(e₂, c′ₐ[m])
            f₁, f₂ = fₑ(e₁, xʳₐ[rₒ]), fₑ(e₂, xʳₐ[rₒ])
            Δ = (c₂ - c₁)/(c′₁ + c′₂)
            if isnan(Δ) δ = 0.0
            elseif Δ ≥ 0 δ = min(Δ, f₂)
            else δ = max(Δ, -f₁) end

            for (n,i) ∈ enumerate(e₁[1:end-1])
                j = e₁[n+1]
                k = findfirst(isequal(j), A[i])::Int64
                xʳₐ[rₒ][i][k] += δ
                xₐ[i][k] += δ
                for m ∈ M
                    cₐ[m][i][k] = cᵢⱼ(i, k, m, xₐ[i][k], assignment)
                    c′ₐ[m][i][k] = c′ᵢⱼ(i, k, m, xₐ[i][k], assignment)
                end
            end
            for (n,i) ∈ enumerate(e₂[1:end-1])
                j = e₂[n+1]
                k = findfirst(isequal(j), A[i])::Int64
                xʳₐ[rₒ][i][k] -= δ
                xₐ[i][k] -= δ
                for m ∈ M
                    cₐ[m][i][k] = cᵢⱼ(i, k, m, xₐ[i][k], assignment)
                    c′ₐ[m][i][k] = c′ᵢⱼ(i, k, m, xₐ[i][k], assignment)
                end
            end
        end

        # PAS identification for arc a wrt origin r using Maximum Cost Search
        function MCS(a, r)
            depth, maxdepth = 1, 2
            flag = true

            i, j = a
            m = Mᵣ[r]
            e₁, e₂ = Int64[], Int64[]
            pᵣⱼ = path(Lᵣ[r], r, j)

            while flag
                # Intialize
                lₖ = [if k ∈ a 1 elseif k ∉ pᵣⱼ 0 else -1 end for k ∈ N]
                L  = [if isequal(k,j) i else -1 end for k ∈ N]

                # Iterate
                t, h = i, j
                while true
                    v = t
                    h = t

                    f = 0.0
                    for p ∈ A′[v]
                        k = findfirst(isequal(v), A[p])
                        x = xʳₐ[r][p][k]
                        c = cᵢⱼ(p, k, m, x, assignment)
                        #if x > ϵ && x > f f, t = x, p end
                        if x > ϵ && c > f f, t = c, p end
                    end
                    
                    L[v] = t
                    if isequal(lₖ[t], -1)   # PAS found
                        e₁ = path(Lᵣ[r], t, j)
                        e₂ = path(L, t, j)
                        shift((e₁, e₂), r, 0)
                        bool,_ = ispotential(a, r)
                        if !bool || isequal(depth, maxdepth) flag = false
                        else depth += 1 end
                        break
                    elseif isone(lₖ[t])     # Cycle found
                        if isequal(depth, maxdepth) flag = false
                        else
                            if isequal(v, t) pᵥₜ = Int64[]
                            else
                                pᵥₜ = path(L, v, t)
                                push!(pᵥₜ, v)
                                δ = fₑ(pᵥₜ, xʳₐ[r])
                            end
                            for (n,i) ∈ enumerate(pᵥₜ[1:end-1])
                                k = findfirst(isequal(pᵥₜ[n+1]), A[i])::Int64
                                xʳₐ[r][i][k] -= δ
                                xₐ[i][k] -= δ
                                for m ∈ M
                                    cₐ[m][i][k] = cᵢⱼ(i, k, m, xₐ[i][k], assignment)
                                    c′ₐ[m][i][k] = c′ᵢⱼ(i, k, m, xₐ[i][k], assignment)
                                end
                            end
                            depth += 1
                        end
                        break
                    else                    # Continue
                        lₖ[t] = 1
                    end
                end
            end
            p = (e₁, e₂)
            return p
        end

        if isequal(log, :on)
            print("\n iter: iteration,  RG:Relative Gap,  TF:Total Flow,  TC: Total Cost,  WT: Wall Time (s)")
            print("\n iter  | logRG      | TF          | TC          | WT (s) ")
            print("\n ------|------------|-------------|-------------|--------")
        end
        
        ## Step 0: Intialization - AON assignment
        T =  Dates.format(now(), "HH:MM:SS:sss")
        tₒ = parse.(Int64, [T[1:2], T[4:5], T[7:8], T[10:12]])
        wt = 0.0
        iter = 0
        for r ∈ R
            m = Mᵣ[r]
            Lᵣ[r] = djk(cₐ[m], r)
            for s ∈ Sᵣ[r]
                qᵣₛ = qᵣ[r,s]
                pᵣₛ = path(Lᵣ[r], r, s)
                for (n,i) ∈ enumerate(pᵣₛ[1:end-1])
                    j = pᵣₛ[n+1]
                    k = findfirst(isequal(j), A[i])
                    xʳₐ[r][i][k] += qᵣₛ
                    xₐ[i][k] += qᵣₛ
                    for m ∈ M
                        cₐ[m][i][k] = cᵢⱼ(i, k, m, xₐ[i][k], assignment)
                        c′ₐ[m][i][k] = c′ᵢⱼ(i, k, m, xₐ[i][k], assignment)
                    end
                end
            end
        end

        # Iterate
        while true
            # Run Time
            T =  Dates.format(now(), "HH:MM:SS:sss")
            tₙ = parse.(Int64, [T[1:2], T[4:5], T[7:8], T[10:12]])
            wt = sum((tₙ - tₒ) .* [3600, 60, 1, 1/1000])

            # Relative Gap
            for r ∈ R Lᵣ[r] = djk(cₐ[Mᵣ[r]], r) end
            num, den = 0.0, 0.0
            for r ∈ R for s ∈ Sᵣ[r] num += qᵣ[r,s] * cₑ(path(Lᵣ[r], r, s), cₐ[Mᵣ[r]]) end end
            for r ∈ R for i ∈ N for k ∈ eachindex(A[i]) den += xʳₐ[r][i][k] * cₐ[Mᵣ[r]][i][k] end end end
            rg = 1 - num/den

            # Total network flow and cost
            tf = sum(sum.(xₐ))
            tc = 0.0
            for r ∈ R for i ∈ N for k ∈ eachindex(A[i]) tc += xʳₐ[r][i][k] * cᵢⱼ(i, k, Mᵣ[r], xₐ[i][k], :UE) end end end

            # Miscellaneous
            push!(report["RG"], log10(abs(rg)))
            push!(report["TF"], tf)
            push!(report["TC"], tc)
            push!(report["WT"], wt)
            if isequal(log, :on)
                if iter < 10 @printf("\n #%.0f    | %.3E | %.5E | %.5E | %.3f  ", iter, log10(abs(rg)), tf, tc, wt)
                else @printf("\n #%.0f   | %.3E | %.5E | %.5E |%.3f ", iter, log10(abs(rg)), tf, tc, wt) end
            end
            
            # Convergence Test
            if log10(abs(rg)) ≤ log10(tol) || iter ≥ maxiters || wt ≥ maxruntime break end
            iter += 1

            ## Step 1
            for r ∈ R
                m = Mᵣ[r]
                Lᵣ[r] = djk(cₐ[m], r)
                Tᵣ = tree(Lᵣ[r], r)
                ## Step 1.1: Indentify potential arcs
                for i ∈ N
                    for (k,j) ∈ enumerate(A[i])
                        if j ∈ Tᵣ[i] continue end
                        a = (i, j)
                        bool, πʳᵢⱼ = false, 0.0
                        bool, πʳᵢⱼ = ispotential(a, r)
                        πʳₐ[r][i][k] = πʳᵢⱼ
                        ## Step 1.2: Flow shift for potential arc
                        if bool
                            p = MCS(a,r)
                            e₁, e₂ = p
                            rₒ = r
                            if !isempty(e₁) && !isempty(e₂) && p ∉ P push!(P, p), push!(rₚ, rₒ) end
                        end
                    end
                end
                # Step 1.3: Local shift
                for k ∈ sample(eachindex(P), length(P) ÷ 4) shift(P[k], rₚ[k], rg/1000) end
            end
            ## Step 2
            for _ ∈ 1:40
                for (k,p) ∈ enumerate(P)
                    if isbad(p, rₚ[k]) deleteat!(P, k), deleteat!(rₚ, k)
                    else shift(p, rₚ[k], rg/1000) end
                end
            end
        end

        # Writing out files
        if writeout
            df₁ = DataFrame(from = Int64[], to = Int64[])
            for m ∈ M df₁[!, Symbol("flow class $m")] = Float64[] end
            for m ∈ M df₁[!, Symbol("cost class $m")] = Float64[] end
            for i ∈ N
                for (k,j) ∈ enumerate(A[i])
                    if isone(ϕ[i][k])
                        push!(df₁[!, :from], i)
                        push!(df₁[!, :to], j)
                        for m ∈ M
                            xᵐ = 0.0
                            for r ∈ R if isequal(Mᵣ[r], m) xᵐ += xʳₐ[r][i][k] end end
                            push!(df₁[!, Symbol("flow class $m")], xᵐ)
                        end
                        for m ∈ M push!(df₁[!, Symbol("cost class $m")], cᵢⱼ(i, k, m, xₐ[i][k], :UE)) end
                    end
                end
            end
            df₂ = DataFrame(ITER = [i for i ∈ eachindex(report["TF"])],
                            TF = report["TF"], TC = report["TC"],
                            LOGRG = report["RG"], WT = report["WT"])
            CSV.write(joinpath(dir,"output-$assignment-MCS.csv"), df₁)
            CSV.write(joinpath(dir,"report-$assignment-MCS.csv"), df₂)
        end
        
        println("\n")
        println("   Total run time: $wt")
        println("   Total network flow: $(sum([xʳₐ[r][i][k] * ϕ[i][k] for r ∈ R for i ∈ N for k ∈ eachindex(A[i])]))")
        println("   Total network cost: $(sum([xʳₐ[r][i][k] * cᵢⱼ(i, k, Mᵣ[r], xₐ[i][k], :UE) * ϕ[i][k] for r ∈ R for i ∈ N for k ∈ eachindex(A[i])]))")
        return
    end

    build()
    mTAPAS(1e-12, 1e-16, true)
end
# ────────────────────────────────────────────────────────────────────────────────
traffic_assignment(network="SCAG", tol=1e-30, maxiters=60, maxruntime=1800, log=:on)