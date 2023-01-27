using Plots
using JLD
using CSV
include("ssp.jl")
plotly()
Random.seed!(1104)


function main()
    # ────────────────────────────────────────────────────────────────────────────────
    network = "SCAG"
    loc = "From POLA"
    
    cv(x) = std(x)/mean(x)
    idr(x) = percentile(x, 90) - percentile(x, 10)


    # ───────────────────────────────────────────────────────────────────────────────
    P = ["TD", "TT", "FC", "CH4", "CO", "CO2", "N2O", "NOx", "PM", "ROG", "SOx"]                                # main parameters
    C = ["TD", "TT", "FC", "CH4", "CO", "CO2", "N2O", "NOx", "PM", "ROG", "SOx", "TC", "GHG", "CP"]             # Extended paramters/Criterion

    C₂P = Dict{String, Array{String,1}}("TD" => ["TD"], "TT" => ["TT"], "FC" => ["FC"],
        "CH₄" => ["CH4"], "CO"  => ["CO"], "CO₂" => ["CO2"], "N₂O" => ["N2O"],
        "NOₓ" => ["NOx"], "PM"  => ["PM"], "ROG" => ["ROG"], "SOₓ" => ["SOx"] ,
        "TC" => ["TD", "TT", "FC"], "GHG" => ["CH4", "CO2", "N2O", "ROG"], "CP" => ["CO", "NOx", "PM", "SOx"],
        "E" => ["CH4", "CO", "CO2", "NOx", "PM", "ROG"])  # Criteria to parameters
    
    C₂L = Dict{String, String}("TD" => "Travel Distance (TD)", "TT" => "Travel Time (TT)",
            "FC" => "Fuel Consumption (FC)", "CH₄" => "Methane emissions (CH₄)",
            "CO"  => "Carbon Monoxide emissions (CO)", "CO₂" => "Carbon Dioxide emissions (CO₂)",
            "N₂O" => "Nitrous Oxide emissions (N₂O)", "NOₓ" => "Nitrogen Oxide emissions (NOₓ)",
            "PM"  => "Particulate Matter emissions (PM)", "ROG" => "Reactive Organic Gases (ROG)",
            "SOₓ" => "Sulphur Oxide emissions (SOₓ)", "TC" => "Travel Cost (TC)")                               # Criteria to Label

    ℿ = Float64[]                           # Cost parameters
    ODs = Array{Int64,1}[]                  # Origin-Destination
    Q = Float64[]                           # OD demandusing

    # Coefficients file
    coefFile = joinpath(@__DIR__, "Network\\$network\\coef.csv")
    csv₁ = CSV.File(coefFile)
    df₁ = DataFrame(csv₁)
    for r ∈ 1:nrow(df₁) append!(ℿ, df₁[r,5]) end

    # Demand file
    dmndFile = joinpath(@__DIR__, "Network\\$network\\demand.csv")
    csv₂ = CSV.File(dmndFile)
    df₂ = DataFrame(csv₂)
    origins = df₂[!,1]
    destinations = parse.(Int64, String.(names(df₂))[2:end])
    for r ∈ origins for s ∈ destinations push!(ODs, [r, s]) end end
    for r ∈ 1:nrow(df₂) for c ∈ 2:ncol(df₂) append!(Q, df₂[r,c]) end end

    # Geofence file
    geofFile = joinpath(@__DIR__, "Network\\$network\\geofence - SELA.csv")
    csv₃ = CSV.File(geofFile)
    df₃ = DataFrame(csv₃)
    geofence = Dict{Tuple{Int64, Int64}, Int64}()
    for r ∈ 1:nrow(df₃) geofence[(df₃[r,1], df₃[r,2])] = df₃[r,3] end
    
    # ODPair [POLA, SB]: [4800, 4322]
    # LCP: RGB 64, 31, 127
    # FP : RGB 127, 64, 31
    # SP : RGB 127, 31, 94

    # Simulator ────────────────────────────────────────────────────────────────────────────────
    """
        sim(;criteria)

    Simulates path minimizing the expected value of the given criteria between all OD pairs

    # Arguments
    - `criteria::String`    : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP

    # Example
    - `sim(criteria = "TT")` returns simulated paramters for fastest paths between OD pairs
    """
    function sim(;criteria)
        X = [Array{Float64,1}[] for _ ∈ eachindex(ODs)]   # Simulated data from ssp
        for (k, OD) ∈ enumerate(ODs)
            r, s = OD
            Z, _, _ = ssp(r, s, network=network, parameter=C₂P[criteria], numsims=1000)
            append!(X[k], Z)
        end
        save(joinpath(@__DIR__, "Results\\$loc - $criteria.jld"), "X", X)
    end

    function arcs(;criteria)
        X = Dict{Tuple{Int64, Int64}, Array{Any, 1}}()
        for (k, OD) ∈ enumerate(ODs)
            r, s = OD
            _, _, uniquePaths = ssp(r, s, network=network, parameter=C₂P[criteria], numsims=1)
            for path ∈ uniquePaths
                i = r
                for l ∈ eachindex(path)
                    if isone(l) continue end
                    j = path[l]
                    if (i,j) ∉ keys(X) X[(i,j)] = [0, 0.0] end
                    X[(i,j)][1] += 1
                    X[(i,j)][2] += Q[k]
                    i = j
                end
            end
        end
        df = DataFrame(FROM = Int64[], TO = Int64[], GEOFENCE = Int64[], COUNT = Int64[], WTDCOUNT = Float64[])
        K = [k for k ∈ keys(X)]
        for (i,j) ∈ K
            push!(df[!, :FROM], i)
            push!(df[!, :TO], j)
            push!(df[!, :GEOFENCE], geofence[(i,j)])
            push!(df[!, :COUNT], X[(i,j)][1])
            push!(df[!, :WTDCOUNT], X[(i,j)][2])
        end
        CSV.write(joinpath(@__DIR__, "Results\\$loc - $criteria - arcs.csv"), df)
    end


    # Measure ────────────────────────────────────────────────────────────────────────────────
    """
        measure(parameters; criteria, metric, weighted=false)

    Returns metric on parameter from parameters for paths minimizing expected value for given criteria 
        
    # Arguments
    - `parameters::Array{String,1}`  : susbset of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria::String`             : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `metric::function`             : statistical function (mean, var, iqr, idr, cv etc.)
    - `weighted::Bool=false`         : if true results are weighted by demand between OD pairs

    # Example
    - `measure(["TT"], criteria="CO₂", metric=mean)` returns mean travel time path minimizing expected carbon dioxide
      emissions between all OD pairs.
    """
    function measure(parameters; criteria, metric, weighted=false)
        X = load(joinpath(@__DIR__, "Results\\$loc - $criteria.jld"))["X"]
        V = Float64[]
        for parameter ∈ parameters
            Y = [0.0 for _ ∈ eachindex(ODs)]
            w = [0.0 for _ ∈ eachindex(ODs)]
            for (k, _) ∈ enumerate(ODs)
                Z = [0.0 for _ ∈ 1:1000]
                for p ∈ C₂P[parameter]
                    i = findfirst(isequal(p), P)
                    Z += X[k][i] #* ℿ[i]
                end
                Y[k] = metric(Z)
                w[k] = weighted ? Q[k] : 1.0
            end
            println("$metric($parameter) on least $criteria path: ", (sum(Y, weights(w))))
            push!(V, mean(Y, weights(w)))
        end
        return V
    end



    # Δ plots ────────────────────────────────────────────────────────────────────────────────
    """
        Δ(parameter; criteria₁, criteria₂, metric, weighted=false)

    Returns %change in parameter's metric value for paths minimizing expected value for criteria₁ and criteria₂

    # Arguments
    - `parameter::String`       : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria₁::String`       : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria₂::String`       : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `metric::function`        : statistical function (mean, var, iqr, idr, cv etc.)
    - `weighted::Bool=false`    : if true results are weighted by demand between OD pairs
    
    # Example
    - `Δ("CO₂", criteria₁="TD", crieteria₂="TT", metric=mean)` returns % change in mean carbon dioxide emissions on 
      the shortest and fastest path between all OD pairs.
    """
    function Δ(parameter; criteria₁, criteria₂, metric, weighted=false, remove_zeros=false)
        X₁ = load(joinpath(@__DIR__, "Results\\$loc - $criteria₁.jld"))["X"]
        X₂ = load(joinpath(@__DIR__, "Results\\$loc - $criteria₂.jld"))["X"]

        Δ = [0.0 for _ ∈ ODs]
        w = [0.0 for _ ∈ ODs]
        for (k, _) ∈ enumerate(ODs)
            Z₁, Z₂ = X₁[k], X₂[k]

            P₁ = [0.0 for _ ∈ 1:1000]
            P₂ = [0.0 for _ ∈ 1:1000]
            for p ∈ C₂P[parameter]
                i = findfirst(x -> (x == p), P)
                P₁ += Z₁[i] * ℿ[i]
                P₂ += Z₂[i] * ℿ[i]
            end
            
            δ = (metric(P₂)/metric(P₁) - 1) * 100
            Δ[k] = δ
            w[k] = weighted ? Q[k] : 1.0
        end
        if remove_zeros
            index = findall(x -> x != 0, Δ)
            println("Re-routed trips: $(length(index)/length(Δ))")
            Δ = Δ[index]
            w = w[index] 
        end
        println("% change in $metric($parameter): ", mean(Δ, weights(w)))
        return mean(Δ, weights(w))
    end

    """
        distributeΔ(parameter; criteria₁, criteria₂, ODpair)
    
    Distribution of the parameter on path minimizing expected value for criteria₁ and criteria₂ between a given OD pair

    # Arguments
    - `parameter::String`    : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria₁::String`    : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria₂::String`    : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `ODpair::Array{Int,1}` : one of the OD pairs

    # Example
    - `distributeΔ("TT"; criteria₁="CO₂", criteria₂="NOₓ", ODpair=ODs[1])` returns distribution for travel time
      on path minimizing expected CO₂ emissions and path minimizing expected NOₓ emissions for the first OD pair
      in the list.

    """
    function distributeΔ(parameter; criteria₁, criteria₂, ODpair)
        k = findfirst(isequal(ODpair), ODs)
        X₁ = load(joinpath(@__DIR__, "Results\\$loc - $criteria₁.jld"))["X"]
        X₂ = load(joinpath(@__DIR__, "Results\\$loc - $criteria₂.jld"))["X"]

        Z₁, Z₂ = X₁[k], X₂[k]
        P₁ = [0.0 for _ ∈ 1:1000]
        P₂ = [0.0 for _ ∈ 1:1000]
        for p ∈ C₂P[parameter]
            i = findfirst(isequal(p), P)
            P₁ += Z₁[i] * ℿ[i]
            P₂ += Z₂[i] * ℿ[i]
        end

        fig = plot(xlab="$(C₂L[parameter]) (\$)", ylab="density",
            tickfont=(14,"calibri"), guidefont=(16,"calibri"), legendfont=(14, "calibri"),
            xlim=(0.95*minimum(minimum.([P₁, P₂])), 1.05*maximum(maximum.([P₁, P₂]))))
        fig = histogram!(P₁, color=RGBA(31/255,127/255,64/255,100/255),
            normalize=:pdf, linewidth=0.1, linecolor=:white, label="")
        fig = histogram!(P₂, color=RGBA(127/255,31/255,94/255,100/255),
            normalize=:pdf, linewidth=0.1, linecolor=:white, label="")
        fig = vline!([mean(P₁)], color=RGBA(31/255,127/255,64/255), linewidth=2, label="")
        fig = vline!([mean(P₂)], color=RGBA(127/255,31/255,94/255), linewidth=2, label="")
        display(fig)
    end

    """
        sprayΔ(parameter₁, parameter₂; criteria₁, criteria₂, metric, weighted=false)

    Scatter plot of % Δmetric(parameter₂) vs % Δmetric(parameter₁) where Δ represents the comparison 
    between the path minimzing the expected value for criteria₁ and the path minimzing the expected 
    value of criteria₂.

    # Arguments
    - `parameter₁::String`      : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `parameter₂::String`      : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria₁::String`       : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria₂::String`       : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `metric::function`        : statistical function (mean, var, iqr, idr, cv etc.)
    - `weighted::Bool=false`    : if true, distribution is weighted by OD flows

    # Example
    - `sprayΔ("CO₂", "PM"; criteria₁="TD", criteria₂="TT", metric=mean)` returns %Δmean(PM) vs. %Δmean(CO₂) 
       comparing shortest path with the fastest path between OD pairs.
    """
    function sprayΔ(parameter₁, parameter₂; criteria₁, criteria₂, metric, weighted=false)
        X₁ = load(joinpath(@__DIR__, "Results\\$loc - $criteria₁.jld"))["X"]
        X₂ = load(joinpath(@__DIR__, "Results\\$loc - $criteria₂.jld"))["X"]

        Δˣ, Δʸ = Float64[], Float64[]
        for (k, OD) ∈ enumerate(ODs)
            Z₁, Z₂ = X₁[k], X₂[k]

            Pˣ₁ = [0.0 for _ ∈ 1:1000]
            Pˣ₂ = [0.0 for _ ∈ 1:1000]
            for p ∈ C₂P[parameter₁]
                i = findfirst(isequal(p), P)
                Pˣ₁ += Z₁[i] * ℿ[i]
                Pˣ₂ += Z₂[i] * ℿ[i]
            end
            append!(Δˣ, (metric(Pˣ₁)/metric(Pˣ₂)-1)*100)

            Pʸ₁ = [0.0 for _ ∈ 1:1000]
            Pʸ₂ = [0.0 for _ ∈ 1:1000]
            for p ∈ C₂P[parameter₂]
                i = findfirst(isequal(p), P)
                Pʸ₁ += Z₁[i] * ℿ[i]
                Pʸ₂ += Z₂[i] * ℿ[i]
            end
            append!(Δʸ, (metric(Pʸ₁)/metric(Pʸ₂)-1)*100)
        end

        ε = Δʸ./Δˣ
        filter!(x -> (!isnan(x)), ε)
        println("Elasticity wrt $parameter₁")
        println("   ols: $(sum(Δʸ.*Δˣ)/sum(Δˣ.*Δˣ))")
        println("   avg₁: $(mean(ε))")
        println("   avg₂: $(mean(Δʸ)/mean(Δˣ))")

        w = [if weighted Q[i] else 1.0 end for i ∈ eachindex(Q)]
        fig = plot(Δˣ, [Δ*sum(w.*Δʸ.*Δˣ)/sum(w.*Δˣ.*Δˣ) for Δ ∈ Δˣ], linewidth=2,
            color=RGBA(180/255,120/255,70/255), label="")
        fig = scatter!(Δˣ, Δʸ, color=RGBA(70/255,130/255,180/255,50/255),
            markersize=2.5, markerstrokewidth=0.1, label="",
            xlab="% Δ$metric($parameter₁)", ylab="% Δ$metric($parameter₂)",
            tickfont="calibri", guidefont="calibri", legendfont="calibri")
        display(fig)
    end

    
    
    # Reliability plots ────────────────────────────────────────────────────────────────────────────────
    """
        crossreliability(; criteria, parameter, weighted=false)
        
    Returns reliability for the parameter on path minimizing the expected value of the given criteria.

    ### Arguments
    - `criteria::String`        : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `parameter::String`       : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `weighted::Bool=false`    : if true, returns results weighted by demand between OD pairs

    # Example
    - `crossreliability(criterion="TT", parameter="CO₂")` returns CO₂ reliability on the 
      simulated fastest paths between the  OD pairs.
    """
    function crossreliability(; criteria, parameter, weighted=false)
        X = load(joinpath(@__DIR__, "Results\\$loc - $criteria.jld"))["X"]

        l = -0.25
        u = 0.25
        pRange = l:(u-l)/100:u

        Y = [[0.0 for _ ∈ ODs] for _ ∈ 0:100]
        
        for (k, OD) ∈ enumerate(ODs)
            Z = X[k]
            V = [0.0 for _ ∈ 1:1000]
            for p ∈ C₂P[parameter]
                j = findfirst(isequal(p), P)
                V += Z[j] * ℿ[j]
            end

            for (j, p) ∈ enumerate(pRange)
                t = mean(V) * (p + 1)
                Y[j][k] = length(findall(x -> (x ≤ t), V))/length(V)
            end
        end

        w = [if weighted Q[i] else 1.0 end for i ∈ eachindex(Q)]
        fig = plot(xlab="% deviation from the avg. $(C₂L[parameter])",
            ylab="$parameter reliability")
        fig = plot!(tickfont=(14,"calibri"), guidefont=(16,"calibri"),
            legendfont=(16,"calibri"), legend=(0.75, 0.8))
        
        fig = plot!(pRange.*100, [mean(y, weights(w)) for y ∈ Y],
            label="", linewidth=2.5, ylims=(0,1), yticks=0:0.2:1.0)
        display(fig)
    end


    
    # Cost-Benefit analysis ────────────────────────────────────────────────────────────────────────────────
    """
        costbenefit(cost, pollutant, weighted=false)

    Returns log₁₀(Cost/Benefit) value
    Returns ΔC [\$] vs. Δpollutant-cost [\$] scatter plot and ΔC/Δpollutant-cost histogram

    # Arguments
    - `cost::String`            : one of TD, TT, TC
    - `pollutant::String`       : one of CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `weighted::Bool=false`    : if true, returns results weighted by demand between OD pairs
    
    # Example
    - `costbenefit("TT", "CO₂")` returns cost benefit value comparing fastest path with path 
      minimizing expected CO₂ emissions, with cost being increase in travel time and benefits 
      being reduction in CO₂ emissions owing to eco-routing for the pollutant.
    """
    function costbenefit(cost, pollutant; weighted=false)
        X₁ = load(joinpath(@__DIR__, "Results\\$loc - $pollutant.jld"))["X"]
        X₂ = load(joinpath(@__DIR__, "Results\\$loc - $cost.jld"))["X"]


        eᵢ = findfirst(isequal(C₂P[pollutant][1]), P)
        Δᵉ, Δᶜ = Float64[], Float64[]
        for (k, OD) ∈ enumerate(ODs)
            Zᵉ, Zᶜ = X₁[k], X₂[k]

            Eᵉ = mean(Zᵉ[eᵢ] * ℿ[eᵢ])
            Eᶜ = mean(Zᶜ[eᵢ] * ℿ[eᵢ])

            Cᵉ = 0.0
            Cᶜ = 0.0
            for p ∈ C₂P[cost]
                i = findfirst(isequal(p), P)
                Cᵉ += mean(Zᵉ[i] * ℿ[i])
                Cᶜ += mean(Zᶜ[i] * ℿ[i])
            end

            append!(Δᵉ, Eᶜ - Eᵉ)
            append!(Δᶜ, Cᵉ - Cᶜ)
        end

        w = [if weighted Q[i] else 1.0 end for i ∈ eachindex(ODs)]
        Z = [Δᶜ[k]/Δᵉ[k] for k ∈ eachindex(ODs)]
        println("Cost-Benefit: ", log10(sum(Δᶜ, weights(w))/sum(Δᵉ, weights(w))))
        index = findall(x -> (x ≥ 0), Z)
        Z = Z[index]
        w = w[index]
        Δᶜ= Δᶜ[index]
        Δᵉ= Δᵉ[index]
        println("Cost-Benefit: ", mean(log10.(Z), weights(w)))
        #=
        fig = plot(Δᵉ, [Δ*sum(w.*Δᶜ.*Δᵉ)/sum(w.*Δᵉ.*Δᵉ) for Δ ∈ Δᵉ], linewidth=2,
            color=RGBA(180/255,120/255,70/255), label="Avg. Cost/Benefit")
        fig = scatter!(Δᵉ, Δᶜ, color=RGBA(70/255,130/255,180/255,50/255),
            markersize=2.5, markerstrokewidth=0.1, label="")
        fig = plot!(xlab="Δ$(C₂L[pollutant]) (\$)", ylab="Δ$(C₂L[cost]) (\$)",
            tickfont="calibri", guidefont="calibri", legendfont="calibri",
            legend=(0.75, 0.5))
        display(fig)

        fig = histogram(log10.(Z), weights=w, normalize=:pdf, bins=50, linewidth=0.1,
            linecolor=:white, labels="log₁₀(C/B)",
            color=RGBA(70/255,130/255,180/255,200/255))
        fig = vline!([mean(log10.(Z), weights(w))], linewidth=2, labels="mean",
            color=RGBA(180/255,120/255,70/255))
        fig = plot!(xlab="Δ$(C₂L[cost]) (\$)/Δ$(C₂L[pollutant]) (\$)", ylab="density",
            tickfont="calibri", guidefont="calibri", legendfont="calibri",
            legend=(0.75, 0.75))
        display(fig)
        =#
    end



    # Proximity analysis ────────────────────────────────────────────────────────────────────────────────
    """
        proximity_analysis(; pollutants, criteria)
    
    Returns scatter plot for % Δpollutant vs. shortest path travel distance between origin-destination 
    (spatial distance) and a scatter plot for % Δpollutant vs. fastest path travel time between origin-
    destination (temporal distance). Δ represents comparison of shortest or fastest path with the path 
    minimizing expected value of the criteria.

    # Arguments
    - `pollutants::Array{String,1}` : subset of CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria::String`            : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP

    # Example
    - `proximity_analysis(pollutants=["CO₂"], criteria="TC")` returns %change in CO₂ emissions between shortest/fastest
       and least cost path between the OD pairs, plotted against the minimum travel distance/time between these origin-
       destination pairs
    """
    function proximity_analysis(;pollutants, criteria)
        Xₒ = load(joinpath(@__DIR__, "Results\\$loc - $criteria.jld"))["X"]
        Xᵈ = load(joinpath(@__DIR__, "Results\\$loc - TD.jld"))["X"]
        Xᵗ = load(joinpath(@__DIR__, "Results\\$loc - TT.jld"))["X"]
        𝐗 = []
        for (i, pollutant) ∈ enumerate(pollutants)
            push!(𝐗, load(joinpath(@__DIR__, "Results\\$loc - $pollutant.jld"))["X"])
        end

        D = Float64[]
        T = Float64[]
        for (k, OD) ∈ enumerate(ODs)
            Zᵈ = Xᵈ[k]
            Zᵗ = Xᵗ[k]
            append!(D, mean(Zᵈ[1]))
            append!(T, mean(Zᵗ[2]))
        end

        Δ = [Float64[] for _ ∈ pollutants]
        for (i, pollutant) ∈ enumerate(pollutants)
            j = findfirst(isequal(C₂P[pollutant][1]), P)
            for (k, OD) ∈ enumerate(ODs)
                Zₒ, Zᵢ = Xₒ[k], 𝐗[i][k]
                append!(Δ[i], (mean(Zᵢ[j])/mean(Zₒ[j])-1)*100)
            end
        end

        fig = plot()
        for (i, pollutant) ∈ enumerate(pollutants)
            X = 0:1:Int(ceil(maximum(D)))
            Y = [mean(Δ[i][findall(x -> (X[j-1] < x ≤ X[j]), D)]) for j ∈ 2:length(X)]
            fig = plot!(X[2:end], Y, label=pollutant)
        end
        display(fig)

        fig = plot()
        for (i, pollutant) ∈ enumerate(pollutants)
            X = 0:1/60:Int(ceil(maximum(T)))
            Y = [mean(Δ[i][findall(x -> (X[j-1] < x ≤ X[j]), T)]) for j ∈ 2:length(X)]
            fig = plot!(X[2:end], Y, label=pollutant)
        end
        display(fig)
    end
    # ────────────────────────────────────────────────────────────────────────────────
    
    measure(["TD", "TT", "FC", "CH₄", "CO", "CO₂", "NOₓ", "PM", "ROG"], criteria="ROG", metric=mean, weighted=true)
    #costbenefit("TC", "CH₄"; weighted=true)
    return
end
main()