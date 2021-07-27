using Plots
#using StatsPlots
using JLD
#using ColorSchemes
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
        "TC" => ["TD", "TT", "FC"], "GHG" => ["CH4", "CO2", "N2O", "ROG"], "CP" => ["CO", "NOx", "PM", "SOx"])  # Criteria to parameters
    
    C₂L = Dict{String, String}("TD" => "Travel Distance (TD)", "TT" => "Travel Time (TT)",
            "FC" => "Fuel Consumption (FC)", "CH₄" => "Methane emissions (CH₄)",
            "CO"  => "Carbon Monoxide emissions (CO)", "CO₂" => "Carbon Dioxide emissions (CO₂)",
            "N₂O" => "Nitrous Oxide emissions (N₂O)", "NOₓ" => "Nitrogen Oxide emissions (NOₓ)",
            "PM"  => "Particulate Matter emissions (PM)", "ROG" => "Reactive Organic Gases (ROG)",
            "SOₓ" => "Sulphur Oxide emissions (SOₓ)", "TC" => "Travel Cost (TC)")                               # Criteria to Label

    ℿ = Float64[]                           # Cost parameters
    ODs = Array{Int64,1}[]                  # Origin-Destination
    Q = Float64[]                           # OD demand

    # Coefficients file
    coefFile = joinpath(@__DIR__, "Network\\$network\\coef.csv")
    csv₁ = CSV.File(coefFile)
    df₁ = DataFrame(csv₁)
    for r in 1:nrow(df₁) append!(ℿ, df₁[r,5]) end

    # Demand file
    dmndFile = joinpath(@__DIR__, "Network\\$network\\demand.csv")
    csv₂ = CSV.File(dmndFile)
    df₂ = DataFrame(csv₂)
    origins = df₂[!,1]
    destinations = parse.(Int64, String.(names(df₂))[2:end])
    for r in origins for s in destinations push!(ODs, [r, s]) end end
    for r in 1:nrow(df₂) for c in 2:ncol(df₂) append!(Q, df₂[r,c]) end end


    # Simulator ────────────────────────────────────────────────────────────────────────────────
    """
        sim(;criteria)

    Finds path between all OD pairs minimizing the expected value of the given criteria

    ### Arguments
    - `criteria::String`    : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    """
    function sim(;criteria)
        X = [Array{Float64,1}[] for _ in 1:length(ODs)]   # Simulated data from ssp
        for (k, OD) in enumerate(ODs)
            r, s = OD
            Z, _, _ = ssp(r, s, network=network, parameter=C₂P[criteria], numsims=1000)
            append!(X[k], Z)
        end
        save(joinpath(@__DIR__, "Results\\$loc\\$loc - $criteria - HCT.jld"), "X", X)
    end



    # Measure ────────────────────────────────────────────────────────────────────────────────
    """
        measure(parameters; criteria, metric, weighted=false)

    Returns metric on parameter from parameters for paths established for given criteria 
        
    ### Arguments
    - `parameters::Array{String,1}`  : vector on TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria::String`             : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `metric::function`             : statistical function (mean, var, iqr, idr, cv etc.)
    - `weighted::Bool=false`         : if true results are weighted by demand between OD pairs
    """
    function measure(parameters; criteria, metric, weighted=false)
        X = load(joinpath(@__DIR__, "Results\\SP\\$loc\\$loc - $criteria.jld"))["X"]
        for parameter ∈ parameters
            Y = [0.0 for _ in 1:length(ODs)]
            w = [0.0 for _ in 1:length(ODs)]
            for (k, _) in enumerate(ODs)
                Z = [0.0 for _ in 1:1000]
                for p in C₂P[parameter]
                    i = findfirst(x -> (x == p), P)
                    Z += X[k][i] * ℿ[i]
                end
                Y[k] = metric(Z)
                w[k] = weighted ? Q[k] : 1.0
            end
            println("$metric($parameter) on least $criteria path: ", (mean(Y), weights(w)))
        end
    end

    # Δ plots ────────────────────────────────────────────────────────────────────────────────
    """
        Δ(parameter; criteria₁, criteria₂, metric, weighted=false)

    Returns %change in parameter metric value in criteria₂ paths relative to criteria₁ paths 

    ### Arguments
    - `parameter::String`       : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria₁::String`       : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria₂::String`       : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `metric::function`        : statistical function (mean, var, iqr, idr, cv etc.)
    - `weighted::Bool=false`    : if true results are weighted by demand between OD pairs
    """
    function Δ(parameter; criteria₁, criteria₂, metric, weighted=false)
        X₁ = load(joinpath(@__DIR__, "Results\\$loc - $criteria₁.jld"))["X"]
        X₂ = load(joinpath(@__DIR__, "Results\\$loc - $criteria₂.jld"))["X"]
        
        Δ = [0.0 for _ in ODs]
        w = [0.0 for _ in ODs]
        for (k, _) in enumerate(ODs)
            Z₁, Z₂ = X₁[k], X₂[k]

            P₁ = [0.0 for _ in 1:1000]
            P₂ = [0.0 for _ in 1:1000]
            for p in C₂P[parameter]
                i = findfirst(x -> (x == p), P)
                P₁ += Z₁[i] * ℿ[i]
                P₂ += Z₂[i] * ℿ[i]
            end
            
            δ = (metric(P₂)/metric(P₁) - 1) * 100
            Δ[k] = δ
            w[k] = weighted ? Q[k] : 1.0
        end
        println("% change in $metric($parameter): ", mean(Δ, weights(w)))
    end

    """
        distributeΔ(parameter; criteria₁, criteria₂, ODpair)
    
    Absolute distribution of the parameter given OD pair and path criteria

    ### Arguments
    - `parameter::String`    : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria₁::String`    : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria₂::String`    : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `ODpair::Array{Int,1}` : Vector of OD pairs
    """
    function distributeΔ(parameter; criteria₁, criteria₂, ODpair)
        k = findfirst(x -> (x == ODpair), ODs)
        X₁ = load(joinpath(@__DIR__, "Results\\$loc\\$loc - $criteria₁.jld"))["X"]
        X₂ = load(joinpath(@__DIR__, "Results\\$loc\\$loc - $criteria₂.jld"))["X"]

        Z₁, Z₂ = X₁[k], X₂[k]
        P₁ = [0.0 for _ in 1:1000]
        P₂ = [0.0 for _ in 1:1000]
        for p in C₂P[parameter]
            i = findfirst(x -> (x == p), P)
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
        fig = vline!([mean(P₂)], color=RGBA(64/255,31/255,127/255), linewidth=2, label="")
        display(fig)

        # LCP: RGB 64, 31, 127
        # FP : RGB 127, 64, 31
        # SP : RGB 127, 31, 94
    end

    """
        sprayΔ(X, Y; criteria₁, criteria₂, metric, weighted=false)

    Scatter plot of % Δmetric(Y) vs % Δmetric(X) between paths with different criteria

    ### Arguments
    - `X::String`               : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `Y::String`               : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria₁::String`       : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria₂::String`       : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `metric::function`        : statistical function (mean, var, iqr, idr, cv etc.)
    - `weighted::Bool=false`    : if true, distribution is weighted by OD flows
    """
    function sprayΔ(X, Y; criteria₁, criteria₂, metric, weighted=false)
        X₁ = load(joinpath(@__DIR__, "Results\\$loc\\$loc - $criteria₁.jld"))["X"]
        X₂ = load(joinpath(@__DIR__, "Results\\$loc\\$loc - $criteria₂.jld"))["X"]

        Δˣ, Δʸ = Float64[], Float64[]
        for (k, OD) in enumerate(ODs)
            Z₁, Z₂ = X₁[k], X₂[k]

            Pˣ₁ = [0.0 for _ in 1:1000]
            Pˣ₂ = [0.0 for _ in 1:1000]
            for p in C₂P[parameter₁]
                i = findfirst(x -> (x == p), P)
                Pˣ₁ += Z₁[i] * ℿ[i]
                Pˣ₂ += Z₂[i] * ℿ[i]
            end
            append!(Δˣ, (metric(Pˣ₁)/metric(Pˣ₂)-1)*100)

            Pʸ₁ = [0.0 for _ in 1:1000]
            Pʸ₂ = [0.0 for _ in 1:1000]
            for p in C₂P[parameter₂]
                i = findfirst(x -> (x == p), P)
                Pʸ₁ += Z₁[i] * ℿ[i]
                Pʸ₂ += Z₂[i] * ℿ[i]
            end
            append!(Δʸ, (metric(Pʸ₁)/metric(Pʸ₂)-1)*100)
        end

        ε = Δʸ./Δˣ
        filter!(x -> (!isnan(x)), ε)
        println("Elasticity wrt $X")
        println("   ols: $(sum(Δʸ.*Δˣ)/sum(Δˣ.*Δˣ))")
        println("   avg₁: $(mean(ε))")
        println("   avg₂: $(mean(Δʸ)/mean(Δˣ))")

        w = [if weighted Q[i] else 1.0 end for i in 1:length(Q)]
        fig = plot(Δˣ, [Δ*sum(w.*Δʸ.*Δˣ)/sum(w.*Δˣ.*Δˣ) for Δ in Δˣ], linewidth=2,
            color=RGBA(180/255,120/255,70/255), label="")
        fig = scatter!(Δˣ, Δʸ, color=RGBA(70/255,130/255,180/255,50/255),
            markersize=2.5, markerstrokewidth=0.1, label="",
            xlab="% Δ$metric($X)", ylab="% Δ$metric($Y)",
            tickfont="calibri", guidefont="calibri", legendfont="calibri")
        display(fig)
    end


    # Reliability plots ────────────────────────────────────────────────────────────────────────────────
    """
        crossreliability(; criterion, parameter, reliability="parameter", weighted=false)
        
    ???

    ### Arguments
    - `parameter::String`               : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criterion::Array{String, 1}`     : subset of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `reliability::String`             : criterion or parameter
    - `weighted::Bool=false`            : if true, returns results weighted by demand between OD pairs
    """
    function crossreliability(; criterion, parameter, reliability="parameter", weighted=false)
        Xₒ = load(joinpath(@__DIR__, "Results\\$loc\\$loc - $parameter.jld"))["X"]
        𝐗 = []
        for (i, criteria) in enumerate(criterion)
            push!(𝐗, load(joinpath(@__DIR__, "Results\\$loc\\$loc - $criteria.jld"))["X"])
        end

        l = -0.25
        u = 0.25
        pRange = l:(u-l)/100:u

        𝐘 = [[[0.0 for _ in ODs] for _ in 0:100] for _ in criterion]
        for (i, criteria) in enumerate(criterion)
            if reliability == "criterion" x, X = criteria, Xₒ end
            if reliability == "parameter" x, X = parameter, 𝐗[i] end

            for (k, OD) in enumerate(ODs)
                Z = X[k]
                P = [0.0 for _ in 1:1000]
                for p in C₂P[x]
                    j = findfirst(x -> (x == p), P)
                    P += Z[j] * ℿ[j]
                end

                for (j, p) in enumerate(pRange)
                    t = mean(P) * (p + 1)
                    𝐘[i][j][k] = length(findall(x -> (x <= t), P))/length(P)
                end
            end
        end
        #legend = ["SP", "FP", "LCP", "LEP - $(criterion[1])"]
        #rgba = [RGBA(127/255, 31/255, 94/255, 100/255), RGBA(127/255, 64/255, 31/255, 100/255),
        #        RGBA(64/255, 31/255, 127/255, 100/255), RGBA(31/255, 127/255, 64/255, 100/255)]
        legend = ["SP", "LEP - $(criterion[2])", "LEP - $(criterion[3])", "LEP - $(criterion[4])",
                "LEP - $(criterion[5])", "LEP - $(criterion[6])", "LEP - $(criterion[7])"]

        w = [if weighted Q[i] else 1.0 end for i in 1:length(Q)]
        fig = plot(xlab="% deviation from the avg. $(C₂L[parameter])",
            ylab="$parameter reliability")
        fig = plot!(tickfont=(14,"calibri"), guidefont=(16,"calibri"),
            legendfont=(16,"calibri"), legend=(0.75, 0.8))
        for (k, Y) in enumerate(𝐘)
            fig = plot!(pRange.*100, [mean(y, weights(w)) for y in Y],
                label=legend[k], color=palette(:Paired_7)[k], linewidth=2.5,
                ylims=(0,1), yticks=0:0.2:1.0)
        end
        display(fig)
    end


    # Cost-Benefit analysis ────────────────────────────────────────────────────────────────────────────────
    """
        costbenefit(cost, pollutant, weighted=false)

    Returns log₁₀(Cost/Benefit) value
    Returns ΔC [\$] vs. Δpollutant-cost [\$] scatter plot and ΔC/Δpollutant-cost histogram

    ### Arguments
    - `cost::String`            : one of TD, TT, TC
    - `pollutant::String`       : one of CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `weighted::Bool=false`    : if true, returns results weighted by demand between OD pairs
    """
    function costbenefit(cost, pollutant; weighted=false)
        X₁ = load(joinpath(@__DIR__, "Results\\$loc\\$loc - $pollutant.jld"))["X"]
        X₂ = load(joinpath(@__DIR__, "Results\\$loc\\$loc - $cost.jld"))["X"]


        eᵢ = findfirst(x -> (x == C₂P[pollutant][1]), P)
        Δᵉ, Δᶜ = Float64[], Float64[]
        for (k, OD) in enumerate(ODs)
            Zᵉ, Zᶜ = X₁[k], X₂[k]

            Eᵉ = mean(Zᵉ[eᵢ] * ℿ[eᵢ])
            Eᶜ = mean(Zᶜ[eᵢ] * ℿ[eᵢ])

            Cᵉ = 0.0
            Cᶜ = 0.0
            for p in parameterize[cost]
                i = findfirst(x -> (x == p), P)
                Cᵉ += mean(Zᵉ[i] * ℿ[i])
                Cᶜ += mean(Zᶜ[i] * ℿ[i])
            end

            append!(Δᵉ, Eᶜ - Eᵉ)
            append!(Δᶜ, Cᵉ - Cᶜ)
        end

        w = [if weighted Q[i] else 1.0 end for i in 1:length(Q)]
        Z = [Δᶜ[k]/Δᵉ[k] for k in 1:length(ODs)]
        index = findall(x -> (x ≥ 0), Z)
        Z = Z[index]
        w = w[index]
        println("Cost-Benefit: ", mean(log10.(Z), weights(w)))


        fig = plot(Δᵉ, [Δ*sum(w.*Δᶜ.*Δᵉ)/sum(w.*Δᵉ.*Δᵉ) for Δ in Δᵉ], linewidth=2,
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
    end


    # Proximity analysis ────────────────────────────────────────────────────────────────────────────────
    """
        proximity_analysis(; pollutants, criteria)
    
    Returns a scatter plot for % Δpollutant vs. shortest path travel distance between origin-destination
    and a scatter plot for % Δpollutant vs. fastest path travel time between origin-destination

    ### Arguments
    - `pollutants::Array{String,1}` : subset of CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    - `criteria::String`            : one of TD, TT, FC, CH₄, CO, CO₂, N₂O, NOₓ, PM, ROG, SOₓ, TC, GHG, CP
    """
    function proximity_analysis(;pollutants, criteria)
        Xₒ = load(joinpath(@__DIR__, "Results\\$loc\\$loc - $criteria.jld"))["X"]
        Xᵈ = load(joinpath(@__DIR__, "Results\\$loc\\$loc - TD.jld"))["X"]
        Xᵗ = load(joinpath(@__DIR__, "Results\\$loc\\$loc - TT.jld"))["X"]
        𝐗 = []
        for (i, pollutant) in enumerate(pollutants)
            push!(𝐗, load("Network\\$networkName\\analysis\\$loc\\$loc - $pollutant.jld")["X"])
        end

        D = Float64[]
        T = Float64[]
        for (k, OD) in enumerate(ODs)
            Zᵈ = Xᵈ[k]
            Zᵗ = Xᵗ[k]
            append!(D, mean(Zᵈ[1]))
            append!(T, mean(Zᵗ[2]))
        end

        Δ = [Float64[] for _ in pollutants]
        for (i, pollutant) in enumerate(pollutants)
            j = findfirst(x -> (x == C₂P[pollutant][1]), P)
            for (k, OD) in enumerate(ODs)
                Zₒ, Zᵢ = Xₒ[k], 𝐗[i][k]
                append!(Δ[i], (mean(Zᵢ[j])/mean(Zₒ[j])-1)*100)
            end
        end

        fig = plot()
        for (i, pollutant) in enumerate(pollutants)
            X = 0:1:Int(ceil(maximum(D)))
            Y = [mean(Δ[i][findall(x -> (X[j-1] < x <= X[j]), D)]) for j in 2:length(X)]
            fig = plot!(X[2:end], Y, label=pollutant)
        end
        display(fig)

        fig = plot()
        for (i, pollutant) in enumerate(pollutants)
            X = 0:1/60:Int(ceil(maximum(T)))
            Y = [mean(Δ[i][findall(x -> (X[j-1] < x <= X[j]), T)]) for j in 2:length(X)]
            fig = plot!(X[2:end], Y, label=pollutant)
        end
        display(fig)
    end
    # ────────────────────────────────────────────────────────────────────────────────

    pmt  = "TC"    
    c₁   = "TT"
    c₂   = "ROG"
    Δ(pmt; criteria₁=c₁, criteria₂=c₂, metric=mean, weighted=true)
end
main()

#= ────────────────────────────────────────────────────────────────────────────────
# TODO: 
1. Complete doc strings of some unfinished functions
2. Update file directories
──────────────────────────────────────────────────────────────────────────────── =#
