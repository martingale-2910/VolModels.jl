module Options

export arbitrage, arbitrage_bounds

"""
    arbitrage(V:Float64, S::Float64, K::Float64, D::Float64, DF::Float64, is_call::Bool)

Compute whether a european option with a given price `V` constitutes an arbitrage given its strike `K`, the future value of its dividends `D`, the underlying's spot price `S` and the to-maturity discount factor `DF`.

# Examples
```julia-repl
julia> arbitrage(101., 100., 110., 0., exp(-0.05*1.0), true)
true
```
"""
function arbitrage(V::Float64, S::Float64, K::Float64, D::Float64, DF::Float64, is_call::Bool)
    bounds = arbitrage_bounds(S, K, D, DF, is_call)
    return !(bounds[1] <= V && V <= bounds[2])
end

"""
    arbitrage_bounds(V:Float64, S::Float64, K::Float64, D::Float64, DF::Float64, is_call::Bool)

Compute european option arbitrage bounds for an option with strike `K`, future value of dividends `D`, on an underlying with spot price `S` and given the to-maturity discount factor `DF`.

# Examples
```julia-repl
julia> arbitrage_bounds(5., 100., 110., 0., exp(-0.05*1.0), true)
(0.0, 100.0)
```
"""
function arbitrage_bounds(S::Float64, K::Float64, D::Float64, DF::Float64, is_call::Bool)
    if is_call
        return (max(S - K*DF - D*DF, 0.0), S)
    else
        return (max(K*DF - S + D*DF, 0.0), K*DF)
    end
end

end # module Options