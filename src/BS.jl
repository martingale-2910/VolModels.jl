module BS

export price, delta, vega, theta, rho, gamma, vanna, volga, implied_vol, d1, d2, df, arbitrage

using Logging

include("Utils.jl")
using .Utils: Φ, φ

include("Options.jl")
using .Options: arbitrage_bounds

"""
    d1(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64)

Compute the `d₁` score used to determine the risk-adjusted probability that the option will be in-the-money.

# Examples
```julia-repl
julia> d1(100., 0., 0.05, 0.2, 110, 1.)
-0.12655089902162442
```
"""
function d1(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64)
    return (log(S/K) + ((r - q) + 0.5*vol^2)*T)/(vol*sqrt(T))
end

"""
    d2(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64)

Compute the `d₂` score used to determine the probability that the option will `not` be in-the-money.

# Examples
```julia-repl
julia> d2(100., 0., 0.05, 0.2, 110, 1.)
-0.32655089902162443
```
"""
function d2(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64)
    return (log(S/K) + ((r - q) - 0.5*vol^2)*T)/(vol*sqrt(T))
end

"""
    df(r::Float64, t::Float64)

Compute the time `t` discount factor in the Black-Scholes model with interest rate `r`.

# Examples
```julia-repl
julia> df(0.05, 1.)
0.951229424500714
```
"""
function df(r::Float64, t::Float64)
    return exp(-r*t)
end

"""
    price(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)

Compute the Black-Scholes price `V` of a European option with Strike `K` expiring in `T` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> price(100., 0., 0.05, 0.2, 110., 1., true)
6.040088129724232
```
"""
function price(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    if is_call
        return S*df(q, T)*Φ(d1(S, q, r, vol, K, T)) - K*df(r, T)*Φ(d2(S, q, r, vol, K, T))
    else
        return K*df(r, T)*Φ(-d2(S, q, r, vol, K, T)) - S*df(q, T)*Φ(-d1(S, q, r, vol, K, T))
    end
end

"""
    delta(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)

Compute the Black-Scholes delta `Δ` (`∂V/∂S`) of a European option with Strike `K` expiring in `T` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> delta(100., 0., 0.05, 0.2, 110., 1., true)
0.44964793063717595
```
"""
function delta(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    if is_call
        return df(q, T)*Φ(d1(S, q, r, vol, K, T))
    else
        return -df(q, T)*Φ(-d1(S, q, r, vol, K, T))
    end
end

"""
    vega(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)

Compute the Black-Scholes vega `ν` (`∂V/∂σ`) of a European option with Strike `K` expiring in `T` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> vega(100., 0., 0.05, 0.2, 110., 1., true)
39.57604803881934
```
"""
function vega(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    if is_call  # Both expressions are actually equal
        return S*df(q, T)*sqrt(T)*φ(d1(S, q, r, vol, K, T))
    else
        return K*df(r, T)*sqrt(T)*φ(d2(S, q, r, vol, K, T))
    end
end

"""
    theta(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)

Compute the Black-Scholes theta `θ` (`∂V/∂T`) of a European option with Strike `K` expiring in `T` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> theta(100., 0., 0.05, 0.2, 110., 1., true)
-5.903840050581602
```
"""
function theta(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    if is_call
        return -df(q, T)*S*φ(d1(S, q, r, vol, K, T))*vol/(2*sqrt(T)) - r*K*df(r, T)*Φ(d2(S, q, r, vol, K, T)) + q*S*df(q, T)*Φ(d1(S, q, r, vol, K, T))
    else
        return -df(q, T)*S*φ(d1(S, q, r, vol, K, T))*vol/(2*sqrt(T)) + r*K*df(r, T)*Φ(d2(S, q, r, vol, K, T)) - q*S*df(q, T)*Φ(d1(S, q, r, vol, K, T))
    end
end

"""
    rho(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)

Compute the Black-Scholes rho `ρ` (`∂V/∂r`) of a European option with Strike `K` expiring in `T` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> rho(100., 0., 0.05, 0.2, 110., 1., true)
38.92470493399336
```
"""
function rho(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    if is_call
        return K*T*df(r, T)*Φ(d2(S, q, r, vol, K, T))
    else
        return -K*T*df(r, T)*Φ(-d2(S, q, r, vol, K, T))
    end
end

"""
    gamma(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)

Compute the Black-Scholes gamma `Γ` (`∂²V/∂S²`) of a European option with Strike `K` expiring in `T` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> gamma(100., 0., 0.05, 0.2, 110., 1., true)
0.019788024019409666
```
"""
function gamma(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    if is_call  # Both expressions are actually equal
        return df(q, T)*φ(d1(S, q, r, vol, K, T))/(S*vol*sqrt(T))
    else
        return df(r, T)*K*φ(d2(S, q, r, vol, K, T))/(S^2*vol*sqrt(T))
    end
end

"""
    vanna(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)

Compute the Black-Scholes vanna (`∂²V/∂S∂σ`) of a European option with Strike `K` expiring in `T` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> vanna(100., 0., 0.05, 0.2, 110., 1., true)
0.6461797033399724
```
"""
function vanna(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    return -df(q, T)*φ(d1(S, q, r, vol, K, T))*d2(S, q, r, vol, K, T)/vol
end

"""
    volga(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)

Compute the Black-Scholes volga (`∂²V/∂σ²`) of a European option with Strike `K` expiring in `T` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> volga(100., 0., 0.05, 0.2, 110., 1., true)
8.17746223872001
```
"""
function volga(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    return df(q, T)*S*sqrt(T)*φ(d1(S, q, r, vol, K, T))*d1(S, q, r, vol, K, T)*d2(S, q, r, vol, K, T)/vol
end

"""
    implied_vol(S::Float64, q::Float64, r::Float64, V::Float64, K::Float64, T::Float64, is_call::Bool; max_iter::Int64=100, tol::Float64=1e-2, eps::Float64=1e-4, verbose::Bool=true)

Compute the Black-Scholes implied volatility of a European option with strike `K`, expiring in `T` years and quoted at price `V` on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

Uses the bisection method to solve the equation `V = compute_value(S, q, r, σ, K, T, is_call)` for `σ` .

Computes at most `max_iter` (200 by default) iterations with a tolerance of `tol` (1e-2 by default) for the option value and `eps` (1e-4 by default) for the volatility.

# Examples
```julia-repl
julia> implied_vol(100., 0., 0.05, 6.040088129724232, 110., 1., true)
0.19999999999999993
```
"""
function implied_vol(S::Float64, q::Float64, r::Float64, V::Float64, K::Float64, T::Float64, is_call::Bool; max_iter::Int64=200, tol::Float64=1e-2, eps::Float64=1e-4, verbose::Bool=true)
    arbitrage = arbitrage_bounds(S, K, df(-q, T), df(r, T), is_call)

    if V < arbitrage[1] || arbitrage[2] < V
        if verbose
            @warn "Option arbitrage boundaries violated."
        end
        return V < arbitrage[1] ? 0.0 : Inf64
    end

    vol_low = 0.0
    vol_high = 2.0

    for _ in 0:max_iter
        vol_mid = 0.5*(vol_high + vol_low)
        V_mid = price(S, q, r, vol_mid, K, T, is_call)

        if abs(V_mid - V) < tol || abs(0.5*(vol_high - vol_low)) < eps
            return vol_mid
        elseif V_mid < V
            vol_low = vol_mid
        else
            vol_high = vol_mid
        end
    end

    if !(abs(V_mid - V) < tol || abs(0.5*(vol_high - vol_low)) < eps) && verbose
        @warn "Implied volatility did not converge." "abs(V_mid - V) = abs($V_mid - $V) = $(abs(V_mid - V))" "abs(vol_high - vol_low) = abs($vol_high - $vol_low) = $(abs(vol_high - vol_low))"
    end
    return vol_mid
end

end # module BS