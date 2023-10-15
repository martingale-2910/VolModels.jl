module QuantTools

using Distributions: Normal, cdf, pdf

Φ(x) = cdf(Normal(), x)
φ(x) = pdf(Normal(), x)

@doc raw"""
```math
    bs_value(St, r, vol, K, t, T, is_call)

Compute the Black-Scholes fair value at time `t` of an option with Strike `K` expiring at time `T` on an underlying with spot value `St` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_value(100., 0.05, 0.2, 110., 0., 1., true)
6.040088129724232
```
"""
function bs_value(St::Float64, r::Float64, vol::Float64, K::Float64, t::Float64, T::Float64, is_call::Bool)
    ttm = T - t
    d1 = (log(St/K) + (r + 0.5*vol^2)*ttm)/(vol*sqrt(ttm))
    d2 = d1 - vol*sqrt(ttm)
    df = exp(-r*ttm)
    if is_call
        return St*Φ(d1) - K*df*Φ(d2)
    else
        return K*df*Φ(-d2) - St*Φ(-d1)
    end
end

@doc raw"""
```math
    bs_delta(St, r, vol, K, t, T, is_call)

Compute the Black-Scholes delta (`$\frac{\partial V_t}{\partial S_t}$`) at time `t` of an option with Strike `K` expiring at time `T` on an underlying with spot value `St` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_delta(100., 0.05, 0.2, 110., 0., 1., true)
0.44964793063717595
```
"""
function bs_delta(St::Float64, r::Float64, vol::Float64, K::Float64, t::Float64, T::Float64, is_call::Bool)
    ttm = T - t
    d1 = (log(St/K) + (r + 0.5*vol^2)*ttm)/(vol*sqrt(ttm))
    if is_call
        return Φ(d1)
    else
        return -Φ(-d1)
    end
end

@doc raw"""
```math
    bs_vega(St, r, vol, K, t, T, is_call)

Compute the Black-Scholes vega (`$\frac{\partial V_t}{\partial \sigma}$`) at time `t` of an option with Strike `K` expiring at time `T` on an underlying with spot value `St` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_vega(100., 0.05, 0.2, 110., 0., 1., true)
39.57604803881934
```
"""
function bs_vega(St::Float64, r::Float64, vol::Float64, K::Float64, t::Float64, T::Float64, is_call::Bool)
    ttm = T - t
    df = exp(-r*ttm)
    if is_call  # Both expressions are actually equal
        d1 = (log(St/K) + (r + 0.5*vol^2)*ttm)/(vol*sqrt(ttm))
        return St*φ(d1)*sqrt(ttm)
    else
        d2 = (log(St/K) + (r - 0.5*vol^2)*ttm)/(vol*sqrt(ttm))
        return K*df*φ(d2)*sqrt(ttm)
    end
end

@doc raw"""
```math
    bs_theta(St, r, vol, K, t, T, is_call)

Compute the Black-Scholes theta (`$\frac{\partial V_t}{\partial (t - T)}$`) at time `t` of an option with Strike `K` expiring at time `T` on an underlying with spot value `St` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_theta(100., 0.05, 0.2, 110., 0., 1., true)
39.06095301313599
```
"""
function bs_theta(St::Float64, r::Float64, vol::Float64, K::Float64, t::Float64, T::Float64, is_call::Bool)
    ttm = T - t
    d1 = (log(St/K) + (r + 0.5*vol^2)*ttm)/(vol*sqrt(ttm))
    d2 = d1 - vol*sqrt(ttm)
    df = exp(-r*ttm)
    if is_call
        return -(St*φ(d1)*vol)/(2*sqrt(ttm)) - r*K*df*Φ(d2) + St*Φ(d1)
    else
        return -(St*φ(d1)*vol)/(2*sqrt(ttm)) + r*K*df*Φ(-d2) - St*Φ(-d1)
    end
end

@doc raw"""
```math
    bs_rho(St, r, vol, K, t, T, is_call)

Compute the Black-Scholes rho (`$\frac{\partial V_t}{\partial r}$`) at time `t` of an option with Strike `K` expiring at time `T` on an underlying with spot value `St` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_rho(100., 0.05, 0.2, 110., 0., 1., true)
38.92470493399336
```
"""
function bs_rho(St::Float64, r::Float64, vol::Float64, K::Float64, t::Float64, T::Float64, is_call::Bool)
    ttm = T - t
    d2 = (log(St/K) + (r - 0.5*vol^2)*ttm)/(vol*sqrt(ttm))
    df = exp(-r*ttm)
    if is_call
        return K*ttm*df*Φ(d2)
    else
        return -K*ttm*df*Φ(-d2)
    end
end

@doc raw"""
```math
    bs_gamma(St, r, vol, K, t, T, is_call)

Compute the Black-Scholes gamma (`$\frac{\partial^2 V_t}{\partial \left(S_t\right)^2}$`) at time `t` of an option with Strike `K` expiring at time `T` on an underlying with spot value `St` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_gamma(100., 0.05, 0.2, 110., 0., 1., true)
0.019788024019409666
```
"""
function bs_gamma(St::Float64, r::Float64, vol::Float64, K::Float64, t::Float64, T::Float64, is_call::Bool)
    ttm = T - t
    df = exp(-r*ttm)
    if is_call  # Both expressions are actually equal
        d1 = (log(St/K) + (r + 0.5*vol^2)*ttm)/(vol*sqrt(ttm))
        return φ(d1)/(St*vol*sqrt(ttm))
    else
        d2 = (log(St/K) + (r - 0.5*vol^2)*ttm)/(vol*sqrt(ttm))
        return K*df*(φ(d2))/(St^2*vol*sqrt(ttm))
    end
end

@doc raw"""
```math
    bs_vanna(St, r, vol, K, t, T, is_call)

Compute the Black-Scholes vanna (`$\frac{\partial^2 V_t}{\partial S_t \partial \sigma}$`) at time `t` of an option with Strike `K` expiring at time `T` on an underlying with spot value `St` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_vanna(100., 0.05, 0.2, 110., 0., 1., true)
0.6461797033399725
```
"""
function bs_vanna(St::Float64, r::Float64, vol::Float64, K::Float64, t::Float64, T::Float64, is_call::Bool)
    ttm = T - t
    d1 = (log(St/K) + (r + 0.5*vol^2)*ttm)/(vol*sqrt(ttm))
    vega = bs_vega(St, r, vol, K, t, T, is_call)
    return (1 - d1/(vol*sqrt(ttm)))*vega/St
end

@doc raw"""
```math
    bs_volga(St, r, vol, K, t, T, is_call)

Compute the Black-Scholes volga (`$\frac{\partial^2 V_t}{\partial \sigma^2}$`) at time `t` of an option with Strike `K` expiring at time `T` on an underlying with spot value `St` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_volga(100., 0.05, 0.2, 110., 0., 1., true)
8.17746223872001
```
"""
function bs_volga(St::Float64, r::Float64, vol::Float64, K::Float64, t::Float64, T::Float64, is_call::Bool)
    ttm = T - t
    d1 = (log(St/K) + (r + 0.5*vol^2)*ttm)/(vol*sqrt(ttm))
    d2 = d1 - vol*sqrt(ttm)
    vega = bs_vega(St, r, vol, K, t, T, is_call)
    return vega*d1*d2/vol
end

@doc raw"""
```math
    bs_implied_vol(St, r, Vt, K, t, T, is_call)

Compute the Black-Scholes implied volatility at time `t` of an option with strike `K` and expiration in `T` quoted at price `Vt` on an underlying with spot price `St` and interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

Uses the Newton-Rhapson method to solve the equation `$V_t = bs_value(S_t, r, \sigma_{impl}, K, t, T, is_call)$` and converge the the correct value.

# Examples
```julia-repl
julia> bs_implied_vol(100., 0.05, 6.040088129724232, 110., 0., 1., true)
0.19999999999999993
```
"""
function bs_implied_vol(St::Float64, r::Float64, Vt::Float64, K::Float64, t::Float64, T::Float64, is_call::Bool)
    vol = 0.3
    tol = 1e-8
    max_iter = 100
    for _ in 1:max_iter
        price = bs_value(St, r, vol, K, t, T, is_call)
        vega = bs_vega(St, r, vol, K, t, T, is_call)
        diff = Vt - price
        if abs(diff) < tol
            return vol
        end
        vol = vol + diff/vega
    end
    return vol
end

end # module QuantTools