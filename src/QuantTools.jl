module QuantTools

using Distributions: Normal, cdf, pdf

Φ(x) = cdf(Normal(), x)
φ(x) = pdf(Normal(), x)
d1(S, r, vol, K, T) = (log(S/K) + (r + 0.5*vol^2)*T)/(vol*sqrt(T))
d2(S, r, vol, K, T) = (log(S/K) + (r - 0.5*vol^2)*T)/(vol*sqrt(T))
df(r, t) = exp(-r*t)

@doc raw"""
    bs_value(S, r, vol, K, T, is_call)

Compute the Black-Scholes fair value `V` of a European option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_value(100., 0.05, 0.2, 110., 1., true)
6.040088129724232
```
"""
function bs_value(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    if is_call
        return S*Φ(d1(S, r, vol, K, T)) - K*df(r, T)*Φ(d2(S, r, vol, K, T))
    else
        return K*df(r, T)*Φ(-d2(S, r, vol, K, T)) - S*Φ(-d1(S, r, vol, K, T))
    end
end

@doc raw"""
    bs_delta(S, r, vol, K, T, is_call)

Compute the Black-Scholes delta (```math $\frac{\partial V}{\partial S}$```) of a European option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_delta(100., 0.05, 0.2, 110., 1., true)
0.44964793063717595
```
"""
function bs_delta(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    if is_call
        return Φ(d1(S, r, vol, K, T))
    else
        return -Φ(-d1(S, r, vol, K, T))
    end
end

@doc raw"""
    bs_vega(S, r, vol, K, T, is_call)

Compute the Black-Scholes vega (```math $\frac{\partial V}{\partial \sigma}$```) of a European option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_vega(100., 0.05, 0.2, 110., 1., true)
39.57604803881934
```
"""
function bs_vega(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    if is_call  # Both expressions are actually equal
        return S*sqrt(T)*φ(d1(S, r, vol, K, T))
    else
        return K*df(r, T)*φ(d2(S, r, vol, K, T))*sqrt(T)
    end
end

@doc raw"""
    bs_theta(S, r, vol, K, T, is_call)

Compute the Black-Scholes theta (```math $\frac{\partial V}{\partial T}$```) of a European option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_theta(100., 0.05, 0.2, 110., 1., true)
-5.903840050581602
```
"""
function bs_theta(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    if is_call
        return -(S*φ(d1(S, r, vol, K, T))*vol)/(2*sqrt(T)) - r*K*df(r, T)*Φ(d2(S, r, vol, K, T))
    else
        return -(S*φ(d1(S, r, vol, K, T))*vol)/(2*sqrt(T)) + r*K*df(r, T)*Φ(-d2(S, r, vol, K, T))
    end
end

@doc raw"""
    bs_rho(S, r, vol, K, T, is_call)

Compute the Black-Scholes rho (```math $\frac{\partial V}{\partial r}$```) of a European option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_rho(100., 0.05, 0.2, 110., 1., true)
38.92470493399336
```
"""
function bs_rho(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    if is_call
        return K*T*df(r, T)*Φ(d2(S, r, vol, K, T))
    else
        return -K*T*df(r, T)*Φ(-d2(S, r, vol, K, T))
    end
end

@doc raw"""
    bs_gamma(S, r, vol, K, T, is_call)

Compute the Black-Scholes gamma (```math $\frac{\partial^2 V}{\partial \left(S\right)^2}$`) of a European option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_gamma(100., 0.05, 0.2, 110., 1., true)
0.019788024019409666
```
"""
function bs_gamma(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    if is_call  # Both expressions are actually equal
        return φ(d1(S, r, vol, K, T))/(S*vol*sqrt(T))
    else
        return K*df(r, T)*(φ(d2(S, r, vol, K, T)))/(S*vol*sqrt(T))
    end
end

@doc raw"""
    bs_vanna(S, r, vol, K, T, is_call)

Compute the Black-Scholes vanna (```math $\frac{\partial^2 V}{\partial S \partial \sigma}$```) of a European option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_vanna(100., 0.05, 0.2, 110., 1., true)
0.6461797033399724
```
"""
function bs_vanna(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    return -φ(d1(S, r, vol, K, T))*d2(S, r, vol, K, T)/vol
end

@doc raw"""
    bs_volga(S, r, vol, K, T, is_call)

Compute the Black-Scholes volga (```math $\frac{\partial^2 V}{\partial \sigma^2}$```) of a European option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_volga(100., 0.05, 0.2, 110., 1., true)
8.17746223872001
```
"""
function bs_volga(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    return S*sqrt(T)*φ(d1(S, r, vol, K, T))*d1(S, r, vol, K, T)*d2(S, r, vol, K, T)/vol
end

@doc raw"""
    bs_implied_vol(S, r, V, K, T, is_call)

Compute the Black-Scholes implied volatility of a European option with strike `K`, expiring in `T` years and quoted at price `V` on an underlying with spot price `S` and interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

Uses the Newton-Rhapson method to solve the equation ```math $V = \text{bs_value}(S, r, \sigma_{impl}, K, T, \text{is_call})$``` and converge the the correct value.

Assumes an initial starting volatility `vol` of 30% and computes at most `max_iter` (100 by default) iterations with a tolerance of `tol` (set to 1e-8 by default).

# Examples
```julia-repl
julia> bs_implied_vol(100., 0.05, 6.040088129724232, 110., 1., true)
0.19999999999999993
```
"""
function bs_implied_vol(S::Float64, r::Float64, V::Float64, K::Float64, T::Float64, is_call::Bool, vol::Float64=0.2, max_iter::Int64=100, tol::Float64=1e-8)
    for _ in 1:max_iter
        price = bs_value(S, r, vol, K, T, is_call)
        vega = bs_vega(S, r, vol, K, T, is_call)
        diff = V - price
        if abs(diff) < tol
            return vol
        end
        vol = vol + diff/vega
    end
    return vol
end

end # module QuantTools