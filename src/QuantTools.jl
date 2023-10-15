module QuantTools

using Distributions: Normal, cdf, pdf

Φ(x) = cdf(Normal(), x)
φ(x) = pdf(Normal(), x)

@doc raw"""
```math
    bs_value(S, r, vol, K, T, is_call)

Compute the Black-Scholes fair value `V` of an option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_value(100., 0.05, 0.2, 110., 1., true)
6.040088129724232
```
"""
function bs_value(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    d1 = (log(S/K) + (r + 0.5*vol^2)*T)/(vol*sqrt(T))
    d2 = d1 - vol*sqrt(T)
    df = exp(-r*T)
    if is_call
        return S*Φ(d1) - K*df*Φ(d2)
    else
        return K*df*Φ(-d2) - S*Φ(-d1)
    end
end

@doc raw"""
```math
    bs_delta(S, r, vol, K, T, is_call)

Compute the Black-Scholes delta (`$\frac{\partial V}{\partial S}$`) of an option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_delta(100., 0.05, 0.2, 110., 1., true)
0.44964793063717595
```
"""
function bs_delta(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    d1 = (log(S/K) + (r + 0.5*vol^2)*T)/(vol*sqrt(T))
    if is_call
        return Φ(d1)
    else
        return -Φ(-d1)
    end
end

@doc raw"""
```math
    bs_vega(S, r, vol, K, T, is_call)

Compute the Black-Scholes vega (`$\frac{\partial V}{\partial \sigma}$`) of an option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_vega(100., 0.05, 0.2, 110., 1., true)
39.57604803881934
```
"""
function bs_vega(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    df = exp(-r*T)
    if is_call  # Both expressions are actually equal
        d1 = (log(S/K) + (r + 0.5*vol^2)*T)/(vol*sqrt(T))
        return S*φ(d1)*sqrt(T)
    else
        d2 = (log(S/K) + (r - 0.5*vol^2)*T)/(vol*sqrt(T))
        return K*df*φ(d2)*sqrt(T)
    end
end

@doc raw"""
```math
    bs_theta(S, r, vol, K, T, is_call)

Compute the Black-Scholes theta (`$\frac{\partial V}{\partial T}$`) of an option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_theta(100., 0.05, 0.2, 110., 1., true)
39.06095301313599
```
"""
function bs_theta(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    d1 = (log(S/K) + (r + 0.5*vol^2)*T)/(vol*sqrt(T))
    d2 = d1 - vol*sqrt(T)
    df = exp(-r*T)
    if is_call
        return -(S*φ(d1)*vol)/(2*sqrt(T)) - r*K*df*Φ(d2) + S*Φ(d1)
    else
        return -(S*φ(d1)*vol)/(2*sqrt(T)) + r*K*df*Φ(-d2) - S*Φ(-d1)
    end
end

@doc raw"""
```math
    bs_rho(S, r, vol, K, T, is_call)

Compute the Black-Scholes rho (`$\frac{\partial V}{\partial r}$`) of an option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_rho(100., 0.05, 0.2, 110., 1., true)
38.92470493399336
```
"""
function bs_rho(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    d2 = (log(S/K) + (r - 0.5*vol^2)*T)/(vol*sqrt(T))
    df = exp(-r*T)
    if is_call
        return K*T*df*Φ(d2)
    else
        return -K*T*df*Φ(-d2)
    end
end

@doc raw"""
```math
    bs_gamma(S, r, vol, K, T, is_call)

Compute the Black-Scholes gamma (`$\frac{\partial^2 V}{\partial \left(S\right)^2}$`) of an option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_gamma(100., 0.05, 0.2, 110., 1., true)
0.019788024019409666
```
"""
function bs_gamma(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    df = exp(-r*T)
    if is_call  # Both expressions are actually equal
        d1 = (log(S/K) + (r + 0.5*vol^2)*T)/(vol*sqrt(T))
        return φ(d1)/(S*vol*sqrt(T))
    else
        d2 = (log(S/K) + (r - 0.5*vol^2)*T)/(vol*sqrt(T))
        return K*df*(φ(d2))/(S^2*vol*sqrt(T))
    end
end

@doc raw"""
```math
    bs_vanna(S, r, vol, K, T, is_call)

Compute the Black-Scholes vanna (`$\frac{\partial^2 V}{\partial S \partial \sigma}$`) of an option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_vanna(100., 0.05, 0.2, 110., 1., true)
0.6461797033399725
```
"""
function bs_vanna(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    T = T - t
    d1 = (log(S/K) + (r + 0.5*vol^2)*T)/(vol*sqrt(T))
    vega = bs_vega(S, r, vol, K, T, is_call)
    return (1 - d1/(vol*sqrt(T)))*vega/S
end

@doc raw"""
```math
    bs_volga(S, r, vol, K, T, is_call)

Compute the Black-Scholes volga (`$\frac{\partial^2 V}{\partial \sigma^2}$`) of an option with Strike `K` expiring in `T` years on an underlying with spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_volga(100., 0.05, 0.2, 110., 1., true)
8.17746223872001
```
"""
function bs_volga(S::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool)
    d1 = (log(S/K) + (r + 0.5*vol^2)*T)/(vol*sqrt(T))
    d2 = d1 - vol*sqrt(T)
    vega = bs_vega(S, r, vol, K, T, is_call)
    return vega*d1*d2/vol
end

@doc raw"""
```math
    bs_implied_vol(S, r, V, K, T, is_call)

Compute the Black-Scholes implied volatility of an option with strike `K`, expiring in `T` years and quoted at price `V` on an underlying with spot price `S` and interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

Uses the Newton-Rhapson method to solve the equation `$V = \text{bs_value}(S, r, \sigma_{impl}, K, T, \text{is_call})$` and converge the the correct value.

# Examples
```julia-repl
julia> bs_implied_vol(100., 0.05, 6.040088129724232, 110., 1., true)
0.19999999999999993
```
"""
function bs_implied_vol(S::Float64, r::Float64, V::Float64, K::Float64, T::Float64, is_call::Bool)
    vol = 0.3
    tol = 1e-8
    max_iter = 100
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