module QuantTools

using Distributions: Normal, cdf, pdf

Φ(x) = cdf(Normal(), x)
φ(x) = pdf(Normal(), x)
d1(S, q, r, vol, K, tau) = (log(S/K) + ((r - q) + 0.5*vol^2)*tau)/(vol*sqrt(tau))
d2(S, q, r, vol, K, tau) = (log(S/K) + ((r - q) - 0.5*vol^2)*tau)/(vol*sqrt(tau))
df(r, t) = exp(-r*t)

"""
    bs_value(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)

Compute the Black-Scholes fair value `V` of a European option with Strike `K` expiring in `tau` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_value(100., 0.05, 0.2, 110., 1., true)
6.040088129724232
```
"""
function bs_value(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)
    if is_call
        return S*df(q, tau)*Φ(d1(S, q, r, vol, K, tau)) - K*df(r, tau)*Φ(d2(S, q, r, vol, K, tau))
    else
        return K*df(r, tau)*Φ(-d2(S, q, r, vol, K, tau)) - S*df(q, tau)*Φ(-d1(S, q, r, vol, K, tau))
    end
end

"""
    bs_delta(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)

Compute the Black-Scholes delta (`` \\frac{\\partial V}{\\partial S}``) of a European option with Strike `K` expiring in `tau` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_delta(100., 0.05, 0.2, 110., 1., true)
0.44964793063717595
```
"""
function bs_delta(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)
    if is_call
        return df(q, tau)*Φ(d1(S, q, r, vol, K, tau))
    else
        return -df(q, tau)*Φ(-d1(S, q, r, vol, K, tau))
    end
end

"""
    bs_vega(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)

Compute the Black-Scholes vega (``\\frac{\\partial V}{\\partial \\sigma}``) of a European option with Strike `K` expiring in `tau` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_vega(100., 0.0, 0.05, 0.2, 110., 1., true)
39.57604803881934
```
"""
function bs_vega(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)
    if is_call  # Both expressions are actually equal
        return S*df(q, tau)*sqrt(tau)*φ(d1(S, q, r, vol, K, tau))
    else
        return K*df(r, tau)*sqrt(tau)*φ(d2(S, q, r, vol, K, tau))
    end
end

"""
    bs_theta(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)

Compute the Black-Scholes theta (``\\frac{\\partial V}{\\partial T}``) of a European option with Strike `K` expiring in `tau` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_theta(100., 0.0, 0.05, 0.2, 110., 1., true)
-5.903840050581602
```
"""
function bs_theta(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)
    if is_call
        return -df(q, tau)*S*φ(d1(S, q, r, vol, K, tau))*vol/(2*sqrt(tau)) - r*K*df(r, tau)*Φ(d2(S, q, r, vol, K, tau)) + q*S*df(q, tau)*Φ(d1(S, q, r, vol, K, tau))
    else
        return -df(q, tau)*S*φ(d1(S, q, r, vol, K, tau))*vol/(2*sqrt(tau)) + r*K*df(r, tau)*Φ(d2(S, q, r, vol, K, tau)) - q*S*df(q, tau)*Φ(d1(S, q, r, vol, K, tau))
    end
end

"""
    bs_rho(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)

Compute the Black-Scholes rho (``\\frac{\\partial V}{\\partial r}``) of a European option with Strike `K` expiring in `tau` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_rho(100., 0.0, 0.05, 0.2, 110., 1., true)
38.92470493399336
```
"""
function bs_rho(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)
    if is_call
        return K*T*df(r, tau)*Φ(d2(S, q, r, vol, K, tau))
    else
        return -K*T*df(r, tau)*Φ(-d2(S, q, r, vol, K, tau))
    end
end

"""
    bs_gamma(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)

Compute the Black-Scholes gamma (``\\frac{\\partial^2 V}{\\partial \\left(S\\right)^2}``) of a European option with Strike `K` expiring in `tau` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_gamma(100., 0.0, 0.05, 0.2, 110., 1., true)
0.019788024019409666
```
"""
function bs_gamma(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)
    if is_call  # Both expressions are actually equal
        return df(q, tau)*φ(d1(S, q, r, vol, K, tau))/(S*vol*sqrt(tau))
    else
        return df(r, tau)*K*φ(d2(S, q, r, vol, K, tau))/(S^2*vol*sqrt(tau))
    end
end

"""
    bs_vanna(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)

Compute the Black-Scholes vanna (``\\frac{\\partial^2 V}{\\partial S \\partial \\sigma}``) of a European option with Strike `K` expiring in `tau` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_vanna(100., 0.0, 0.05, 0.2, 110., 1., true)
0.6461797033399724
```
"""
function bs_vanna(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)
    return -df(q, tau)*φ(d1(S, q, r, vol, K, tau))*d2(S, q, r, vol, K, tau)/vol
end

"""
    bs_volga(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)

Compute the Black-Scholes volga (``\\frac{\\partial^2 V}{\\partial \\sigma^2}``) of a European option with Strike `K` expiring in `tau` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

# Examples
```julia-repl
julia> bs_volga(100., 0.0, 0.05, 0.2, 110., 1., true)
8.17746223872001
```
"""
function bs_volga(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, tau::Float64, is_call::Bool)
    return df(q, tau)*S*sqrt(tau)*φ(d1(S, q, r, vol, K, tau))*d1(S, q, r, vol, K, tau)*d2(S, q, r, vol, K, tau)/vol
end

"""
    bs_implied_vol(S::Float64, q::Float64, r::Float64, V::Float64, K::Float64, tau::Float64, is_call::Bool, vol::Float64=0.3, max_iter::Int64=100, tol::Float64=1e-8)

Compute the Black-Scholes implied volatility of a European option with strike `K`, expiring in `tau` years and quoted at price `V` on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

Uses the Newton-Rhapson method to solve the equation ``V = \\text{bs_value}(S, q, r, \\sigma_{impl}, K, tau, \\text{is_call})`` and converge to the correct value.

Assumes an initial starting volatility `vol` (30% by default) and computes at most `max_iter` (100 by default) iterations with a tolerance of `tol` (1e-8 by default).

# Examples
```julia-repl
julia> bs_implied_vol(100., 0.0, 0.05, 6.040088129724232, 110., 1., true)
0.19999999999999993
```
"""
function bs_implied_vol(S::Float64, q::Float64, r::Float64, V::Float64, K::Float64, tau::Float64, is_call::Bool, vol::Float64=0.3, max_iter::Int64=100, tol::Float64=1e-8)
    for _ in 1:max_iter
        price = bs_value(S, q, r, vol, K, tau, is_call)
        vega = bs_vega(S, q, r, vol, K, tau, is_call)
        diff = V - price
        if abs(diff) < tol
            return vol
        end
        vol = vol + diff/vega
    end
    return vol
end

end # module QuantTools