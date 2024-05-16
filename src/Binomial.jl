module Binomial

export Model, compute_value

# preissner_pratt_inv_fn(x) = 

# TODO: Add support for Leisen-Reimer
"""
    Model

Exposes available Binomial models implementations.
"""
@enum Model begin
    cox_ross_rubinstein
    jarrow_rud
end

"""
    calibrate(r::Float64, q::Float64, vol::Float64, dt::Float64, model::Model)

Calibrates the Binomial model based on interest rate `r`, dividend rate `q` and volatility `vol`.

Computes and returns magnitudes of up `u` and down `d` movements in the model.

Predicates based on desired Binomial model type `model`.

# Examples
```julia-repl
julia> calibrate(0.05, 0., 0.2, 0.1, cox_ross_rubinstein)
(1.0652883920946086, 0.9387129414165152)
```
"""
function calibrate(r::Float64, q::Float64, vol::Float64, dt::Float64, model::Model)
    if model == cox_ross_rubinstein
        u = exp(vol*sqrt(dt))
        d = exp(-vol*sqrt(dt))
    elseif model == jarrow_rud
        u = exp(((r - q) - 0.5*vol^2)*dt + vol*sqrt(dt))
        d = exp(((r - q) - 0.5*vol^2)*dt - vol*sqrt(dt))
    else 
        error("Unsupported Binomial model implementation.")
    end
    return u, d
end

"""
    compute_value(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool, n::Int64; model::Model=cox_ross_rubinstein)

Compute the Binomial model fair value `V` of a European option with Strike `K` expiring in `T` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

The argument `n` specifies the size of the binomial tree.

The optional argument `model` specifies which Binomial model to use for computing the tree and is by default set to `cox_ross_rubinstein` (Cox-Ross-Rubinstein).

# Examples
```julia-repl
julia> compute_value(100., 0., 0.05, 0.2, 110., 1., true, 10)
6.099184904016942
```
"""
function compute_value(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool, n::Int64; model::Model=cox_ross_rubinstein)
    dt = T/n
    u, d = calibrate(r, q, vol, dt, model)

    ST = S*d.^(n:-1:0).*u.^(0:n)

    V = is_call ? max.(ST .- K, 0.0) : max.(K .- ST, 0.0)

    q = (exp((r - q)*dt) - d)/(u - d)

    for i=n:-1:1
        V = q*V[2:i + 1] + (1 - q)*V[1:i]
    end

    return exp(-r*T)*V[1]
end

"""
    compute_delta(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool, n::Int64; model::Model=cox_ross_rubinstein)

Compute the Binomial model delta `Δ` (`∂V/∂S`) of a European option with Strike `K` expiring in `T` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

The argument `n` specifies the size of the binomial tree.

The optional argument `model` specifies which Binomial model to use for computing the tree and is by default set to `cox_ross_rubinstein`` (Cox-Ross-Rubinstein).

# Examples
```julia-repl
julia> compute_delta(100., 0., 0.05, 0.2, 110., 1., true, 10)
0.44072597450896794
```
"""
function compute_delta(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool, n::Int64; model::Model=cox_ross_rubinstein)
    dt = T/n
    u, d = calibrate(r, q, vol, dt, model)

    ST = S*d.^(n:-1:0).*u.^(0:n)

    V = is_call ? max.(ST .- K, 0.0) : max.(K .- ST, 0.0)

    q = (exp((r - q)*dt) - d)/(u - d)

    for i=n:-1:2
        V = q*V[2:i + 1] + (1 - q)*V[1:i]
    end

    V = exp(-r*dt*(n - 1))*V

    return (V[2] - V[1])/(S*(u - d))
end

"""
    compute_gamma(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool, n::Int64; model::Model=cox_ross_rubinstein)

Compute the Binomial model gamma `Γ` (`∂²V/∂S²`) of a European option with Strike `K` expiring in `T` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

The argument `n` specifies the size of the binomial tree.

The optional argument `model` specifies which Binomial model to use for computing the tree and is by default set to `cox_ross_rubinstein`` (Cox-Ross-Rubinstein).

# Examples
```julia-repl
julia> compute_gamma(100., 0., 0.05, 0.2, 110., 1., true, 10)
0.02050590337391038
```
"""
function compute_gamma(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool, n::Int64; model::Model=cox_ross_rubinstein)
    dt = T/n
    u, d = calibrate(r, q, vol, dt, model)

    ST = S*d.^(n:-1:0).*u.^(0:n)

    V = is_call ? max.(ST .- K, 0.0) : max.(K .- ST, 0.0)

    q = (exp((r - q)*dt) - d)/(u - d)

    for i=n:-1:3
        V = q*V[2:i + 1] + (1 - q)*V[1:i]
    end

    V = exp(-r*dt*(n - 2))*V

    return ((V[3] - V[2])/(S*u*(u - d)) - (V[2] - V[1])/(S*d*(u - d)))/(0.5*(S*(u + d)*(u - d)))
end

end # module Binomial