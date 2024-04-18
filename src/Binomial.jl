module Binomial

export Model, compute_value

# preissner_pratt_inv_fn(x) = 

# TODO: Add support for Leisen-Reimer
"""
    Model

Enum encoding available binomial models.
"""
@enum Model begin
    cox_ross_rubinstein
    jarrow_rud
end

"""
    compute_value(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool, n::Int64; model::Model=cox_ross_rubinstein)

Compute the Binomial model fair value `V` of a European option with Strike `K` expiring in `T` years on an underlying with dividend rate `q`, spot value `S` and volatility `vol` given the interest rate `r`.

The argument `is_call` specifies whether the contract is a Call or a Put.

The argument `n` specifies the size of the binomial tree.

The optional argument `model` specifies which Binomial model to use for computing the tree and is by default set to cox_ross_rubinstein (Cox-Ross-Rubinstein).

# Examples
```julia-repl
julia> compute_value(100., 0., 0.05, 0.2, 110., 1., true, 10)
6.099184904016942
```
"""
function compute_value(S::Float64, q::Float64, r::Float64, vol::Float64, K::Float64, T::Float64, is_call::Bool, n::Int64; model::Model=cox_ross_rubinstein)
    dt = T/n
    if model == cox_ross_rubinstein
        u = exp(vol*sqrt(dt))
        d = exp(-vol*sqrt(dt))
    elseif model == jarrow_rud
        u = exp(((r - q) - 0.5*vol^2)*dt + vol*sqrt(dt))
        d = exp(((r - q) - 0.5*vol^2)*dt - vol*sqrt(dt))
    else 
        error("Unsupported binomial model implementation.")
    end

    ST = S*d.^(n:-1:0).*u.^(0:n)

    if is_call
        V = max.(ST .- K, 0.0)
    else
        V = max.(K .- ST, 0.0)
    end

    q = (exp((r - q)*dt) - d)/(u - d)

    for i=n:-1:1
        V = q*V[2:i + 1] + (1 - q)*V[1:i]
    end

    return exp(-r*T)*V[1]
end

end # module Binomial