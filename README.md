# VolModels.jl - Volatility models in Julia

## Introduction
A just-for-fun project, **VolModels.jl** is a `Julia` package that provides simple functionality for working with volatility models when considering pricing options, hedging them and analyzing their PnL.

The inspiration for this project was [The Volatility Smile](https://www.amazon.com/Volatility-Smile-Wiley-Finance/dp/1118959167).

It exposes the following functionalities:
- `Black-Scholes` model:
    - pricing;
    - greeks;
    - implied volatility.
- `Binomial` model:
    - pricing;
    - greeks (delta, gamma).


## Example usage

- ### Option valuation
```Julia
import VolModels.BS.price as bs_price, VolModels.Binomial.price as bin_price

r = 0.05 # riskless rate of return

S = 100. # spot price of underlying
q = 0. # dividend rate of underlying
vol = 0.2 # volatility of underlying

K = 110. # option strike
T = 1. # time-to-maturity
is_call = true # contract right

println("BS price of a European option is $(bs_price(S, q, r, vol, K, T, is_call)).")

n = 200

println("Binomial price of a European option is $(bin_price(S, q, r, vol, K, T, is_call, n)).")
```

For more complex examples, please refer to the `/examples/` directory. 

## Dependencies

- `Distributions`
- `Logging`