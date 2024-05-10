# QuantTools - Tools for Quantitative Finance

A just-for-fun project on quantitative finance in Julia.

## Features

**QuantTools** is a `Julia` package that provides simple functionality for **quantitative finance**.

It exposes the following functionalities:
- `Black-Scholes` model:
    - valuation;
    - implied volatility;
    - greeks.
- `Binomial` model:
    - valuation.

## Example usage

- ### Black-Scholes option valuation
```Julia
import QuantTools.BS.compute_value as compute_bs_value, QuantTools.Binomial.compute_value as compute_bin_value

r = 0.05 # riskless rate of return

S = 100. # spot price of underlying
q = 0. # dividend rate of underlying
vol = 0.2 # volatility of underlying

K = 110. # option strike
T = 1. # time-to-maturity
is_call = true # contract right

println("BS fair price of a European option is $(compute_bs_value(S, q, r, vol, K, T, is_call)).")

n = 10

println("Binomial fair price of a European option is $(compute_bin_value(S, q, r, vol, K, T, is_call, n)).")
```

For more complex examples, please refer to the `/examples/` directory. 

## Dependencies

- `Distributions`
- `Logging`