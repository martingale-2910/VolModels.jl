# QuantTools - Tools for Quantitative Finance

A just-for-fun project on quantitative finance in Julia.

## Features

**QuantTools** is a `Julia` package that provides simple functionality for **quantitative finance**.

## Example usage

```Julia
using Printf
using QuantTools: bs_value

r = 0.05 # riskless rate of return

S = 100. # spot price of underlying
q = 0.0 # dividend rate of underlying
vol = 0.2 # volatility of underlying

K = 110. # option strike
T = 1. # time-to-maturity
is_call = true # contract right

price = bs_value(S, q, r, vol, K, T, is_call)

@printf('BS fair price of a European option is %f.', price)
```

For more complex examples, please refer to the `/examples/` directory. 

## Dependencies

- `Distributions` package
