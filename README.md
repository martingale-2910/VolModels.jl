# QuantTools - Tools for Quantitative Finance

A just-for-fun project on quantitative finance in Julia.

## Features

**QuantTools** is a `Julia` package that provides simple functionality for **quantitative finance**.

## Usage

```Julia
using Printf
using QuantTools: bs_value

r = 0.05 # riskless rate of return
vol = 0.2 # volatility of underlying

S = 100. # spot price of underlying

K = 110. # option strike
T = 1. # time-to-maturity
is_call = true


price = bs_value(S, r, vol, K, T, is_call)

@printf('BS fair price is %f.', price)
```

## Dependencies

- `Distributions` package
