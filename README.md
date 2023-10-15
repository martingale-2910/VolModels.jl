# QuantTools - Tools for Quantitative Finance

A just-for-fun project on quantitative finance in Julia.

## Features

**QuantTools** is a `Julia` package that provides simple functionality for **quantitative finance**.

## Usage

```Julia
using Printf
using QuantTools: bs_value

r = 0.05 # riskless rate of return
vol = 0.2 # volatility of underlying security

K = 110. # option strike
t = 0.
T = 1. # option time-to-maturity (TTM)
is_call = true

St = 100. # initial underlying price

price = bs_value(St, r, vol, K, t, T, is_call)

@printf('BS fair price is %f.', price)
```

## Dependencies

- `Distributions` package
