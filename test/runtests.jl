using Test
using Printf
using QuantTools: bs_value, bs_delta, bs_theta, bs_vega, bs_gamma, bs_vanna, bs_volga, d1, d2

S = 100.
q = 0.0
r = 0.05
vol = 0.2
K = 110.
T = 1.

@testset "BS formulas test" begin
    # Call option value
    @test begin
        expected_value = 6.040088129724232
        actual_value = bs_value(S, q, r, vol, K, T, true)
        isapprox(actual_value, expected_value; atol=1e-7)
    end

    # Call-Put parity
    @test begin
        call = bs_value(S, q, r, vol, K, T, true)
        put = bs_value(S, q, r, vol, K, T, false)
        isapprox(call - put, S - exp(-r*T)*K; atol=1e-7)
    end

    # BS Delta, Theta, Gamma
    @test begin
        call = bs_value(S, q, r, vol, K, T, true)
        delta = bs_delta(S, q, r, vol, K, T, true)
        theta = bs_theta(S, q, r, vol, K, T, true)
        gamma = bs_gamma(S, q, r, vol, K, T, true)
        isapprox(theta + 0.5*vol^2*S^2*gamma + r*S*delta, r*call; atol=1e-7)
    end

    # BS Vega, Vanna, Volga
    @test begin
        gamma = bs_gamma(S, q, r, vol, K, T, true)
        vega = bs_vega(S, q, r, vol, K, T, true)
        isapprox(vega, gamma*S^2*T*vol; atol=1e-7)
        vanna = bs_vanna(S, q, r, vol, K, T, true)
        isapprox(vanna, vega*d2(S, q, r, vol, K, T)/(S*vol*sqrt(T)); atol=1e-7)
        volga = bs_volga(S, q, r, vol, K, T, true)
        isapprox(volga, vega*d1(S, q, r, vol, K, T)*d2(S, q, r, vol, K, T)/(vol); atol=1e-7)
    end
end