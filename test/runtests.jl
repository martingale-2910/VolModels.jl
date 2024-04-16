using Test
using Printf
using QuantTools.BS: compute_value, compute_delta, compute_theta, compute_vega, compute_gamma, compute_vanna, compute_volga, d1, d2, df

S = 100.
q = 0.0
r = 0.05
vol = 0.2
K = 110.
T = 1.

@testset "Black-Scholes model tests" verbose=true begin
    # Call option value
    @test begin
        expected_value = 6.040088129724232
        actual_value = compute_value(S, q, r, vol, K, T, true)
        isapprox(actual_value, expected_value; atol=1e-7)
    end

    # Call-Put parity
    @test begin
        call = compute_value(S, q, r, vol, K, T, true)
        put = compute_value(S, q, r, vol, K, T, false)
        isapprox(call - put, S - df(r, T)*K; atol=1e-7)
    end

    # BS Delta, Theta, Gamma
    @test begin
        call = compute_value(S, q, r, vol, K, T, true)
        delta = compute_delta(S, q, r, vol, K, T, true)
        theta = compute_theta(S, q, r, vol, K, T, true)
        gamma = compute_gamma(S, q, r, vol, K, T, true)
        isapprox(theta + 0.5*vol^2*S^2*gamma + r*S*delta, r*call; atol=1e-7)
    end

    # BS Vega, Vanna, Volga
    @test begin
        gamma = compute_gamma(S, q, r, vol, K, T, true)
        vega = compute_vega(S, q, r, vol, K, T, true)
        isapprox(vega, gamma*S^2*T*vol; atol=1e-7)
        vanna = compute_vanna(S, q, r, vol, K, T, true)
        isapprox(vanna, vega*d2(S, q, r, vol, K, T)/(S*vol*sqrt(T)); atol=1e-7)
        volga = compute_volga(S, q, r, vol, K, T, true)
        isapprox(volga, vega*d1(S, q, r, vol, K, T)*d2(S, q, r, vol, K, T)/(vol); atol=1e-7)
    end
end