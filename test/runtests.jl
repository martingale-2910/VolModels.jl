using Test
using Printf
import QuantTools.BS as BS, QuantTools.Binomial as Bin

@testset "Black-Scholes model tests" verbose=true begin
    S = 100.
    q = 0.0
    r = 0.05
    vol = 0.2
    K = 110.
    T = 1.
    # Call option value
    @test begin
        expected_value = 6.040088129724232
        actual_value = BS.compute_value(S, q, r, vol, K, T, true)
        isapprox(actual_value, expected_value; atol=1e-7)
    end

    # Call-Put parity
    @test begin
        call = BS.compute_value(S, q, r, vol, K, T, true)
        put = BS.compute_value(S, q, r, vol, K, T, false)
        isapprox(call - put, S - BS.df(r, T)*K; atol=1e-7)
    end

    # BS Delta, Theta, Gamma
    @test begin
        call = BS.compute_value(S, q, r, vol, K, T, true)
        delta = BS.compute_delta(S, q, r, vol, K, T, true)
        theta = BS.compute_theta(S, q, r, vol, K, T, true)
        gamma = BS.compute_gamma(S, q, r, vol, K, T, true)
        isapprox(theta + 0.5*vol^2*S^2*gamma + r*S*delta, r*call; atol=1e-7)
    end

    # BS Vega, Vanna, Volga
    @test begin
        gamma = BS.compute_gamma(S, q, r, vol, K, T, true)
        vega = BS.compute_vega(S, q, r, vol, K, T, true)
        isapprox(vega, gamma*S^2*T*vol; atol=1e-7)
        vanna = BS.compute_vanna(S, q, r, vol, K, T, true)
        isapprox(vanna, vega*BS.d2(S, q, r, vol, K, T)/(S*vol*sqrt(T)); atol=1e-7)
        volga = BS.compute_volga(S, q, r, vol, K, T, true)
        isapprox(volga, vega*BS.d1(S, q, r, vol, K, T)*BS.d2(S, q, r, vol, K, T)/(vol); atol=1e-7)
    end
end

@testset "Binomial model tests" verbose=true begin
    S = 100.
    q = 0.0
    r = 0.05
    vol = 0.2
    K = 110.
    T = 1.
    n = 140
    # BS vs Binomial
    @test begin
        bs_price = BS.compute_value(S, q, r, vol, K, T, true)
        bin_price = Bin.compute_value(S, q, r, vol, K, T, true, n)
        isapprox(bin_price, bs_price; atol=1e-2)
    end
end