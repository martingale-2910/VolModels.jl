using Test
using Logging
import QuantTools.BS as BS, QuantTools.Binomial as Bin

@testset "Black-Scholes model tests" verbose=true begin
    S = 100.
    q = 0.0
    r = 0.05
    vol = 0.2
    K = 110.
    T = 1.

    # Call option price
    @test begin
        V_expected = 6.040088129724232
        V_actual = BS.price(S, q, r, vol, K, T, true)
        isapprox(V_actual, V_expected; atol=1e-7)
    end

    # Call-Put parity
    @test begin
        V_call = BS.price(S, q, r, vol, K, T, true)
        V_put = BS.price(S, q, r, vol, K, T, false)
        isapprox(V_call - V_put, S - BS.df(r, T)*K; atol=1e-7)
    end

    # BS Delta, Theta, Gamma
    @test begin
        V_call = BS.price(S, q, r, vol, K, T, true)
        delta = BS.delta(S, q, r, vol, K, T, true)
        theta = BS.theta(S, q, r, vol, K, T, true)
        gamma = BS.gamma(S, q, r, vol, K, T, true)
        isapprox(theta + 0.5*vol^2*S^2*gamma + r*S*delta, r*V_call; atol=1e-7)
    end

    # BS Vega, Vanna, Volga
    @test begin
        gamma = BS.gamma(S, q, r, vol, K, T, true)
        vega = BS.vega(S, q, r, vol, K, T, true)
        isapprox(vega, gamma*S^2*T*vol; atol=1e-7)
        vanna = BS.vanna(S, q, r, vol, K, T, true)
        isapprox(vanna, vega*BS.d2(S, q, r, vol, K, T)/(S*vol*sqrt(T)); atol=1e-7)
        volga = BS.volga(S, q, r, vol, K, T, true)
        isapprox(volga, vega*BS.d1(S, q, r, vol, K, T)*BS.d2(S, q, r, vol, K, T)/(vol); atol=1e-7)
    end

    # BS Implied Vol
    @test begin
        V_actual = BS.price(S, q, r, vol, K, T, true)
        vol_impl = BS.implied_vol(S, q, r, V_actual, K, T, true)
        isapprox(vol_impl, vol; atol=1e-8)
        V_impl = BS.price(S, q, r, vol_impl, K, T, true)
        isapprox(V_impl, V_actual; atol=1e-2)
        @test_logs (:warn, "Option arbitrage boundaries violated.") min_level=Logging.Warn vol_impl = BS.implied_vol(S, q, r, S + 1, K, T, true)
        isequal(vol_impl, NaN64)
        @test_logs (:warn, "Option arbitrage boundaries violated.") min_level=Logging.Warn vol_impl = BS.implied_vol(S, q, r, max(S - K*BS.df(r, T), 0.0) - 1, K, T, true)
        isequal(vol_impl, NaN64)
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

    # Binomial-BS convergence
    @test begin
        bs_price = BS.price(S, q, r, vol, K, T, true)
        bin_price = Bin.price(S, q, r, vol, K, T, true, n)
        isapprox(bin_price, bs_price; atol=1e-2)
        bs_delta = BS.delta(S, q, r, vol, K, T, true)
        bin_delta = Bin.delta(S, q, r, vol, K, T, true, n)
        isapprox(bin_delta, bs_delta; atol=1e-4)
        bs_gamma = BS.gamma(S, q, r, vol, K, T, true)
        bin_gamma = Bin.gamma(S, q, r, vol, K, T, true, n)
        isapprox(bin_gamma, bs_gamma; atol=1e-4)
    end
end