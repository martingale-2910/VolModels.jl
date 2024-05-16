using Test
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
        V_expected = 6.040088129724232
        V_actual = BS.compute_value(S, q, r, vol, K, T, true)
        isapprox(V_actual, V_expected; atol=1e-7)
    end

    # Call-Put parity
    @test begin
        V_call = BS.compute_value(S, q, r, vol, K, T, true)
        V_put = BS.compute_value(S, q, r, vol, K, T, false)
        isapprox(V_call - V_put, S - BS.df(r, T)*K; atol=1e-7)
    end

    # BS Delta, Theta, Gamma
    @test begin
        V_call = BS.compute_value(S, q, r, vol, K, T, true)
        delta = BS.compute_delta(S, q, r, vol, K, T, true)
        theta = BS.compute_theta(S, q, r, vol, K, T, true)
        gamma = BS.compute_gamma(S, q, r, vol, K, T, true)
        isapprox(theta + 0.5*vol^2*S^2*gamma + r*S*delta, r*V_call; atol=1e-7)
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

    # BS Implied Vol
    @test begin
        V_actual = BS.compute_value(S, q, r, vol, K, T, true)
        vol_impl = BS.compute_implied_vol(S, q, r, V_actual, K, T, true)
        isapprox(vol_impl, vol; atol=1e-8)
        V_impl = BS.compute_value(S, q, r, vol_impl, K, T, true)
        isapprox(V_impl, V_actual; atol=1e-2)
        vol_impl = BS.compute_implied_vol(S, q, r, S + 1, K, T, true)
        isequal(vol_impl, NaN64)
        vol_impl = BS.compute_implied_vol(S, q, r, max(S - K*BS.df(r, T), 0.0) - 1, K, T, true)
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
        bs_price = BS.compute_value(S, q, r, vol, K, T, true)
        bin_price = Bin.compute_value(S, q, r, vol, K, T, true, n)
        isapprox(bin_price, bs_price; atol=1e-2)
        bs_delta = BS.compute_delta(S, q, r, vol, K, T, true)
        bin_delta = Bin.compute_delta(S, q, r, vol, K, T, true, n)
        isapprox(bin_delta, bs_delta; atol=1e-4)
        bs_gamma = BS.compute_gamma(S, q, r, vol, K, T, true)
        bin_gamma = Bin.compute_gamma(S, q, r, vol, K, T, true, n)
        isapprox(bin_gamma, bs_gamma; atol=1e-4)
    end
end