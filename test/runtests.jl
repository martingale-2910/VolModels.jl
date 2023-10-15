using Test
using Printf
using QuantTools: bs_value

St = 100.
r = 0.05
vol = 0.2
K = 110.
t = 0.
T = 1.
is_call = true

@testset "BS value test" begin
    @test begin
        expected_value = 6.040088129724232
        actual_value = bs_value(St, r, vol, K, t, T, is_call)
        isapprox(actual_value, expected_value; atol=1e-7)
    end
end