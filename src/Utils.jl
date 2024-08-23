module Utils

export Φ, φ

using Distributions: Normal, cdf, pdf

N01 = Normal()

Φ(x) = cdf(N01, x)
φ(x) = pdf(N01, x)

end # module Utils