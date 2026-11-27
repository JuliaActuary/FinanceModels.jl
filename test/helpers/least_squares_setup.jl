# Shared by the least-squares test files (least_squares*.jl).
using ForwardDiff, Random, LinearAlgebra

# A bond index: semiannual bonds on a monthly maturity grid, priced off a smooth curve plus a spread per
# issuer, so that no curve reprices them all.
function lsq_index(nb; seed = 11)
    rng = Random.Xoshiro(seed)
    truth = Yield.NelsonSiegel(2.0, 0.058, -0.02, 0.01)
    mats = round.(12 .* (0.5 .+ 19.5 .* rand(rng, nb))) ./ 12
    cpns = round.(800 .* (0.02 .+ 0.05 .* rand(rng, nb))) ./ 800
    spreads = 0.003 .* randn(rng, nb)
    bonds = [Bond.Fixed(c, Periodic(2), m) for (c, m) in zip(cpns, mats)]
    return [Quote(pv(truth + Yield.Constant(Continuous(s)), b), b) for (b, s) in zip(bonds, spreads)]
end

# The knot tenors of the curves fitted to the index
const lsq_knots = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]

# A contract type that `fit` knows nothing about: it wraps a contract and prices as that contract does.
struct LsqOpaque{C} <: FinanceCore.AbstractContract
    inner::C
end
FinanceCore.maturity(c::LsqOpaque) = FinanceCore.maturity(c.inner)
FinanceCore.present_value(m, c::LsqOpaque) = FinanceCore.present_value(m, c.inner)

zr(c) = rate.(knot_rates(c))
# the pricing errors of the quotes `q` on a curve with knot tenors `kn`, written independently of the fit
errors(spline, q, kn) = z -> [pv(ZeroRateCurve(z, kn, spline), x.instrument) - x.price for x in q]
# The fit's certificate, recomputed with ForwardDiff's own Hessian of half the weighted sum of squares:
# positive definite, and a Newton correction of any one knot rate moves no quote's price by more than
# 1e-10 of the quote's size.
function certified(e, z, q, w = ones(length(q)))
    ϕ(x) = sum(w .* abs2.(e(x))) / 2
    H, J, f = ForwardDiff.hessian(ϕ, z), ForwardDiff.jacobian(e, z), e(z)
    δ = H \ ForwardDiff.gradient(ϕ, z)
    size(i) = max(abs(q[i].price), abs(f[i] + q[i].price))
    return isposdef(Symmetric(H)) && all(abs(J[i, k] * δ[k]) <= 1.0e-10 * size(i) for i in axes(J, 1), k in axes(J, 2))
end
