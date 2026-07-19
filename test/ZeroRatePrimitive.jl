# Tests for the continuous-zero-rate primitive refactor (FinanceModels 5.7):
# every curve type's `discount` is consistent with its `zero`, composition happens
# in zero-rate space, `zero(c, 0)` is finite (regression for a former 0/0 → NaN),
# and `forward` matches its log-discount definition.

@testset "Zero-rate primitive" begin
    mats = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    zrates = [0.01, 0.015, 0.02, 0.025, 0.03, 0.032]

    spl = Yield.Spline(Spline.Linear(), mats, zrates)           # Yield.Spline (linear)
    splc = Yield.Spline(Spline.Cubic(), mats, zrates .+ 0.004)  # Yield.Spline (cubic)
    mc = Yield.MonotoneConvex(zrates, mats)                     # Yield.MonotoneConvex
    cc = Yield.Constant(Continuous(0.03))
    cper = Yield.Constant(0.04)                                  # Periodic under the hood
    ns = Yield.NelsonSiegel(2.0, 0.03, -0.01, 0.01)
    nss = Yield.NelsonSiegelSvensson(1.5, 3.0, 0.03, -0.02, 0.01, 0.005)
    cpr = Yield.CairnsPritchard(0.5, 3.0, 0.03, -0.01, -0.005)
    comp_add = spl + splc
    comp_sub = spl - splc
    scaled = spl * 0.79
    shifted = spl + ((z, t) -> z + Continuous(0.01))

    curves = [
        ("Spline(linear)", spl), ("Spline(cubic)", splc), ("MonotoneConvex", mc),
        ("Constant(cont)", cc), ("Constant(per)", cper),
        ("NelsonSiegel", ns), ("NSS", nss), ("CairnsPritchard", cpr),
        ("Composite(+)", comp_add), ("Composite(-)", comp_sub),
        ("Scaled", scaled), ("TenorShift", shifted),
    ]

    @testset "discount ≡ exp(-z·t)  ($name)" for (name, c) in curves
        for t in (0.25, 0.5, 1.0, 2.0, 3.7, 5.0, 10.0, 20.0, 30.0)
            z = FinanceCore.rate(zero(c, t))
            @test exp(-z * t) ≈ discount(c, t) rtol = 1e-12
            @test discount(zero(c, t), t) ≈ discount(c, t) rtol = 1e-12
        end
    end

    @testset "discount(c, 0) == 1 exactly  ($name)" for (name, c) in curves
        @test discount(c, 0.0) == 1.0
    end

    @testset "zero(c, 0) is finite (regression: was 0/0 → NaN)" begin
        for (name, c) in (("Spline", spl), ("Constant(cont)", cc),
            ("Constant(per)", cper), ("Composite(+)", comp_add), ("Scaled", scaled),
            ("MonotoneConvex", mc), ("NelsonSiegel", ns), ("NSS", nss), ("CairnsPritchard", cpr))
            @test !isnan(FinanceCore.rate(zero(c, 0.0)))
        end
        # flat curve: zero(·, 0) is its own continuous rate
        @test FinanceCore.rate(zero(cc, 0.0)) ≈ 0.03
        # continuity of the zero curve into t = 0
        @test FinanceCore.rate(zero(spl, 1e-9)) ≈ FinanceCore.rate(zero(spl, 0.0)) atol = 1e-6
        @test FinanceCore.rate(zero(comp_add, 1e-9)) ≈ FinanceCore.rate(zero(comp_add, 0.0)) atol = 1e-6
    end

    @testset "NS/NSS array tenors preserve scalar semantics" begin
        ts = [0.0, 0.5, 2.0, 10.0]
        for c in (ns, nss)
            @test zero(c, ts) == zero.(Ref(c), ts)
            @test discount(c, ts) == discount.(Ref(c), ts)
        end
    end

    @testset "discount(c, 0) preserves curve numeric type" begin
        ns_big = Yield.NelsonSiegel(big"2.0", big"0.03", big"-0.01", big"0.01")
        scaled_big = Yield.Constant(big"0.03") * big"0.79"
        @test discount(ns_big, 0.0) isa BigFloat
        @test discount(scaled_big, 0.0) isa BigFloat
    end

    @testset "composition is additive in zero-rate space" begin
        for t in (0.5, 1.0, 5.0, 20.0)
            z1 = FinanceCore.rate(zero(spl, t))
            z2 = FinanceCore.rate(zero(splc, t))
            @test FinanceCore.rate(zero(comp_add, t)) ≈ z1 + z2
            @test FinanceCore.rate(zero(comp_sub, t)) ≈ z1 - z2
            @test FinanceCore.rate(zero(scaled, t)) ≈ z1 * 0.79
            @test discount(comp_add, t) ≈ discount(spl, t) * discount(splc, t)
        end
    end

    @testset "forward matches log-discount definition  ($name)" for (name, c) in curves
        for (from, to) in ((0.0, 1.0), (1.0, 2.0), (2.0, 5.0), (5.0, 10.0))
            expected = Continuous(log(discount(c, from) / discount(c, to)) / (to - from))
            @test forward(c, from, to) ≈ expected rtol = 1e-9
        end
    end
end
