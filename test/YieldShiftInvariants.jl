@testset "Yield shift internal-consistency invariants" begin
    base_const = Yield.Constant(Continuous(0.05))
    zqs = ZCBYield.([0.01, 0.018, 0.024, 0.029, 0.032], [1.0, 3.0, 5.0, 10.0, 20.0])
    base_spline = fit(Spline.Linear(), zqs, Fit.Bootstrap())
    bases = (("Constant", base_const), ("Spline", base_spline))
    tenors = [1.0, 5.0, 10.0, 20.0]

    @testset "compounding-convention equivalence ($name base)" for (name, base) in bases
        @testset "Continuous(c) ≡ Periodic(exp(c)-1, 1) ≡ Periodic(m=4)" begin
            c = 0.01
            s_cont = base + ((z, t) -> z + Continuous(c))
            s_per1 = base + ((z, t) -> z + Periodic(exp(c) - 1, 1))
            s_per4 = base + ((z, t) -> z + Periodic(4 * (exp(c / 4) - 1), 4))
            for t in tenors
                @test discount(s_cont, t) ≈ discount(s_per1, t)
                @test discount(s_cont, t) ≈ discount(s_per4, t)
                @test zero(s_cont, t) ≈ zero(s_per1, t)
                @test zero(s_cont, t) ≈ zero(s_per4, t)
            end
        end

        @testset "Periodic(r, m) ≡ Continuous(m * log(1 + r/m))" begin
            r, m = 0.04, 4
            c = m * log(1 + r / m)
            s_per = base + ((z, t) -> z + Periodic(r, m))
            s_cont = base + ((z, t) -> z + Continuous(c))
            for t in tenors
                @test discount(s_per, t) ≈ discount(s_cont, t)
                @test zero(s_per, t) ≈ zero(s_cont, t)
            end
        end
    end

    @testset "zero-magnitude shift = base ($name base)" for (name, base) in bases
        shift = base + ((z, t) -> z + Continuous(0.0))
        for t in tenors
            @test discount(shift, t) ≈ discount(base, t)
            @test zero(shift, t) ≈ zero(base, t)
        end
    end

    @testset "inverse symmetry ($name base, $sname shift)" for (name, base) in bases,
                                                                (sname, Δ) in (("Continuous", Continuous(0.0075)),
                                                                                ("Periodic(m=2)", Periodic(0.01, 2)))

        up = base + ((z, t) -> z + Δ)
        round_trip = up + ((z, t) -> z - Δ)
        for t in tenors
            @test discount(round_trip, t) ≈ discount(base, t)
            @test zero(round_trip, t) ≈ zero(base, t)
        end
    end

    @testset "additivity of stacked parallel shifts ($name base)" for (name, base) in bases
        a = Continuous(0.005)
        b = Continuous(0.003)
        stacked = (base + ((z, t) -> z + a)) + ((z, t) -> z + b)
        combined = base + ((z, t) -> z + a + b)
        for t in tenors
            @test discount(stacked, t) ≈ discount(combined, t)
            @test zero(stacked, t) ≈ zero(combined, t)
        end
    end

    @testset "stacked-shift order independence ($name base)" for (name, base) in bases
        a = Continuous(0.005)
        b = Continuous(0.003)
        ab = (base + ((z, t) -> z + a)) + ((z, t) -> z + b)
        ba = (base + ((z, t) -> z + b)) + ((z, t) -> z + a)
        for t in tenors
            @test discount(ab, t) ≈ discount(ba, t)
        end
    end

    @testset "discount/zero roundtrip with tenor-dependent rule ($name base)" for (name, base) in bases
        shift = base + ((z, t) -> z + Continuous(0.01 * t / 10))
        for t in tenors
            z = zero(shift, t)
            @test discount(shift, t) ≈ exp(-z.continuous_value * t)
        end
    end

    @testset "ProjectedShift τ=0 / fully-phased boundaries" begin
        base = Yield.Constant(Continuous(0.05))
        phase_in = (τ, z, t) -> z + Continuous(-0.015 * min(τ, 10) / 10)
        ps_at_0 = Yield.ProjectedShift(base, phase_in, 0.0)
        ps_full = Yield.ProjectedShift(base, phase_in, 10.0)
        ps_past = Yield.ProjectedShift(base, phase_in, 15.0)
        ts_full = Yield.TenorShift(base, (z, t) -> z + Continuous(-0.015))
        for t in tenors
            @test discount(ps_at_0, t) ≈ discount(base, t)
            @test discount(ps_full, t) ≈ discount(ts_full, t)
            @test discount(ps_past, t) ≈ discount(ps_full, t)
        end
    end

    @testset "ProjectedShift inverse symmetry across τ" begin
        base = Yield.Constant(Continuous(0.05))
        rule_up = (τ, z, t) -> z + Continuous(0.001 * τ)
        rule_dn = (τ, z, t) -> z - Continuous(0.001 * τ)
        for τ in [0.0, 1.0, 5.0, 10.0]
            up = Yield.ProjectedShift(base, rule_up, τ)
            rt = Yield.ProjectedShift(up, rule_dn, τ)
            for t in tenors
                @test discount(rt, t) ≈ discount(base, t)
                @test zero(rt, t) ≈ zero(base, t)
            end
        end
    end
end
