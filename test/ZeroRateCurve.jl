using ForwardDiff

@testset "ZeroRateCurve" begin
    rates = [0.02, 0.03, 0.035, 0.04]
    tenors = [1.0, 2.0, 5.0, 10.0]

    @testset "returns the concrete knot curve" begin
        @test !(ZeroRateCurve isa Type)
        for spline in (
                Spline.Linear(), Spline.Quadratic(), Spline.Cubic(),
                Spline.PCHIP(), Spline.Akima(), Spline.MonotoneConvex(), Spline.BSpline(1),
            )
            flat = ZeroRateCurve(fill(0.05, 4), tenors, spline)
            zrc = ZeroRateCurve(rates, tenors, spline)
            direct = spline isa Spline.MonotoneConvex ? Yield.MonotoneConvex(rates, tenors) :
                Yield.Spline(spline, tenors, rates)
            @test zrc isa Yield.AbstractInterpolatedZeroCurve
            @test zrc isa (spline isa Spline.MonotoneConvex ? Yield.MonotoneConvex : Yield.Spline)
            @test zrc == direct && typeof(zrc) == typeof(direct)
            for t in (0.0, 1.0e-20, 0.5, 3.0, 10.0, 2.0e4)
                @test rate(zero(flat, t)) ≈ 0.05
                @test zero(zrc, t) == zero(direct, t)
                @test rate(zero(flat + Yield.Constant(Continuous(0.01)), t)) ≈ 0.06
                @test rate(zero(flat + ((z, t) -> z + Continuous(0.01)), t)) ≈ 0.06
            end
            @test_throws DomainError zero(zrc, -1.0)
            @test_throws DomainError discount(zrc, -1.0)
        end
        @test_throws DomainError Yield.instantaneous_forward(ZeroRateCurve(rates, tenors), -1.0)
        # `Yield.Spline` covers the DataInterpolations methods only
        @test_throws ArgumentError Yield.Spline(Spline.MonotoneConvex(), tenors, rates)

        zrc = ZeroRateCurve(rates, tenors, Spline.Linear())
        for t in (0.0, 1.0e-20, 3.0, 2.0e4)
            g = ForwardDiff.gradient(rs -> rate(zero(ZeroRateCurve(rs, tenors, Spline.Linear()), t)), rates)
            @test g ≈ ForwardDiff.gradient(rs -> rate(zero(reconstruct(zrc; rates = rs), t)), rates)
        end
    end

    @testset "flat zero rate before the first knot" begin
        # DataInterpolations-backed curves hold the first knot's zero rate on [0, t₁]
        for spline in (
                    Spline.Linear(), Spline.Quadratic(), Spline.Cubic(),
                    Spline.PCHIP(), Spline.Akima(), Spline.BSpline(1), Spline.BSpline(3),
                ),
                extrapolation in (:flat_forward, :flat_zero, :linear, :extension)
            zrc = ZeroRateCurve(rates, tenors, spline; extrapolation)
            for t in (0.0, 0.1, 0.5, 0.999, 1.0)
                @test rate(zero(zrc, t)) == first(rates)
            end
            @test discount(zrc, 0.5) == exp(-first(rates) * 0.5)
        end
        # MonotoneConvex interpolates from t = 0 as part of the Hagan-West construction:
        # its zero rate at the origin is the instantaneous forward f(0), not z₁
        mc = ZeroRateCurve(rates, tenors)
        @test rate(zero(mc, 0.0)) == Yield.instantaneous_forward(mc, 0.0)
        @test rate(zero(mc, 0.0)) != first(rates)
        @test rate(zero(mc, 1.0)) ≈ first(rates)
    end

    @testset "MonotoneConvex (default)" begin
        zrc = ZeroRateCurve(rates, tenors)

        @testset "t=0 returns 1.0" begin
            @test discount(zrc, 0.0) == 1.0
        end

        @testset "exact tenor points" begin
            for (r, t) in zip(rates, tenors)
                @test discount(zrc, t) ≈ exp(-r * t)
            end
        end

        @testset "callable interface" begin
            @test zrc(1.0) ≈ exp(-0.02 * 1.0)
            @test zrc(0.0) == 1.0
        end

        @testset "inherited methods" begin
            # zero rate extraction
            z = zero(zrc, 2.0)
            @test FinanceCore.rate(z) ≈ 0.03 atol = 1.0e-10
        end
    end

    @testset "Linear" begin
        zrc = ZeroRateCurve(rates, tenors, Spline.Linear())

        @testset "exact tenor points" begin
            for (r, t) in zip(rates, tenors)
                @test discount(zrc, t) ≈ exp(-r * t)
            end
        end

        @testset "interpolation between tenors" begin
            # t=3.5 is between tenors 2.0 (r=0.03) and 5.0 (r=0.035)
            t = 3.5
            w = (3.5 - 2.0) / (5.0 - 2.0)
            r_interp = 0.03 + w * (0.035 - 0.03)
            @test discount(zrc, t) ≈ exp(-r_interp * t) atol = 1.0e-10
        end
    end

    @testset "from AbstractYieldModel" begin
        @testset "round-trip with Constant" begin
            c = Yield.Constant(0.05)
            tenors = [1.0, 2.0, 5.0, 10.0]
            zrc = ZeroRateCurve(c, tenors)
            for t in tenors
                @test discount(zrc, t) ≈ discount(c, t) atol = 1.0e-10
            end
        end

        @testset "round-trip with NelsonSiegel" begin
            ns = Yield.NelsonSiegel(1.0, 0.04, -0.02, 0.01)
            tenors = [1.0, 2.0, 5.0, 10.0, 20.0]
            zrc = ZeroRateCurve(ns, tenors)
            for t in tenors
                @test discount(zrc, t) ≈ discount(ns, t) atol = 1.0e-10
            end
        end

        @testset "explicit spline kwarg" begin
            c = Yield.Constant(0.04)
            tenors = [1.0, 5.0, 10.0]
            zrc = ZeroRateCurve(c, tenors; spline = Spline.Linear())
            for t in tenors
                @test discount(zrc, t) ≈ discount(c, t) atol = 1.0e-10
            end
        end

        @testset "unsorted tenors are sorted automatically" begin
            c = Yield.Constant(0.04)
            zrc_sorted = ZeroRateCurve(c, [1.0, 5.0, 10.0])
            zrc_unsorted = ZeroRateCurve(c, [10.0, 1.0, 5.0])
            for t in [1.0, 3.0, 5.0, 7.0, 10.0]
                @test discount(zrc_sorted, t) ≈ discount(zrc_unsorted, t) atol = 1.0e-10
            end
        end

        @testset "error on non-positive tenors" begin
            c = Yield.Constant(0.05)
            @test_throws ArgumentError ZeroRateCurve(c, [0.0, 1.0, 2.0])
            @test_throws ArgumentError ZeroRateCurve(c, [-1.0, 1.0, 2.0])
        end
    end

    @testset "Cubic" begin
        cubic_rates = [0.02, 0.025, 0.03, 0.035, 0.04]
        cubic_tenors = [1.0, 2.0, 5.0, 7.0, 10.0]
        zrc = ZeroRateCurve(cubic_rates, cubic_tenors, Spline.Cubic())

        @testset "exact tenor points" begin
            for (r, t) in zip(cubic_rates, cubic_tenors)
                @test discount(zrc, t) ≈ exp(-r * t) atol = 1.0e-8
            end
        end

        @testset "two-point case matches linear" begin
            r2 = [0.03, 0.05]
            t2 = [1.0, 5.0]
            zrc_lin = ZeroRateCurve(r2, t2, Spline.Linear())
            zrc_cub = ZeroRateCurve(r2, t2, Spline.Cubic())
            @test discount(zrc_lin, 3.0) ≈ discount(zrc_cub, 3.0) atol = 1.0e-6
        end
    end

    @testset "eager-build: ForwardDiff pass-through" begin
        # The eager build runs once with the input rate type; for Dual-typed
        # rates, the model itself becomes Dual-typed and propagates through
        # discount. At an exact knot t = tenors[k], discount = exp(-rates[k] * t),
        # so ∂discount/∂rates[k] = -t · discount.
        tenors_ad = [1.0, 2.0, 5.0, 10.0]
        rates_ad = [0.02, 0.03, 0.035, 0.04]
        for spl in (Spline.Linear(), Spline.MonotoneConvex())
            f(r) = discount(ZeroRateCurve([r, rates_ad[2:end]...], tenors_ad, spl), 1.0)
            ad = ForwardDiff.derivative(f, 0.02)
            @test ad ≈ -1.0 * exp(-0.02 * 1.0) atol = 1.0e-12
        end

        # The terminal forward (including the interpolant's endpoint derivative)
        # is computed with the rate's numeric type, so long-end extrapolation must
        # retain AD information rather than freezing the boundary at construction.
        # `:extension` deliberately delegates to DataInterpolations' legacy tail;
        # cubic, B-spline and Akima extension can themselves return NaN sensitivities
        # far beyond the grid, which is one reason it is no longer the default.
        for spl in (
                    Spline.Linear(), Spline.Quadratic(), Spline.Cubic(),
                    Spline.BSpline(3), Spline.PCHIP(), Spline.Akima(),
                ),
                extrapolation in (:flat_forward, :flat_zero, :linear)
            f(r) = rate(
                zero(
                    ZeroRateCurve(
                        [rates_ad[1:(end - 1)]...; r], tenors_ad, spl; extrapolation
                    ), 100.0
                )
            )
            @test isfinite(ForwardDiff.derivative(f, last(rates_ad)))
        end
    end

    @testset "eager-build: structural equality preserved" begin
        # Two ZRCs built from `==`-equal inputs must compare `==` despite their
        # internal interpolation caches being different prebuilt instances.
        r = [0.02, 0.03, 0.04]; t = [1.0, 2.0, 5.0]
        @test ZeroRateCurve(r, t, Spline.Linear()) == ZeroRateCurve(copy(r), copy(t), Spline.Linear())
        @test hash(ZeroRateCurve(r, t, Spline.Linear())) == hash(ZeroRateCurve(copy(r), copy(t), Spline.Linear()))
        # Different rates must compare unequal
        @test ZeroRateCurve(r, t, Spline.Linear()) != ZeroRateCurve(r .+ 0.01, t, Spline.Linear())
        # Different splines must compare unequal — relies on Sp.SplineCurve
        # subtypes implementing equality correctly (singleton splines under
        # `Sp.Linear()`, `Sp.Cubic()` etc. are `===` to other instances of the
        # same type).
        @test ZeroRateCurve(r, t, Spline.Linear()) != ZeroRateCurve(r, t, Spline.Cubic())
        @test hash(ZeroRateCurve(r, t, Spline.Linear())) != hash(ZeroRateCurve(r, t, Spline.Cubic()))
        # The extrapolation policy is value-carrying state too.
        @test ZeroRateCurve(r, t, Spline.Linear()) !=
            ZeroRateCurve(r, t, Spline.Linear(); extrapolation = :flat_zero)
    end

    # ─── Owned, read-only storage ─────────────────────────────────────────────

    @testset "read-only storage" begin
        for spl in (Spline.Linear(), Spline.MonotoneConvex())
            zrc = ZeroRateCurve([0.02, 0.03, 0.04], [1.0, 2.0, 5.0], spl)
            df1, h1 = discount(zrc, 1.0), hash(zrc)
            @test zrc.rates isa AbstractVector{Float64}
            @test zrc.tenors isa AbstractVector{Float64}
            # every mutation path must throw
            @test_throws ArgumentError zrc.rates[1] = 0.2
            @test_throws ArgumentError zrc.rates .= 0.0
            @test_throws ArgumentError fill!(zrc.tenors, 1.0)
            @test_throws ArgumentError reverse!(zrc.tenors)
            @test_throws ArgumentError sort!(zrc.rates; rev = true)
            @test_throws ArgumentError view(zrc.rates, 1:2)[1] = 0.2
            @test_throws ArgumentError view(zrc.rates, :) .= 0.0
            # hidden backing storage; derived caches are listed only as private properties
            @test_throws ArgumentError zrc.rates._data
            @test propertynames(zrc.rates) == ()
            @test propertynames(zrc.rates, true) == ()
            @test propertynames(zrc) == (:spline, :rates, :tenors, :extrapolation)
            @test propertynames(zrc, true) == fieldnames(typeof(zrc))
            @test all(p -> startswith(String(p), "_"), setdiff(propertynames(zrc, true), propertynames(zrc)))
            @test knot_rates(zrc) === zrc.rates && knot_tenors(zrc) === zrc.tenors
            @test_throws ArgumentError knot_rates(zrc)[1] = 0.2
            @test_throws ArgumentError knot_tenors(zrc)[1] = 0.2
            # nothing above changed the curve
            @test discount(zrc, 1.0) == df1
            @test hash(zrc) == h1
            # read paths behave like a Vector
            c = copy(zrc.rates)
            @test c isa Vector{Float64}
            c[1] = 0.9
            @test zrc.rates[1] == 0.02
            @test collect(zrc.tenors) isa Vector{Float64}
            @test hash(zrc.rates) == hash(collect(zrc.rates))
            @test isequal(zrc.rates, collect(zrc.rates))
            @test zrc.rates == [0.02, 0.03, 0.04]
            @test searchsortedlast(zrc.tenors, 3.0) == 2
            @test zrc.rates .+ 0.01 isa Vector{Float64}
            @test collect(zip(zrc.tenors, zrc.rates)) == [(1.0, 0.02), (2.0, 0.03), (5.0, 0.04)]
        end
    end

    @testset "caller-input isolation" begin
        r = [0.02, 0.03, 0.04]; t = [1.0, 2.0, 5.0]
        for spl in (Spline.Linear(), Spline.MonotoneConvex())
            zrc = ZeroRateCurve(r, t, spl)
            ref = ZeroRateCurve(copy(r), copy(t), spl)
            d = Dict(zrc => :found)
            r[2] = 0.1; t[2] = 3.0            # mutate the caller's inputs
            @test discount(zrc, 2.0) ≈ exp(-0.03 * 2.0)
            @test zrc == ref
            @test isequal(zrc, ref)
            @test hash(zrc) == hash(ref)
            @test d[zrc] == :found
            @test d[ref] == :found
            @test zrc != ZeroRateCurve(r, t, spl)
            r[2] = 0.03; t[2] = 2.0            # restore for the next spline
        end
    end

    @testset "signed zero: isequal/hash contract" begin
        a = ZeroRateCurve([-0.0], [1.0])
        a′ = ZeroRateCurve([-0.0], [1.0])
        b = ZeroRateCurve([0.0], [1.0])
        @test a == b                       # `==` follows array `==` (-0.0 == 0.0)
        @test !isequal(a, b)               # `isequal` follows array `isequal`
        @test isequal(a, a′) && hash(a) == hash(a′)
        d = Dict(a => 1)
        @test haskey(d, a′)
        @test !haskey(d, b)                # same behaviour as `Dict([-0.0] => 1)` with key `[0.0]`
    end

    # ─── Accessors / ConstructionBase ─────────────────────────────────────────

    @testset "ConstructionBase round-trips" begin
        CB = Accessors.ConstructionBase
        z = ZeroRateCurve([0.02, 0.03, 0.04], [1.0, 2.0, 5.0], Spline.Linear())
        @test keys(CB.getproperties(z)) == (:spline, :rates, :tenors, :extrapolation)
        @test CB.setproperties(z, CB.getproperties(z)) == z
        raw = CB.constructorof(typeof(z))(CB.getfields(z)...)
        @test raw == z && discount(raw, 3.5) == discount(z, 3.5)
        @test Accessors.mapproperties(identity, z) == z
        @test Accessors.getall(z, Accessors.Properties()) ==
            (z.spline, z.rates, z.tenors, z.extrapolation)
        newrates = [0.03, 0.04, 0.05]
        zp = Accessors.setall(
            z, Accessors.Properties(),
            (z.spline, newrates, z.tenors, :flat_zero)
        )
        @test zp == ZeroRateCurve(
            newrates, [1.0, 2.0, 5.0], Spline.Linear();
            extrapolation = :flat_zero
        )
        @test discount(zp, 5.0) ≈ exp(-0.05 * 5.0)
        mc = ZeroRateCurve([0.02, 0.03, 0.04], [1.0, 2.0, 5.0])
        rawmc = CB.constructorof(typeof(mc))(CB.getfields(mc)...)
        @test rawmc == mc && discount(rawmc, 3.5) == discount(mc, 3.5)
    end

    @testset "Accessors rebuild the cached model" begin
        CB = Accessors.ConstructionBase
        r = [0.02, 0.03, 0.04]; t = [1.0, 2.0, 5.0]
        for spl in (Spline.Linear(), Spline.MonotoneConvex())
            zrc = ZeroRateCurve(r, t, spl)
            new = Accessors.@set zrc.rates[2] = 0.05
            @test new.rates[2] == 0.05
            @test discount(new, t[2]) ≈ exp(-0.05 * t[2])        # cache rebuilt
            @test new == ZeroRateCurve([0.02, 0.05, 0.04], t, spl)
            @test discount(zrc, t[2]) ≈ exp(-0.03 * t[2])        # original untouched
            @test zrc.rates[2] == 0.03
        end
        zrc = ZeroRateCurve(r, t, Spline.MonotoneConvex())
        # spline swap rebuilds the interpolant
        lin = Accessors.@set zrc.spline = Spline.Linear()
        @test lin == ZeroRateCurve(r, t, Spline.Linear())
        @test discount(lin, 3.5) ≈ discount(ZeroRateCurve(r, t, Spline.Linear()), 3.5)
        @test discount(lin, 3.5) != discount(zrc, 3.5)
        flat_zero = Accessors.@set lin.extrapolation = :flat_zero
        @test flat_zero.extrapolation === :flat_zero
        @test rate(zero(flat_zero, 10.0)) ≈ last(r)
        @test_throws ArgumentError Accessors.@set lin.extrapolation = :unknown
        # tenor patches: whole-grid replacement and order-preserving single-element updates work
        g = Accessors.@set zrc.tenors = [1.0, 3.0, 6.0]
        @test g == ZeroRateCurve(r, [1.0, 3.0, 6.0], Spline.MonotoneConvex())
        @test discount(g, 6.0) ≈ exp(-0.04 * 6.0)
        g2 = Accessors.@set zrc.tenors[2] = 1.5
        @test g2.tenors == [1.0, 1.5, 5.0] && discount(g2, 1.5) ≈ exp(-0.03 * 1.5)
        # ... but anything that breaks the invariants throws
        @test_throws ArgumentError Accessors.@set zrc.tenors[2] = 0.5     # breaks ordering
        @test_throws ArgumentError Accessors.@set zrc.tenors[2] = 5.0     # duplicate
        @test_throws ArgumentError Accessors.@set zrc.tenors[1] = -1.0    # negative
        @test_throws ArgumentError Accessors.@set zrc.rates = [NaN, 0.03, 0.04]
        @test_throws ArgumentError Accessors.@set zrc.rates = [0.02, 0.03]  # length mismatch
        # the derived caches cannot be patched
        @test_throws ArgumentError (Accessors.@set zrc._f = Float64[])
        @test_throws ArgumentError CB.setproperties(zrc, (_tail = nothing,))
        @test_throws ArgumentError CB.setproperties(lin, (_fn = nothing,))
        @test_throws ArgumentError CB.setproperties(zrc, (foo = 1,))
        # ConstructionBase reconstruction preserves the public policy and rebuilds the cache.
        sp = CB.setproperties(zrc, (rates = [0.03, 0.04, 0.05],))
        @test sp == ZeroRateCurve([0.03, 0.04, 0.05], t, Spline.MonotoneConvex())
        @test discount(sp, 5.0) ≈ exp(-0.05 * 5.0)
        # the cache is retyped when rates become Duals via Accessors
        dz = Accessors.@set zrc.rates[1] = ForwardDiff.Dual(0.02, 1.0)
        @test eltype(dz.rates) <: ForwardDiff.Dual
        @test discount(dz, 1.0) isa ForwardDiff.Dual
        @test ForwardDiff.partials(discount(dz, 1.0), 1) ≈ -1.0 * exp(-0.02 * 1.0) atol = 1.0e-12
        # The policy is keyword-only and there is no public positional cache constructor.
        @test_throws MethodError ZeroRateCurve(r, t, Spline.Linear(), :flat_zero)
        @test_throws MethodError ZeroRateCurve(r, t, Spline.Linear(), getfield(lin, :_fn))
    end

    # ─── Validation ───────────────────────────────────────────────────────────

    @testset "validation" begin
        @test_throws ArgumentError ZeroRateCurve(
            [0.02, 0.03], [1.0, 2.0];
            extrapolation = :unknown
        )
        @test_throws ArgumentError ZeroRateCurve(
            [0.02, 0.03], [1.0, 2.0];
            extrapolation = "flat_forward"
        )
        @test_throws ArgumentError ZeroRateCurve(
            [0.02, 0.03], [1.0, 2.0];
            extrapolation = :extension
        )  # MonotoneConvex has no polynomial extension
        @test_throws ArgumentError ZeroRateCurve([0.02, 0.03], [1.0])                  # length mismatch
        @test_throws ArgumentError ZeroRateCurve(Float64[], Float64[])                 # empty
        @test_throws ArgumentError ZeroRateCurve([0.02, 0.03], [-1.0, 2.0])            # negative tenor
        @test_throws ArgumentError ZeroRateCurve([0.02, 0.03], [NaN, 2.0])             # NaN tenor
        @test_throws ArgumentError ZeroRateCurve([0.02, 0.03], [1.0, Inf])             # Inf tenor
        @test_throws ArgumentError ZeroRateCurve([NaN, 0.03], [1.0, 2.0])              # NaN rate
        @test_throws ArgumentError ZeroRateCurve([Inf, 0.03], [1.0, 2.0])              # Inf rate
        @test_throws ArgumentError ZeroRateCurve([0.02, -Inf], [1.0, 2.0])             # -Inf rate
        @test_throws ArgumentError ZeroRateCurve(["a"], [1.0])                         # non-numeric
        @test_throws ArgumentError ZeroRateCurve([0.02], ["1"])                        # non-numeric
        # duplicate / unsorted tenors are rejected for every interpolant (Linear/PCHIP/Akima
        # used to accept duplicates silently and produce NaN; MonotoneConvex accepted unsorted)
        for spl in (Spline.Linear(), Spline.Cubic(), Spline.PCHIP(), Spline.Akima(), Spline.MonotoneConvex())
            @test_throws ArgumentError ZeroRateCurve([0.02, 0.03, 0.04], [1.0, 1.0, 2.0], spl)
            @test_throws ArgumentError ZeroRateCurve([0.02, 0.03, 0.04], [2.0, 1.0, 3.0], spl)
        end
        # negative rates are valid
        zneg = ZeroRateCurve([-0.005, 0.01], [1.0, 5.0], Spline.Linear())
        @test discount(zneg, 1.0) ≈ exp(0.005)
        # the sampling form takes its spline as a keyword, not positionally
        @test_throws MethodError ZeroRateCurve(Yield.Constant(0.03), [1.0, 2.0], Spline.Linear())
        # invalid spline descriptors are rejected at descriptor construction
        @test_throws ArgumentError Spline.BSpline(-1)
        @test_throws ArgumentError Spline.BSpline(0)
        @test_throws ArgumentError Spline.PolynomialSpline(0)
        # orders above 3 used to build a cubic spline silently
        @test_throws ArgumentError Spline.PolynomialSpline(4)
        @test_throws ArgumentError Spline.PolynomialSpline(99)
        @test Spline.BSpline(3).order == 3 && Spline.Cubic().order == 3
        @test Spline.PolynomialSpline(1) == Spline.Linear()
        @test Spline.PolynomialSpline(3) == Spline.Cubic()
    end

    @testset "element-type promotion" begin
        z = ZeroRateCurve(Real[0.02f0, 0.03], [1, 2], Spline.Linear())
        @test eltype(z.rates) === Float64 && eltype(z.tenors) === Float64
        @test z == ZeroRateCurve([Float64(0.02f0), 0.03], [1.0, 2.0], Spline.Linear())
        zb = ZeroRateCurve((0.02, 0.03), (1.0f0, big"2"), Spline.Linear())
        @test eltype(zb.tenors) === BigFloat && eltype(zb.rates) === Float64
        @test discount(zb, big"1.5") isa BigFloat
        zd = ZeroRateCurve(Real[0.02, ForwardDiff.Dual(0.03, 1.0)], [1.0, 2.0], Spline.Linear())
        @test isconcretetype(eltype(zd.rates)) && eltype(zd.rates) <: ForwardDiff.Dual
        @test isfinite(ForwardDiff.value(discount(zd, 1.5)))
        zr = ZeroRateCurve((1:3) ./ 100, 1.0:3.0)          # ranges
        @test zr.rates == [0.01, 0.02, 0.03] && zr.tenors == [1.0, 2.0, 3.0]
        @test zr.tenors isa AbstractVector{Float64}
        zi = ZeroRateCurve([0.02, 0.03, 0.04], [1, 2, 5])  # Int tenors
        @test eltype(zi.tenors) === Float64
        @test zi == ZeroRateCurve([0.02, 0.03, 0.04], [1.0, 2.0, 5.0])
        # constructing from another curve's read-only fields copies
        zc = ZeroRateCurve(zi.rates, zi.tenors, Spline.Linear())
        @test zc.rates == zi.rates && zc.rates !== zi.rates
    end

    @testset "minimum knots per interpolant" begin
        expected = Dict(
            Spline.Linear() => 1, Spline.Quadratic() => 1, Spline.Cubic() => 1, Spline.BSpline(3) => 1,
            Spline.PCHIP() => 3, Spline.Akima() => 3, Spline.MonotoneConvex() => 1,
        )
        for (spl, k) in expected
            @test FinanceModels.Yield.__min_knots(spl) == k
            for T in (Float64, BigFloat)
                rates = T[0.02 + 0.005 * i for i in 0:(k - 1)]
                tenors = T[1 + 2i for i in 0:(k - 1)]
                if k > 1
                    err = try
                        ZeroRateCurve(rates[1:(k - 1)], tenors[1:(k - 1)], spl); nothing
                    catch e
                        e
                    end
                    @test err isa ArgumentError
                    @test occursin("requires at least $k knots", sprint(showerror, err))
                end
                z = ZeroRateCurve(rates, tenors, spl)
                for τ in (T(0.5), (first(tenors) + last(tenors)) / 2, last(tenors) + 1)
                    df = discount(z, τ)
                    @test isfinite(df) && 0 < df <= 1
                end
            end
        end
        # BigFloat Akima specifically (2 knots used to hit an UndefRefError deep in DataInterpolations)
        @test_throws ArgumentError ZeroRateCurve(BigFloat[0.02, 0.03], BigFloat[1, 2], Spline.Akima())
        za = ZeroRateCurve(BigFloat[0.02, 0.03, 0.035], BigFloat[1, 2, 5], Spline.Akima())
        @test eltype(za.rates) === BigFloat
        @test isfinite(discount(za, big"1.5"))
        # a single knot is a flat curve for every method that accepts one, under every policy
        for spl in (Spline.MonotoneConvex(), Spline.Linear(), Spline.Cubic(), Spline.BSpline(3)),
                extrapolation in (:flat_forward, :flat_zero, :linear)
            z1 = ZeroRateCurve([0.03], [2.0], spl; extrapolation)
            for t in (0.0, 1.0, 2.0, 4.0, 50.0)
                @test rate(zero(z1, t)) ≈ 0.03
            end
        end
        z1 = ZeroRateCurve([0.03], [2.0], Spline.Cubic(); extrapolation = :extension)
        @test rate(zero(z1, 50.0)) == 0.03
    end

    @testset "zero-tenor knot, interior times" begin
        for spl in (Spline.Linear(), Spline.MonotoneConvex())
            z0 = ZeroRateCurve([0.02, 0.03], [0.0, 1.0], spl)
            @test discount(z0, 0.0) == 1.0
            df = discount(z0, 0.5)
            @test isfinite(df) && 0 < df < 1
            @test discount(z0, 1.0) ≈ exp(-0.03)
            @test isfinite(discount(z0, 2.0))
        end
        @test discount(ZeroRateCurve([0.02, 0.03], [0.0, 1.0], Spline.Linear()), 0.5) ≈ exp(-0.025 * 0.5)
    end

    @testset "sampling form" begin
        c = Yield.Constant(0.05)
        # a stateful iterator is normalised exactly once
        zs = ZeroRateCurve(c, Iterators.Stateful([1.0, 2.0]))
        @test zs == ZeroRateCurve(c, [1.0, 2.0])
        @test zs.tenors == [1.0, 2.0]
        # sampling goes through `zero`, which is stable at extreme tenors
        # (`-log(discount)/t` gave -0.0 at 1e-20 and Inf at 2e4)
        ze = ZeroRateCurve(c, [1.0e-20, 1.0, 2.0e4])
        @test all(r -> r ≈ log(1.05), ze.rates)
        @test_throws ArgumentError ZeroRateCurve(c, [1.0, 1.0])
        @test_throws ArgumentError ZeroRateCurve(c, Float64[])
        @test_throws ArgumentError ZeroRateCurve(c, [0.0, 1.0])
        @test_throws ArgumentError ZeroRateCurve(c, [-1.0, 1.0])
        # non-finite tenors are rejected before the source curve is ever evaluated
        touched = Float64[]
        spy = c + (z, t) -> (push!(touched, t); z)
        @test_throws ArgumentError ZeroRateCurve(spy, [1.0, Inf])
        @test isempty(touched)
        @test_throws ArgumentError ZeroRateCurve(spy, [NaN, 1.0])
        @test isempty(touched)
    end

    @testset "generic fit via __default_optic" begin
        t = [1.0, 2.0, 5.0, 10.0]
        target = [0.02, 0.025, 0.03, 0.035]
        qs = ZCBYield.(Continuous.(target), t)        # continuous quotes: fitted zero rates == target
        zrc0 = ZeroRateCurve(fill(0.01, length(t)), t)
        # one batch optic over all knot rates: a single rebuild per optimizer candidate
        optics = FinanceModels.__default_optic(zrc0)
        @test length(optics) == 1
        o = first(optics).first
        @test o isa FinanceModels.KnotRatesOptic
        @test Accessors.getall(zrc0, o) === Tuple(zrc0.rates)
        zb = Accessors.setall(zrc0, o, target)
        @test zb isa Yield.MonotoneConvex && zb.rates == target && zb.tenors == t
        @test discount(zb, 2.0) ≈ exp(-0.025 * 2.0)
        @test zrc0.rates == fill(0.01, length(t))                  # original untouched
        @test Accessors.modify(x -> 2x, zrc0, o).rates == fill(0.02, length(t))
        @test_throws ArgumentError Accessors.setall(zrc0, o, [NaN, 0.0, 0.0, 0.0])   # still validated
        fitted = fit(zrc0, qs)
        @test fitted isa Yield.MonotoneConvex
        @test fitted.tenors == t
        @test maximum(abs, present_value(fitted, q.instrument) - q.price for q in qs) < 1.0e-6
        fl = fit(ZeroRateCurve(fill(0.01, length(t)), t, Spline.Linear()), qs)
        @test fl.rates ≈ target atol = 1.0e-5
    end

    # ─── Knot-curve interface ─────────────────────────────────────────────────

    @testset "reconstruct" begin
        local rates = [0.02, 0.03, 0.035, 0.04]
        local tenors = [1.0, 2.0, 5.0, 10.0]
        splines = (
            Spline.Linear(), Spline.Quadratic(), Spline.Cubic(), Spline.PCHIP(),
            Spline.Akima(), Spline.BSpline(3), Spline.MonotoneConvex(),
        )
        ts = (0.0, 0.3, 1.0, 2.7, 5.0, 9.9, 10.0, 25.0, 1.0e3)
        for spline in splines, extrapolation in (:flat_forward, :flat_zero, Yield.FlatForwardAt(Continuous(0.03)))
            c = ZeroRateCurve(rates, tenors, spline; extrapolation)
            # the same inputs rebuild an equal curve that prices bitwise identically
            for r in (reconstruct(c), reconstruct(c; rates = collect(knot_rates(c))))
                @test r == c && isequal(r, c) && hash(r) == hash(c) && typeof(r) == typeof(c)
                @test all(discount(r, t) === discount(c, t) for t in ts)
            end
            # each argument replaces only its own input
            up = reconstruct(c; rates = rates .+ 0.01)
            @test knot_rates(up) == rates .+ 0.01 && knot_tenors(up) == tenors
            @test up.spline == spline && up.extrapolation == extrapolation
            @test discount(up, 5.0) ≈ exp(-(0.035 + 0.01) * 5.0)
            @test knot_tenors(reconstruct(c; tenors = tenors .* 2)) == tenors .* 2
            @test reconstruct(c; extrapolation = :linear).extrapolation === :linear
            @test knot_rates(c) == rates                                # original untouched
        end
        # replacing the method can change the concrete type
        mc = ZeroRateCurve(rates, tenors)
        lin = reconstruct(mc; spline = Spline.Linear())
        @test lin isa Yield.Spline && lin == ZeroRateCurve(rates, tenors, Spline.Linear())
        @test reconstruct(lin; spline = Spline.MonotoneConvex()) == mc
        # the requested method survives a short grid (Cubic on two knots evaluates linearly)
        short = ZeroRateCurve(rates[1:2], tenors[1:2], Spline.Cubic())
        @test short.spline == Spline.Cubic()
        @test reconstruct(short; rates, tenors) == ZeroRateCurve(rates, tenors, Spline.Cubic())
        # validation is the constructor's
        @test_throws ArgumentError reconstruct(mc; rates = rates[1:3])
        @test_throws ArgumentError reconstruct(mc; tenors = reverse(tenors))
        @test_throws ArgumentError reconstruct(mc; extrapolation = :extension)
        @test_throws ArgumentError reconstruct(mc; rates = [NaN, 0.03, 0.035, 0.04])
        # knot-rate derivatives through `reconstruct` are exact at the knots
        for spline in splines
            c = ZeroRateCurve(rates, tenors, spline)
            g = ForwardDiff.gradient(z -> discount(reconstruct(c; rates = z), 5.0), collect(knot_rates(c)))
            @test g ≈ [0.0, 0.0, -5.0 * exp(-0.035 * 5.0), 0.0] atol = 1.0e-12
        end
    end

    @testset "knot vectors are read-only and copyable" begin
        local rates = [0.02, 0.03, 0.035, 0.04]
        local tenors = [1.0, 2.0, 5.0, 10.0]
        for spline in (Spline.Linear(), Spline.MonotoneConvex())
            c = ZeroRateCurve(rates, tenors, spline)
            @test_throws ArgumentError knot_rates(c) .= 0.0
            @test_throws ArgumentError sort!(knot_tenors(c); rev = true)
            v = copy(knot_rates(c))
            @test v isa Vector{Float64}
            v[1] = 0.5
            @test knot_rates(c)[1] == 0.02
        end
    end

    @testset "structural equality across concrete types" begin
        local rates = [0.02, 0.03, 0.035, 0.04]
        local tenors = [1.0, 2.0, 5.0, 10.0]
        lin = ZeroRateCurve(rates, tenors, Spline.Linear())
        mc = ZeroRateCurve(rates, tenors)
        @test lin != mc && !isequal(lin, mc) && hash(lin) != hash(mc)
        @test lin != ZeroRateCurve(rates, tenors, Spline.BSpline(1))     # numerically identical, different method
        @test lin == Yield.Spline(Spline.Linear(), tenors, rates)
        @test mc == Yield.MonotoneConvex(rates, tenors)
        d = Dict(lin => 1, mc => 2)
        @test d[ZeroRateCurve(copy(rates), copy(tenors), Spline.Linear())] == 1
        @test d[ZeroRateCurve(copy(rates), copy(tenors))] == 2
    end

    @testset "show prints a construction call that rebuilds the curve" begin
        local rates = [0.02, 0.03, 0.035, 0.04]
        local tenors = [1.0, 2.0, 5.0, 10.0]
        for c in (
                ZeroRateCurve(rates, tenors),
                ZeroRateCurve(rates, tenors, Spline.Cubic(); extrapolation = :linear),
                ZeroRateCurve(rates, tenors, Spline.BSpline(3); extrapolation = :extension),
                ZeroRateCurve(rates, tenors, Spline.PCHIP(); extrapolation = Yield.FlatForwardAt(Periodic(0.04, 2))),
            )
            s = repr(c)
            @test startswith(s, "ZeroRateCurve(")
            rebuilt = Core.eval(Main, Meta.parse(s))
            @test rebuilt == c && typeof(rebuilt) == typeof(c)
        end
    end

    @testset "one monotone convex selector (#272)" begin
        # `Spline.MonotoneConvex()` is the only selector: construction and every fit route
        # through it to the native curve; the former `Yield.MonotoneConvex()` placeholder is gone.
        @test_throws MethodError Yield.MonotoneConvex()
        qs = ZCBYield.([0.02, 0.025, 0.03], [1.0, 2.0, 5.0])
        fitted = fit(Spline.MonotoneConvex(), qs)
        @test fitted isa Yield.MonotoneConvex
        @test fit(Spline.MonotoneConvex(), qs, Fit.Loss(x -> x^2)) == fitted
        @test maximum(abs(present_value(fitted, q.instrument) - q.price) for q in qs) < 1.0e-7
        @test_throws ArgumentError fit(Spline.MonotoneConvex(), qs; extrapolation = :extension)
        @test_throws "fit(Spline.MonotoneConvex(), quotes)" fit(Spline.MonotoneConvex(), qs, Fit.Bootstrap())
        # refitting an existing curve varies only its knot rates
        refit = fit(ZeroRateCurve(fill(0.01, 3), [1.0, 2.0, 5.0]), qs)
        @test refit isa Yield.MonotoneConvex && knot_tenors(refit) == [1.0, 2.0, 5.0]
        @test knot_rates(refit) ≈ knot_rates(fitted) atol = 1.0e-5
    end

    @testset "bootstrap knots are the quote maturities" begin
        qs = CMTYield.([0.03, 0.032, 0.035, 0.037], [0.5, 1.0, 2.0, 5.0])
        for spline in (Spline.Linear(), Spline.BSpline(1))
            c = fit(spline, qs, Fit.Bootstrap())
            @test c isa Yield.Spline && c.spline == spline
            @test knot_tenors(c) == [0.5, 1.0, 2.0, 5.0]
            @test length(knot_rates(c)) == 4
            # the knots alone determine the curve, including its flat short end
            @test reconstruct(c; rates = collect(knot_rates(c))) == c
            @test rate(zero(c, 0.0)) == first(knot_rates(c)) == rate(zero(c, 0.5))
            for q in qs
                @test present_value(c, q.instrument) ≈ q.price atol = 1.0e-12
            end
        end
        one_quote = fit(Spline.Linear(), [ZCBPrice(0.95, 2.0)], Fit.Bootstrap())
        @test knot_tenors(one_quote) == [2.0]
        @test discount(one_quote, 2.0) ≈ 0.95
        @test rate(zero(one_quote, 0.0)) == rate(zero(one_quote, 2.0)) == rate(zero(one_quote, 10.0))
    end
end
