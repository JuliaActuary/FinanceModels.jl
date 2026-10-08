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
        # `Yield.Spline` covers the DataInterpolations methods only, and the direct form takes the
        # method positionally
        @test_throws MethodError Yield.Spline(Spline.MonotoneConvex(), tenors, rates)
        @test_throws MethodError ZeroRateCurve(rates, tenors; spline = Spline.Linear())

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

        @testset "tenors at and below zero" begin
            c = Yield.Constant(0.05)
            # a knot at t = 0 takes the source curve's zero-rate limit there
            z0 = ZeroRateCurve(c, [0.0, 1.0, 2.0])
            @test knot_tenors(z0) == [0.0, 1.0, 2.0] && all(r -> r ≈ log(1.05), knot_rates(z0))
            # a source without its own zero rate also has the limit, its short rate (the generic
            # zero rate L/t was 0/0 there, which the knot grid rejected)
            sw = Yield.SmithWilson(ufr = 0.03, α = 0.1)
            @test knot_rates(ZeroRateCurve(sw, [0.0, 1.0])) ≈ [0.03, 0.03] rtol = 1.0e-14
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
            # every mutation path throws: the vectors have no setindex!
            @test_throws Base.CanonicalIndexError zrc.rates[1] = 0.2
            @test_throws Base.CanonicalIndexError zrc.rates .= 0.0
            @test_throws Base.CanonicalIndexError fill!(zrc.tenors, 1.0)
            @test_throws Base.CanonicalIndexError reverse!(zrc.tenors)
            @test_throws Base.CanonicalIndexError sort!(zrc.rates; rev = true)
            @test_throws Base.CanonicalIndexError view(zrc.rates, 1:2)[1] = 0.2
            @test_throws Base.CanonicalIndexError view(zrc.rates, :) .= 0.0
            # derived caches are listed only as private properties
            @test propertynames(zrc) == (:spline, :rates, :tenors, :extrapolation)
            @test propertynames(zrc, true) == fieldnames(typeof(zrc))
            @test all(p -> startswith(String(p), "_"), setdiff(propertynames(zrc, true), propertynames(zrc)))
            @test knot_rates(zrc) === zrc.rates && knot_tenors(zrc) === zrc.tenors
            @test_throws Base.CanonicalIndexError knot_rates(zrc)[1] = 0.2
            @test_throws Base.CanonicalIndexError knot_tenors(zrc)[1] = 0.2
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
        # the derived caches cannot be patched: they are not `reconstruct` keywords
        @test_throws MethodError (Accessors.@set zrc._f = Float64[])
        @test_throws MethodError CB.setproperties(zrc, (_tail = nothing,))
        @test_throws MethodError CB.setproperties(lin, (_interp = nothing,))
        @test_throws MethodError CB.setproperties(zrc, (foo = 1,))
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
        @test_throws MethodError ZeroRateCurve(r, t, Spline.Linear(), getfield(lin, :_interp))
    end

    # ─── Validation ───────────────────────────────────────────────────────────

    @testset "validation" begin
        @test_throws ArgumentError ZeroRateCurve(
            [0.02, 0.03], [1.0, 2.0];
            extrapolation = :unknown
        )
        # a policy that is not a Symbol or FlatForwardAt has no tail method
        @test_throws MethodError ZeroRateCurve(
            [0.02, 0.03], [1.0, 2.0];
            extrapolation = "flat_forward"
        )
        @test_throws "not Spline.MonotoneConvex()" ZeroRateCurve(
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
        @test_throws MethodError ZeroRateCurve(["a"], [1.0])                           # non-numeric
        @test_throws MethodError ZeroRateCurve([0.02], ["1"])                          # non-numeric
        # duplicate / unsorted tenors are rejected for every interpolant (Linear/PCHIP/Akima
        # used to accept duplicates silently and produce NaN; MonotoneConvex accepted unsorted)
        for spl in (Spline.Linear(), Spline.Cubic(), Spline.PCHIP(), Spline.Akima(), Spline.MonotoneConvex())
            @test_throws ArgumentError ZeroRateCurve([0.02, 0.03, 0.04], [1.0, 1.0, 2.0], spl)
            @test_throws ArgumentError ZeroRateCurve([0.02, 0.03, 0.04], [2.0, 1.0, 3.0], spl)
        end
        # equal tenors are duplicates even as dual numbers with different partials, which
        # ForwardDiff orders by their partials
        dual_tenors = [ForwardDiff.Dual(1.0, 0.0), ForwardDiff.Dual(1.0, 1.0), ForwardDiff.Dual(2.0, 0.0)]
        for spl in (Spline.Linear(), Spline.Cubic(), Spline.MonotoneConvex())
            @test_throws "strictly increasing" ZeroRateCurve([0.02, 0.03, 0.04], dual_tenors, spl)
        end
        @test_throws "strictly increasing" ZeroRateCurve(Yield.Constant(0.05), dual_tenors)
        # the first tenor is also compared by its primal value: a zero tenor carrying a negative
        # partial is a valid knot at t = 0, not a negative tenor
        t0 = [ForwardDiff.Dual(0.0, -1.0), ForwardDiff.Dual(1.0, 0.0), ForwardDiff.Dual(2.0, 0.0)]
        zt0 = ZeroRateCurve([0.02, 0.03, 0.04], t0, Spline.Linear())
        @test ForwardDiff.value(discount(zt0, 1.5)) ≈ discount(ZeroRateCurve([0.02, 0.03, 0.04], [0.0, 1.0, 2.0], Spline.Linear()), 1.5) rtol = 1.0e-15
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
        @test_throws "strictly increasing" ZeroRateCurve(c, [1.0, 1.0])
        @test_throws "at least 1 knots" ZeroRateCurve(c, Float64[])
        @test_throws "≥ 0" ZeroRateCurve(c, [-1.0, 1.0])
        # the sampled grid is validated like any other: non-finite tenors are rejected
        @test_throws "tenors must be finite" ZeroRateCurve(c, [1.0, Inf])
        @test_throws "tenors must be finite" ZeroRateCurve(c, [NaN, 1.0])
    end

    @testset "refitting a knot curve" begin
        t = [1.0, 2.0, 5.0, 10.0]
        target = [0.02, 0.025, 0.03, 0.035]
        qs = ZCBYield.(Continuous.(target), t)        # continuous quotes: fitted zero rates == target
        # The refit keeps the curve's tenors, method and extrapolation policy, and is the spline
        # fit on those knots (same solve, same starting rates), whatever the curve's own rates:
        # flat PCHIP and Akima curves sit on a kink, where the earlier optic-based refit threw.
        for spline in (Spline.MonotoneConvex(), Spline.Linear(), Spline.Cubic(), Spline.PCHIP(), Spline.Akima())
            zrc0 = ZeroRateCurve(fill(0.01, length(t)), t, spline; extrapolation = :flat_zero)
            fitted = fit(zrc0, qs)
            @test fitted isa typeof(ZeroRateCurve(target, t, spline))
            @test fitted.tenors == t && fitted.spline == spline && fitted.extrapolation === :flat_zero
            @test maximum(abs, present_value(fitted, q.instrument) - q.price for q in qs) < 1.0e-6
            @test isequal(fitted, fit(spline, qs; extrapolation = :flat_zero))
            @test zrc0.rates == fill(0.01, length(t))                  # original untouched
        end
        # The curve's own rates are never the start, whether flat, tied, kinked or numerically flat
        # under the price loss: a 100% rate discounts 30 years to 1e-13, where the loss gradient is
        # below the solver's tolerance and a solve starting there stops at once.
        ts = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0]
        zr = 0.02 .+ 0.015 .* (1 .- exp.(-ts ./ 10))
        zcb(z) = ZCBPrice.(exp.(-z .* ts), ts)
        starts = (fill(0.03, length(ts)), [0.03, 0.01, 0.04, 0.005, 0.05, 0.01, 0.04, 0.02, 0.03], fill(1.0, length(ts)))
        splines = (Spline.Linear(), Spline.Quadratic(), Spline.Cubic(), Spline.BSpline(3), Spline.PCHIP(), Spline.Akima(), Spline.MonotoneConvex())
        for spline in splines, z in (zr, zr .+ 0.01 .* (ts ./ 30 .- 0.5))
            direct = fit(spline, zcb(z))
            @test maximum(abs(discount(direct, ti) - exp(-zi * ti)) for (ti, zi) in zip(ts, z)) < 1.0e-8
            @test all(isequal(fit(ZeroRateCurve(z0, ts, spline), zcb(z)), direct) for z0 in starts)
        end
        for (z0, tenor) in ((1.0, 30.0), (0.1, 300.0), (100.0, 10.0))
            @test discount(fit(ZeroRateCurve([z0], [tenor], Spline.Linear()), [ZCBPrice(0.9, tenor)]), tenor) ≈ 0.9 rtol = 1.0e-8
        end
        fl = fit(ZeroRateCurve(fill(0.01, length(t)), t, Spline.Linear()), qs)
        @test fl.rates ≈ target atol = 1.0e-5
        # more quotes than knots: a least-squares fit at the curve's own tenors
        over = fit(ZeroRateCurve(fill(0.01, 3), [1.0, 3.0, 10.0], Spline.Linear()), qs)
        @test over.tenors == [1.0, 3.0, 10.0] && eltype(over.rates) === Float64
        @test all(isfinite, over.rates)
        # knot rates are the variables: a knot curve's fit takes no optics
        zrc0 = ZeroRateCurve(fill(0.01, length(t)), t)
        @test_throws MethodError fit(zrc0, qs; variables = ((@optic(_.rates[1]),),))
        @test fit(zrc0, qs, Fit.Loss(x -> x^2); solve_kwargs = (; g_tol = 1.0e-12)) isa Yield.MonotoneConvex
    end

    @testset "single-pass quote iterators" begin
        # Each fit collects its quotes once, so a stateful iterator gives the vector's result.
        t = [1.0, 3.0, 7.0]
        qs = ZCBPrice.([0.97, 0.91, 0.82], t)
        once() = Iterators.Stateful(qs)
        @test isequal(fit(Spline.Linear(), once()), fit(Spline.Linear(), qs))
        @test isequal(fit(Spline.Cubic(), once(), Fit.Loss(abs2)), fit(Spline.Cubic(), qs, Fit.Loss(abs2)))
        @test isequal(fit(Spline.Linear(), once(), Fit.Bootstrap()), fit(Spline.Linear(), qs, Fit.Bootstrap()))
        c = ZeroRateCurve([0.03, 0.03, 0.03], t, Spline.Linear())
        @test isequal(fit(c, once()), fit(c, qs))
        @test discount(fit(Yield.NelsonSiegel(), once()), 5.0) == discount(fit(Yield.NelsonSiegel(), qs), 5.0)
        sw = Yield.SmithWilson(ufr = 0.04, α = 0.1)
        @test discount(fit(sw, once()), 5.0) == discount(fit(sw, qs), 5.0)
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
            @test_throws Base.CanonicalIndexError knot_rates(c) .= 0.0
            @test_throws Base.CanonicalIndexError sort!(knot_tenors(c); rev = true)
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
