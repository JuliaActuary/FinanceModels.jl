using ForwardDiff

@testset "ZeroRateCurve" begin
    rates = [0.02, 0.03, 0.035, 0.04]
    tenors = [1.0, 2.0, 5.0, 10.0]

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
            @test FinanceCore.rate(z) ≈ 0.03 atol = 1e-10
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
            @test discount(zrc, t) ≈ exp(-r_interp * t) atol = 1e-10
        end
    end

    @testset "from AbstractYieldModel" begin
        @testset "round-trip with Constant" begin
            c = Yield.Constant(0.05)
            tenors = [1.0, 2.0, 5.0, 10.0]
            zrc = ZeroRateCurve(c, tenors)
            for t in tenors
                @test discount(zrc, t) ≈ discount(c, t) atol = 1e-10
            end
        end

        @testset "round-trip with NelsonSiegel" begin
            ns = Yield.NelsonSiegel(1.0, 0.04, -0.02, 0.01)
            tenors = [1.0, 2.0, 5.0, 10.0, 20.0]
            zrc = ZeroRateCurve(ns, tenors)
            for t in tenors
                @test discount(zrc, t) ≈ discount(ns, t) atol = 1e-10
            end
        end

        @testset "explicit spline kwarg" begin
            c = Yield.Constant(0.04)
            tenors = [1.0, 5.0, 10.0]
            zrc = ZeroRateCurve(c, tenors; spline=Spline.Linear())
            for t in tenors
                @test discount(zrc, t) ≈ discount(c, t) atol = 1e-10
            end
        end

        @testset "unsorted tenors are sorted automatically" begin
            c = Yield.Constant(0.04)
            zrc_sorted = ZeroRateCurve(c, [1.0, 5.0, 10.0])
            zrc_unsorted = ZeroRateCurve(c, [10.0, 1.0, 5.0])
            for t in [1.0, 3.0, 5.0, 7.0, 10.0]
                @test discount(zrc_sorted, t) ≈ discount(zrc_unsorted, t) atol = 1e-10
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
                @test discount(zrc, t) ≈ exp(-r * t) atol = 1e-8
            end
        end

        @testset "two-point case matches linear" begin
            r2 = [0.03, 0.05]
            t2 = [1.0, 5.0]
            zrc_lin = ZeroRateCurve(r2, t2, Spline.Linear())
            zrc_cub = ZeroRateCurve(r2, t2, Spline.Cubic())
            @test discount(zrc_lin, 3.0) ≈ discount(zrc_cub, 3.0) atol = 1e-6
        end
    end

    @testset "eager-build: ForwardDiff pass-through" begin
        # The eager build runs once with the input rate type; for Dual-typed
        # rates, the model itself becomes Dual-typed and propagates through
        # discount. At an exact knot t = tenors[k], discount = exp(-rates[k] * t),
        # so ∂discount/∂rates[k] = -t · discount.
        tenors_ad = [1.0, 2.0, 5.0, 10.0]
        rates_ad  = [0.02, 0.03, 0.035, 0.04]
        for spl in (Spline.Linear(), Spline.MonotoneConvex())
            f(r) = discount(ZeroRateCurve([r, rates_ad[2:end]...], tenors_ad, spl), 1.0)
            ad = ForwardDiff.derivative(f, 0.02)
            @test ad ≈ -1.0 * exp(-0.02 * 1.0) atol = 1e-12
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
            # hidden backing storage and cache, in both public and private listings
            @test_throws ArgumentError zrc.rates._data
            @test propertynames(zrc.rates) == ()
            @test propertynames(zrc.rates, true) == ()
            @test propertynames(zrc) == (:rates, :tenors, :spline)
            @test propertynames(zrc, true) == (:rates, :tenors, :spline)
            @test_throws ArgumentError zrc._model
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
            r[2] = 0.10; t[2] = 3.0            # mutate the caller's inputs
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
        @test keys(CB.getproperties(z)) == (:rates, :tenors, :spline)
        @test CB.setproperties(z, CB.getproperties(z)) == z
        raw = CB.constructorof(typeof(z))(CB.getfields(z)...)
        @test raw == z && discount(raw, 3.5) == discount(z, 3.5)
        @test Accessors.mapproperties(identity, z) == z
        @test Accessors.getall(z, Accessors.Properties()) == (z.rates, z.tenors, z.spline)
        newrates = [0.03, 0.04, 0.05]
        zp = Accessors.setall(z, Accessors.Properties(), (newrates, z.tenors, z.spline))
        @test zp == ZeroRateCurve(newrates, [1.0, 2.0, 5.0], Spline.Linear())
        @test discount(zp, 5.0) ≈ exp(-0.05 * 5.0)
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
        # the cache cannot be patched (RHS must not read `zrc._model`: that getter throws first)
        @test_throws ArgumentError zrc._model
        @test_throws ArgumentError (Accessors.@set zrc._model = nothing)
        @test_throws ArgumentError CB.setproperties(zrc, (_model = nothing,))
        @test_throws ArgumentError CB.setproperties(zrc, (foo = 1,))
        # positional 4-arg reconstruction path (what `setproperties`/`fit` call internally)
        sp = CB.setproperties(zrc, (rates = [0.03, 0.04, 0.05],))
        @test sp == ZeroRateCurve([0.03, 0.04, 0.05], t, Spline.MonotoneConvex())
        @test discount(sp, 5.0) ≈ exp(-0.05 * 5.0)
        # the cache is retyped when rates become Duals via Accessors
        dz = Accessors.@set zrc.rates[1] = ForwardDiff.Dual(0.02, 1.0)
        @test eltype(dz.rates) <: ForwardDiff.Dual
        @test discount(dz, 1.0) isa ForwardDiff.Dual
        @test ForwardDiff.partials(discount(dz, 1.0), 1) ≈ -1.0 * exp(-0.02 * 1.0) atol = 1e-12
        # no public 4-arg constructor
        @test_throws MethodError ZeroRateCurve(r, t, Spline.Linear(), getfield(zrc, :_model))
    end

    # ─── Validation ───────────────────────────────────────────────────────────

    @testset "validation" begin
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
        @test Spline.BSpline(3).order == 3 && Spline.Cubic().order == 3
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
            Spline.Linear() => 2, Spline.Quadratic() => 2, Spline.Cubic() => 2, Spline.BSpline(3) => 2,
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
        # a single MonotoneConvex knot is a flat curve
        z1 = ZeroRateCurve([0.03], [2.0])
        @test discount(z1, 1.0) ≈ exp(-0.03 * 1.0)
        @test discount(z1, 4.0) ≈ exp(-0.03 * 4.0)
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
        ze = ZeroRateCurve(c, [1e-20, 1.0, 2e4])
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
        @test zb isa ZeroRateCurve && zb.rates == target && zb.tenors == t
        @test discount(zb, 2.0) ≈ exp(-0.025 * 2.0)
        @test zrc0.rates == fill(0.01, length(t))                  # original untouched
        @test Accessors.modify(x -> 2x, zrc0, o).rates == fill(0.02, length(t))
        @test_throws ArgumentError Accessors.setall(zrc0, o, [NaN, 0.0, 0.0, 0.0])   # still validated
        fitted = fit(zrc0, qs)
        @test fitted isa ZeroRateCurve
        @test fitted.tenors == t
        @test maximum(abs, present_value(fitted, q.instrument) - q.price for q in qs) < 1e-6
        fl = fit(ZeroRateCurve(fill(0.01, length(t)), t, Spline.Linear()), qs)
        @test fl.rates ≈ target atol = 1e-5
    end
end
