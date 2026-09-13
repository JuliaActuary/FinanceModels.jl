module ProjectionModelTests

using Test, FinanceCore, FinanceModels, Transducers, ForwardDiff

struct UndeclaredProjectionContract <: FinanceCore.AbstractContract end
struct DeclaredProjectionContract <: FinanceCore.AbstractContract
    key::Symbol
end
FinanceModels.model_requirements(c::DeclaredProjectionContract) = (c.key => Yield.AbstractYieldModel,)

struct IterableProjectionContract{R} <: FinanceCore.AbstractContract
    requirements::R
end
FinanceModels.model_requirements(c::IterableProjectionContract) = c.requirements

@testset "Projection model requirements" begin
    curve = Yield.Constant(Continuous(0.04))
    credit = Yield.Constant(Continuous(0.06))
    fixed = Bond.Fixed(0.05, Periodic(2), 5.0)
    floating = Bond.Floating(0.001, Periodic(4), 5.0, :index)
    swap = InterestRateSwap(curve, 5.0; model_key = :index)
    transformed = floating |> Map(-) |> Map(cf -> cf * 2)
    forward_start = Forward(1.0, floating)
    portfolio = [swap, fixed, forward_start]
    store = Dict(:index => curve)

    @test model_requirements(fixed) == ()
    @test model_requirements(Cashflow(1.0, 2.0)) == ()
    for c in (floating, swap, transformed, forward_start)
        @test model_requirements(c) == (:index => Yield.AbstractYieldModel,)
        @test collect(Projection(c; index = curve)) == collect(Projection(c, store))
        @test present_value(credit, Projection(c; index = curve)) ≈
            present_value(credit, Projection(c, store))
        auto(r) = present_value(credit, Projection(c; index = Yield.Constant(Continuous(r))))
        explicit(r) = present_value(credit, Projection(c, Dict(:index => Yield.Constant(Continuous(r)))))
        @test ForwardDiff.derivative(auto, 0.04) ≈ ForwardDiff.derivative(explicit, 0.04)
    end
    @test present_value(curve, Projection(swap; index = curve)) ≈ 0.0 atol = 1.0e-12
    @test collect(model_requirements(portfolio)) == [:index => Yield.AbstractYieldModel, :index => Yield.AbstractYieldModel]
    @test present_value(credit, Projection(portfolio; index = curve)) ≈
        sum(present_value(credit, Projection(c, store)) for c in portfolio)
    @test isempty(model_requirements(FinanceCore.AbstractContract[]))
    @test isempty(Projection(FinanceCore.AbstractContract[]; index = curve).model)
    @test Projection(fixed).model isa NullModel
    @test collect(Projection(fixed; index = curve)) == collect(fixed)

    # Distinct index keys may share a curve, but remain explicit in the protocol.
    other = Bond.Floating(0.0, Periodic(2), 3.0, "other")
    both = FinanceCore.Composite(floating, other)
    @test model_requirements(both) == (:index => Yield.AbstractYieldModel, "other" => Yield.AbstractYieldModel)
    @test collect(Projection(both; index = curve)) ==
        collect(Projection(both, Dict(:index => curve, "other" => curve)))

    pair = FX.Pair(:EUR, :USD)
    fx = FX.Forwards(pair, 1.08, credit, curve)
    converted = FX.Converted(transformed, pair, :fx)
    @test model_requirements(converted) == (:fx => FX.AbstractFXModel, :index => Yield.AbstractYieldModel)
    @test_throws "explicit model store" Projection(converted; index = curve)
    @test_throws "explicit model store" Projection(converted; index = fx)
    @test_throws "explicit model store" Projection(floating; index = fx)
    explicit_fx = Projection(converted, Dict(:fx => fx, :index => curve))
    @test present_value(credit, explicit_fx) ≈ 1.08 * present_value(curve, Projection(transformed, store))
    fixed_fx = FX.Converted(fixed, pair, :fx)
    @test collect(Projection(fixed_fx; index = fx)) == collect(Projection(fixed_fx, Dict(:fx => fx)))
    @test model_requirements(FX.BasisSwapLeg(pair, [Cashflow(1.0, 2.0)])) == ()

    # An extension must declare its needs; explicit stores keep working without
    # adopting the convenience protocol, and are not replaced by inferred wiring.
    unknown = UndeclaredProjectionContract()
    @test_throws ArgumentError model_requirements(unknown)
    @test_throws r"UndeclaredProjectionContract.*Define FinanceModels.model_requirements.*returning \(\).*explicit model store" Projection(unknown; index = curve)
    @test Projection(unknown, store).model === store
    @test Projection(DeclaredProjectionContract(:custom); index = curve).model[:custom] === curve

    @testset "Runtime portfolios stay lazy through wrappers" begin
        small = fill(floating, 1)
        large = fill(floating, 10_000)
        yield_requirement = :index => Yield.AbstractYieldModel
        fx_requirement = :fx => FX.AbstractFXModel
        for (wrap, expected) in (
                (identity, (yield_requirement,)),
                (cs -> FinanceCore.Composite(fixed, cs), (yield_requirement,)),
                (cs -> FinanceCore.Composite(cs, floating), (yield_requirement, yield_requirement)),
                (cs -> Forward(1.0, FinanceCore.Composite(fixed, cs)), (yield_requirement,)),
                (cs -> cs |> Map(-), (yield_requirement,)),
                (cs -> FX.Converted(cs, pair, :fx), (fx_requirement, yield_requirement)),
                (cs -> FinanceCore.Composite(FX.Converted(cs, pair, :fx), cs), (fx_requirement, yield_requirement, yield_requirement)),
            )
            small_requirements = @inferred model_requirements(wrap(small))
            large_requirements = @inferred model_requirements(wrap(large))
            @test !(large_requirements isa Tuple)
            @test typeof(small_requirements) === typeof(large_requirements)
            @test Tuple(small_requirements) == expected
        end
        @test count(==(:index => Yield.AbstractYieldModel), model_requirements(large)) == length(large)
        @test Projection(large; index = curve).model == store
        @test collect(model_requirements(reshape(large, 100, 100))) == collect(model_requirements(large))
        nested = [Bond.Floating[], small, Bond.Floating[], small]
        @test collect(model_requirements(nested)) == fill(:index => Yield.AbstractYieldModel, 2)

        # Obtaining requirements must not visit array elements, even through wrappers.
        for wrap in (identity, cs -> FinanceCore.Composite(fixed, cs), cs -> FX.Converted(cs, pair, :fx))
            cs = FinanceCore.AbstractContract[floating, unknown]
            requirements = model_requirements(wrap(cs))
            @test first(requirements) in (:index => Yield.AbstractYieldModel, :fx => FX.AbstractFXModel)
            @test_throws ArgumentError collect(requirements)
        end
    end

    @testset "Every occurrence is validated in one pass" begin
        requirements = Iterators.Stateful(
            [
                :index => Yield.AbstractYieldModel,
                :index => typeof(curve),
                :other => Yield.AbstractYieldModel,
            ]
        )
        custom = IterableProjectionContract(requirements)
        @test Projection(custom; index = curve).model == Dict(:index => curve, :other => curve)
        @test isempty(requirements)

        # Repeated keys may impose incompatible constraints, in either order.
        for types in ((Yield.AbstractYieldModel, FX.AbstractFXModel), (FX.AbstractFXModel, Yield.AbstractYieldModel))
            for index in (curve, fx)
                requirements = Iterators.Stateful([:index => types[1], :index => types[2]])
                custom = IterableProjectionContract(requirements)
                invalid_type = index isa types[1] ? types[2] : types[1]
                @test_throws "key :index requires $invalid_type" Projection(custom; index)
            end
        end
        conflict = FX.Converted(fixed, pair, :index)
        cs = FinanceCore.AbstractContract[fill(floating, 1_000); conflict]
        for wrap in (identity, cs -> FinanceCore.Composite(fixed, cs), cs -> Forward(1.0, FinanceCore.Composite(fixed, cs)), cs -> cs |> Map(-))
            @test_throws "key :index requires $(FX.AbstractFXModel)" Projection(wrap(cs); index = curve)
        end
        @test_throws "key :index requires $(Yield.AbstractYieldModel)" Projection(FX.Converted(floating, pair, :index); index = fx)
    end
end

end # module ProjectionModelTests
