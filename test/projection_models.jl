using ForwardDiff

struct UndeclaredProjectionContract <: FinanceCore.AbstractContract end
struct DeclaredProjectionContract <: FinanceCore.AbstractContract
    key::Symbol
end
FinanceModels.model_requirements(c::DeclaredProjectionContract) = (c.key => Yield.AbstractYieldModel,)

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
    @test model_requirements(portfolio) == (:index => Yield.AbstractYieldModel, :index => Yield.AbstractYieldModel)
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
    @test_throws MethodError model_requirements(unknown)
    @test_throws MethodError Projection(unknown; index = curve)
    @test Projection(unknown, store).model === store
    @test Projection(DeclaredProjectionContract(:custom); index = curve).model[:custom] === curve
end
