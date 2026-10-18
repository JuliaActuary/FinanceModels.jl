using ForwardDiff

@testset "implicit roots" begin
    root = FinanceModels.__implicit_root
    primal = FinanceModels.__primal

    # `g_primal` is `g` with primal parameters, so that the solver can differentiate it in `x`.
    # a·eˣ - b = 0 has the root log(b/a); the residual is curved in x and mixes x with a.
    function log_ratio(a, b; max_order = 1)
        g(x) = a * exp(x) - b
        pa, pb = primal(a), primal(b)
        return root(g, x -> pa * exp(x) - pb, 0.0; who = "test root", scale = x -> 1.0, max_order)
    end

    # x³ + a·x - b = 0, with the implicit derivatives written out: F = ∂g/∂x = 3x² + a.
    function cubic(a, b; max_order = 1)
        g(x) = x^3 + a * x - b
        pa, pb = primal(a), primal(b)
        return root(g, x -> x^3 + pa * x - pb, 1.0; who = "test root", scale = x -> 1.0, max_order)
    end
    function cubic_derivatives(a, b)
        x = cubic(a, b)
        F = 3x^2 + a
        xa, xb = -x / F, 1 / F
        # differentiate xa = -x/F and xb = 1/F along each input, with F' = 6x·x' + ∂F/∂input
        Fa, Fb = 6x * xa + 1, 6x * xb
        xaa = -(xa * F - x * Fa) / F^2
        xab = -(xb * F - x * Fb) / F^2
        xbb = -Fb / F^2
        return [xa, xb], [xaa xab; xab xbb]
    end

    @testset "first order is unchanged" begin
        @test ForwardDiff.gradient(θ -> log_ratio(θ...), [2.0, 3.0]) ≈ [-1 / 2, 1 / 3]
        grad, _ = cubic_derivatives(0.5, 2.0)
        @test ForwardDiff.gradient(θ -> cubic(θ...), [0.5, 2.0]) ≈ grad
        # the default supports first order only
        @test_throws "first-order ForwardDiff derivatives only" ForwardDiff.hessian(θ -> log_ratio(θ...), [2.0, 3.0])
    end

    @testset "second order" begin
        θ = [2.0, 3.0]
        @test ForwardDiff.hessian(θ -> log_ratio(θ...; max_order = 2), θ) ≈ [1 / 4 0; 0 -1 / 9] atol = 1.0e-14
        grad, hess = cubic_derivatives(0.5, 2.0)
        @test ForwardDiff.gradient(θ -> cubic(θ...; max_order = 2), [0.5, 2.0]) ≈ grad
        @test ForwardDiff.hessian(θ -> cubic(θ...; max_order = 2), [0.5, 2.0]) ≈ hess rtol = 1.0e-12
        # nested derivatives with different tags: the mixed partial
        mixed = ForwardDiff.derivative(a -> ForwardDiff.derivative(b -> cubic(a, b; max_order = 2), 2.0), 0.5)
        @test mixed ≈ hess[1, 2] rtol = 1.0e-12
        # one chord step alone would get the second derivative wrong: it is not the answer by accident
        one_step(a, b) = (x = cubic(primal(a), primal(b)); x - ((x^3 + a * x - b) - primal(x^3 + a * x - b)) / (3x^2 + primal(a)))
        @test !isapprox(ForwardDiff.hessian(θ -> one_step(θ...), [0.5, 2.0]), hess; rtol = 1.0e-3)
        # the value of a dual result is the primal root exactly
        h = ForwardDiff.Dual{:outer}(ForwardDiff.Dual{:inner}(0.5, 1.0), ForwardDiff.Dual{:inner}(1.0, 0.0))
        @test primal(cubic(h, 2.0; max_order = 2)) === cubic(0.5, 2.0)
        # deeper nesting than supported throws
        @test_throws "up to order 2" ForwardDiff.derivative(
            a -> ForwardDiff.derivative(b -> ForwardDiff.derivative(c -> cubic(a + b + c, 2.0; max_order = 2), 0.0), 0.0), 0.5
        )
    end

    @testset "a sequence of roots" begin
        # x₁ = log(b/a), then x₂ solves e^x₂ - c·x₁ = 0, so x₂ = log(c·log(b/a)): each solve reads the
        # previous root, carrying its second-order partials
        function sequential(θ)
            a, b, c = θ
            x1 = log_ratio(a, b; max_order = 2)
            g(x) = exp(x) - c * x1
            p = primal(c * x1)
            return root(g, x -> exp(x) - p, 0.0; who = "test root", scale = x -> 1.0, max_order = 2)
        end
        closed(θ) = log(θ[3] * log(θ[2] / θ[1]))
        θ = [1.5, 7.0, 0.8]
        @test sequential(θ) ≈ closed(θ)
        @test ForwardDiff.gradient(sequential, θ) ≈ ForwardDiff.gradient(closed, θ) rtol = 1.0e-12
        @test ForwardDiff.hessian(sequential, θ) ≈ ForwardDiff.hessian(closed, θ) rtol = 1.0e-10
    end
end
