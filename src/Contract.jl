# Extending the Core contracts from FinanceCore

### Bonds

"""
The `Bond` module provide a number of fixed-income contracts and related methods.
"""
module Bond
    import ..FinanceCore: Cashflow, Quote, AbstractContract, maturity, Timepoint
    using ..FinanceCore

    using FinanceCore: Periodic, Continuous, Rate

    export ZCBYield, ZCBPrice, ParSwapYield, ParYield, CMTYield

    abstract type AbstractBond <: AbstractContract end
    maturity(b::AbstractBond) = b.maturity


    """
        ZCBPrice(discount,maturity)
        ZCBPrice(yield::Vector)

    Takes spot/zero discount factors and returns a `Quote` for the cashflow occuring at the given `maturity`.

    Use broadcasting to create a set of quotes given a collection of prices and maturities, e.g. `ZCBPrice.(FinanceModels,maturities)`.

    See also [`ZCBYield`](@ref)

    # Examples

    ```julia-repl

    julia> ZCBPrice(0.5,10)
    Quote{Float64, Cashflow{Float64, Int64}}(0.5, Cashflow{Float64, Int64}(1.0, 10))

    julia> ZCBPrice([0.9,0.8,0.75])
    3-element Vector{Quote{Float64, Cashflow{Float64, Int64}}}:
     Quote{Float64, Cashflow{Float64, Int64}}(0.9, Cashflow{Float64, Int64}(1.0, 1))
     Quote{Float64, Cashflow{Float64, Int64}}(0.8, Cashflow{Float64, Int64}(1.0, 2))
     Quote{Float64, Cashflow{Float64, Int64}}(0.75, Cashflow{Float64, Int64}(1.0, 3))
     
    ```

    """
    ZCBPrice(price, time) = Quote(price, Cashflow(1.0, time))


    """
        ZCBYield(yield,maturity)
        ZCBYield(yield::Vector)

    Returns a `Quote` for the cashflow occuring at the given `maturity` and the quoted value is derived from the given `yield`.

    Takes zero (sometimes called "spot") rates. Assumes annual effective compounding (`Periodic(1)``) unless given a `Rate` with a different compounding frequency.

    Use broadcasting to create a set of quotes given a collection of FinanceModels and maturities, e.g. `ZCBYield.(FinanceModels,maturities)`.

    See also [`ZCBPrice`](@ref)

    # Examples

    ```julia-repl
    julia> ZCBYield(0.05,30)
    Quote{Float64, Cashflow{Float64, Int64}}(0.23137744865585788, Cashflow{Float64, Int64}(1.0, 30))

    julia> ZCBYield(Periodic(0.05,1),30)
    Quote{Float64, Cashflow{Float64, Int64}}(0.23137744865585788, Cashflow{Float64, Int64}(1.0, 30))

    julia> ZCBYield(Continuous(0.05),30)
    Quote{Float64, Cashflow{Float64, Int64}}(0.22313016014842982, Cashflow{Float64, Int64}(1.0, 30))

    julia> ZCBYield([0.04,0.05,0.045])
    3-element Vector{Quote{Float64, Cashflow{Float64, Int64}}}:
     Quote{Float64, Cashflow{Float64, Int64}}(0.9615384615384615, Cashflow{Float64, Int64}(1.0, 1))
     Quote{Float64, Cashflow{Float64, Int64}}(0.9070294784580498, Cashflow{Float64, Int64}(1.0, 2))
     Quote{Float64, Cashflow{Float64, Int64}}(0.8762966040549094, Cashflow{Float64, Int64}(1.0, 3))
    ```
    """
    ZCBYield(yield, time) = Quote(discount(yield, time), Cashflow(1.0, time))


    """
        Bond.Fixed(coupon_rate,frequency<:FinanceCore.Frequency,maturity)

    An object representing a fixed coupon bond. `coupon_rate` / `frequency` is the actual payment amount.

    Note that there are a number of convienience constructors which return a Quote for a `Bond.Fixed`: 

    - [`ParYield`](@ref)
    - [`ParSwapYield`](@ref)
    - [`CMTYield`](@ref)
    - [`OISYield`](@ref)

    See also [`FinanceCore.Quote`](@ref).

    # Examples

    ```julia-repl
    julia> Bond.Fixed(0.05,Periodic(2),3)
    FinanceModels.Bond.Fixed{Periodic, Float64, Int64}(0.05, Periodic(2), 3)

    julia> Bond.Fixed(0.05,Periodic(2),3) |> collect
    6-element Vector{Cashflow{Float64, Float64}}:
     Cashflow{Float64, Float64}(0.025, 0.5)
     Cashflow{Float64, Float64}(0.025, 1.0)
     Cashflow{Float64, Float64}(0.025, 1.5)
     Cashflow{Float64, Float64}(0.025, 2.0)
     Cashflow{Float64, Float64}(0.025, 2.5)
     Cashflow{Float64, Float64}(1.025, 3.0)


    julia> ParYield(0.05,10)
    Quote{Float64, FinanceModels.Bond.Fixed{Periodic, Float64, Int64}}(1.0, FinanceModels.Bond.Fixed{Periodic, Float64, Int64}(0.05, Periodic(2), 10))
    ```

    """
    struct Fixed{F <: FinanceCore.Frequency, N <: Real, M <: Timepoint} <: AbstractBond
        coupon_rate::N # coupon_rate / frequency is the actual payment amount
        frequency::F
        maturity::M
    end

    function Base.isapprox(a::Fixed, b::Fixed)
        return isapprox(a.coupon_rate, b.coupon_rate) && ==(a.frequency, b.frequency) && isapprox(a.maturity, b.maturity)
    end

    """
        Bond.Floating(coupon_rate,frequency<:FinanceCore.Frequency,maturity,model_key)

    An object representing a floating coupon bond. (`coupon_rate` + reference rate) / `frequency` is the actual payment amount, where the reference rate requires a `Projection` with a key/value pair where the key is the `model_key` argument and the value is the model which produces the reference rate.


    See also [`FinanceCore.Quote`](@ref).

    # Examples

    ```julia-repl
    julia> p = Projection(
            Bond.Floating(0.02, Periodic(1), 3.0, "SOFR"),
            Dict("SOFR" => Yield.Constant(0.05)),  # note the key/value store used for the model in the projection
            CashflowProjection(),
        );

    julia> collect(p)
    3-element Vector{Cashflow{Float64, Float64}}:
        Cashflow{Float64, Float64}(0.07000000000000005, 1.0)
        Cashflow{Float64, Float64}(0.07000000000000005, 2.0)
        Cashflow{Float64, Float64}(1.07, 3.0)
    ```
    """
    struct Floating{F <: FinanceCore.Frequency, N <: Real, M <: Timepoint, K} <: AbstractBond
        coupon_rate::N # coupon_rate / frequency is the actual payment amount
        frequency::F
        maturity::M
        key::K
    end

    __coerce_periodic(y::Periodic) = y
    __coerce_periodic(y::T) where {T <: Int} = Periodic(y)

    """
    ParYield(yield, maturity; frequency=Periodic(2))
    ParYield(yield::Vector)

    Takes bond equivalent FinanceModels, and assumes that instruments <= one year maturity pay no coupons and that the rest pay semi-annual. Alternative, you may pass a `Rate` as the yield and the coupon frequency will be inferred from the `Rate`'s frequency. 

    Use broadcasting to create a set of quotes given a collection of FinanceModels and maturities, e.g. `ParYield.(FinanceModels,maturities)`.

    # Examples

    ```julia-repl
    julia> ParYield(0.05,10)
    Quote{Float64, FinanceModels.Bond.Fixed{Periodic, Float64, Int64}}(1.0, FinanceModels.Bond.Fixed{Periodic, Float64, Int64}(0.05, Periodic(2), 10))
    ```
    """
    function ParYield(yield, maturity; frequency = Periodic(2))
        # assume the frequency is two or infer it from the yield
        frequency = __coerce_periodic(frequency)
        price = 1.0 # by definition for a par bond
        coupon_rate = rate(frequency(yield))
        return Quote(price, Fixed(coupon_rate, frequency, maturity))
    end
    function ParYield(yield::Rate{N, T}, maturity; frequency = Periodic(2)) where {T <: Periodic, N}
        frequency = yield.compounding
        price = 1.0 # by definition for a par bond
        coupon_rate = rate(frequency(yield))
        return Quote(price, Fixed(coupon_rate, frequency, maturity))
    end

    """
        ParSwapYield(yield, maturity; frequency=Periodic(4))

    Same as [`ParYield`](@ref), except the `frequency` is four times per period by default.
    """
    function ParSwapYield(yield, maturity; frequency = Periodic(4))
        frequency = __coerce_periodic(frequency)
        return ParYield(yield, maturity; frequency = frequency)
    end

    """
        CMTYield(yield,maturity)
        CMTYield(yield::Vector)

    Returns a `Quote` for the correpsonding bond implied by the given bond equivalent `yield`, and assumes that instruments <= one year `maturity`` pay no coupons and that the rest pay semi-annual.

    Use broadcasting to create a set of quotes given a collection of FinanceModels and maturities, e.g. `CMTYield.(FinanceModels,maturities)`.

    See also [`FinanceCore.Quote`](@ref), [`Bond.Fixed`](@ref)

    # Examples

    ```
    julia> CMTYield(0.05,10)
    Quote{Float64, FinanceModels.Bond.Fixed{Periodic, Float64, Int64}}(1.0, FinanceModels.Bond.Fixed{Periodic, Float64, Int64}(0.05, Periodic(2), 10))
    ```
    """
    function CMTYield(yield, maturity)
        # Assume maturity < 1 don't pay coupons and are therefore discount bonds
        # Assume maturity > 1 pay coupons and are therefore par bonds
        frequency = Periodic(2)
        r, v = if maturity ≤ 1
            Periodic(0.0, 1), discount(yield, maturity)
        else
            # coupon paying par bond
            frequency(yield), 1.0
        end
        return Quote(v, Fixed(rate(r), r.compounding, maturity))
    end

    """
    OISYield(yield, maturity)

    Returns the implied `Quote` for the fixed bond implied by the given `yield` and `maturity`. Assumes that maturities less than or equal to 12 months are settled once (per Hull textbook, 4.7), otherwise quarterly and that the FinanceModels given are bond equivalent.

    Use broadcasting to create a set of quotes given a collection of FinanceModels and maturities, e.g. `OISYield.(FinanceModels,maturities)`.

    See also [`FinanceCore.Quote`](@ref), [`Bond.Fixed`](@ref)

    # Examples

    ```
    julia> OISYield(0.05,10)
    Quote{Float64, FinanceModels.Bond.Fixed{Periodic, Float64, Int64}}(1.0, FinanceModels.Bond.Fixed{Periodic, Float64, Int64}(0.05, Periodic(4), 10))
    ```

    """
    function OISYield(yield, maturity)

        if maturity <= 1
            return Quote(discount(yield, maturity), Fixed(0.0, Periodic(1), maturity))
        else
            frequency = Periodic(4)
            r = frequency(yield)
            return Quote(1.0, Fixed(rate(r), frequency, maturity))
        end
    end

    """
        ForwardYields(yields,times) 

    Returns a vector of `Quote` corresponding to the yield at the given forward times. 
        
    # Examples
    ```julia-repl
    julia> FinanceModels.Bond.ForwardYields([0.01,0.02],[1.,3.])
    2-element Vector{Quote{Float64, Cashflow{Float64, Float64}}}:
     Quote{Float64, Cashflow{Float64, Float64}}(0.9900990099009901, Cashflow{Float64, Float64}(1.0, 1.0))
     Quote{Float64, Cashflow{Float64, Float64}}(0.9423223345470445, Cashflow{Float64, Float64}(1.0, 3.0))
    ```
    """
    function ForwardYields(yields, times = eachindex(yields))
        df = 1.0
        t_prior = 0.0
        return map(zip(yields, times)) do (y, t)
            df *= discount(y, t - t_prior)
            t_prior = t
            Quote(
                df,
                Cashflow(1.0, t)
            )
        end
    end


    # Bond utility funcs

    """
        coupon_times(maturity, frequency)

    Generate coupon times for a bond with the given `maturity` and `frequency`.

    # Arguments
    - `maturity::Real`: The maturity of the bond.
    - `frequency::Real`: The coupon frequency of the bond.

    # Returns
    - An array of coupon times for the bond.

    # Examples
    ```julia-repl
    julia> Bond.coupon_times(10, 2)
    0.5:0.5:10.0
    julia> Bond.coupon_times(Bond.Fixed(0.05,Periodic(4),20))
    0.25:0.25:20.0
    ````
    """
    function coupon_times(maturity, frequency)
        Δt = min(1 / frequency, maturity)
        times = maturity:-Δt:0
        if iszero(last(times))
            return reverse(times[1:(end - 1)])
        else
            return reverse(times)
        end
    end
    coupon_times(b::AbstractBond) = coupon_times(b.maturity, b.frequency.frequency)


    for op in (:ZCBPrice, :ZCBYield, :ParYield, :ParSwapYield, :CMTYield, :ForwardYield)
        eval(
            quote
                $op(x::Vector; kwargs...) = $op.(x, float.(eachindex(x)); kwargs...)
            end
        )
    end


end

"""
    FinanceCore.internal_rate_of_return(q::Quote)

Return the internal rate of return (yield to maturity) implied by the quote's price and cashflows.
"""
function FinanceCore.internal_rate_of_return(q::Quote)
    cashflows = collect(q.instrument)
    time_zero = zero(FinanceCore.timepoint(first(cashflows)))
    rate = FinanceCore.internal_rate_of_return([Cashflow(-q.price, time_zero); cashflows])
    return rate
end

"""
    CommonEquity()

A singleton type representing a unit stock.

See also: [`Option`](@ref).

"""
struct CommonEquity <: FinanceCore.AbstractContract end

"""

"""
module Option
import ..FinanceCore: AbstractContract, Timepoint


"""
    EuroCall(contract,strike,maturity)

A European call option on the given contract with the given strike and maturity.

# Arguments
 - contract::AbstractContract -  The underlying contract.
 - strike::Real -  The strike price.
 - maturity::Union{Real,Date} -  The maturity of the option.

 Supertype Hierarchy
≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    EuroCall{S,K,M} <: FinanceCore.AbstractContract <: Any

"""
struct EuroCall{S <: AbstractContract, K <: Real, M <: Timepoint} <: AbstractContract
    underlying::S
    strike::K
    maturity::M
end

"""
    EuroPut(contract,strike,maturity)

A European put option on the given contract with the given strike and maturity.

# Arguments
 - contract::AbstractContract -  The underlying contract.
 - strike::Real -  The strike price.
 - maturity::Union{Real,Date} -  The maturity of the option.

 Supertype Hierarchy
≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    EuroPut{S,K,M} <: FinanceCore.AbstractContract <: Any

"""
struct EuroPut{S <: AbstractContract, K <: Real, M <: Timepoint} <: AbstractContract
    underlying::S
    strike::K
    maturity::M
end

"""
    ZCBCall(expiry, bond_maturity, strike)

A European call option on a zero-coupon bond.
The holder has the right to buy at time `expiry` a ZCB maturing at
`bond_maturity` for `strike`.
"""
struct ZCBCall{T <: Real, S <: Real, K <: Real} <: AbstractContract
    expiry::T
    bond_maturity::S
    strike::K
end

"""
    ZCBPut(expiry, bond_maturity, strike)

A European put option on a zero-coupon bond.
The holder has the right to sell at time `expiry` a ZCB maturing at
`bond_maturity` for `strike`.
"""
struct ZCBPut{T <: Real, S <: Real, K <: Real} <: AbstractContract
    expiry::T
    bond_maturity::S
    strike::K
end

"""
    Cap(strike, frequency, maturity)

An interest rate cap — a portfolio of caplets that pay
`max(L(Tᵢ₋₁,Tᵢ) - strike, 0) · τ` at each payment date `Tᵢ`,
where `L` is the simply-compounded forward rate and `τ = 1/frequency`.

The first caplet resets at time `τ` (the first period's rate is known).
"""
struct Cap{K <: Real, F, M <: Real} <: AbstractContract
    strike::K
    frequency::F
    maturity::M
end

"""
    Floor(strike, frequency, maturity)

An interest rate floor — a portfolio of floorlets that pay
`max(strike - L(Tᵢ₋₁,Tᵢ), 0) · τ` at each payment date `Tᵢ`.
"""
struct Floor{K <: Real, F, M <: Real} <: AbstractContract
    strike::K
    frequency::F
    maturity::M
end

"""
    Swaption(expiry, swap_maturity, strike, frequency; payer=true)

A European swaption — the right to enter an interest rate swap at `expiry`.
The underlying swap has payment dates from `expiry + 1/frequency` to
`swap_maturity`, paying a fixed rate `strike`.

- `payer=true` (default): right to pay fixed, receive floating
- `payer=false`: right to receive fixed, pay floating
"""
struct Swaption{T <: Real, M <: Real, K <: Real, F} <: AbstractContract
    expiry::T
    swap_maturity::M
    strike::K
    frequency::F
    payer::Bool
end

function Swaption(expiry, swap_maturity, strike, frequency; payer = true)
    return Swaption(expiry, swap_maturity, strike, frequency, payer)
end

import ..FinanceCore: maturity
maturity(c::ZCBCall) = c.bond_maturity
maturity(c::ZCBPut) = c.bond_maturity
maturity(c::Cap) = c.maturity
maturity(c::Floor) = c.maturity
maturity(c::Swaption) = c.swap_maturity

end

"""
Forward(time,instrument)

The instrument is relative to the Forward time.
e.g. if you have a `Forward(1.0, Cashflow(1.0, 3.0))` then the instrument is a cashflow that pays 1.0 at time 4.0
"""
struct Forward{T <: FinanceCore.Timepoint, I <: FinanceCore.AbstractContract} <: FinanceCore.AbstractContract
    time::T
    instrument::I
end


"""
    cashflows_timepoints(contracts)
    cashflows_timepoints(quotes)

Create a matrix of cashflows and a vector of timepoints for a collection of quotes or contracts. Timepoints need not be spaced evenly.

This is used when constructing SmithWilson yield curves.

# Arguments
- `contracts` or `quotes`: A collection of `<:AbstractContract`s or `Quotes`.

# Returns
- A tuple `(m, times)` where `m` is a matrix of cashflows and `times` is a vector of timepoints.

# Examples
```julia-repl
julia> FinanceModels.cashflows_timepoints(ParYield.([0.04,0.02,0.04],[1,4,4]))
([0.02 0.01 0.02; 1.02 0.01 0.02; … ; 0.0 0.01 0.02; 0.0 1.01 1.02], [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
```
"""
function cashflows_timepoints(qs)
    cfs = map(q -> collect(q), qs)
    times = map(cfs) do cf
        map(c -> c.time, cf)
    end |> Iterators.flatten |> unique |> sort!

    m = zeros(length(times), length(qs))

    for t in 1:length(times)
        for q in 1:length(qs)
            for c in 1:length(cfs[q])
                if times[t] == cfs[q][c].time
                    m[t, q] += cfs[q][c].amount
                end
            end
        end
    end
    m
    return m, times
end

function cashflows_timepoints(qs::Vector{Q}) where {Q <: Quote}
    return cashflows_timepoints([q.instrument for q in qs])
end

"""
    InterestRateSwap(curve, tenor; model_key="OIS")

A convenience method for creating an interest rate swap given a curve and a tenor via a `Composite` contract consisting of receiving a [fixed bond](@ref Bond.Fixed) and paying (i.e. the negative of) a [floating bond](@ref Bond.Floating).

The notional is a unit (1.0) amount and assumed to settle four times per period.


A [`Projection`](@ref), with an indexable `model_key` is still needed to project a swap. See examples below for what this looks like.

# Examples

```julia-repl

julia> curve = Yield.Constant(0.05);

julia> swap = InterestRateSwap(curve,10);

julia> Projection(swap,Dict("OIS" => curve),CashflowProjection()) |> collect
80-element Vector{Cashflow{Float64, Float64}}:
Cashflow{Float64, Float64}(0.012272234429039353, 0.25)
Cashflow{Float64, Float64}(0.012272234429039353, 0.5)
⋮
Cashflow{Float64, Float64}(-0.012272234429039353, 9.75)
Cashflow{Float64, Float64}(-1.0122722344290391, 10.0)

```

"""
function InterestRateSwap(curve, tenor; model_key = "OIS")
    fixed_rate = par(curve, tenor; frequency = 4)
    fixed_leg = Bond.Fixed(rate(fixed_rate), Periodic(4), tenor)
    float_leg = Bond.Floating(0.0, Periodic(4), tenor, model_key) |> Map(-)
    return Composite(fixed_leg, float_leg)
end
