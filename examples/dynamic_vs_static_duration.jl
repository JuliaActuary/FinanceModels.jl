#=
Deferred Annuity: Dynamic vs Static Duration (Two-Curve, Rate Reset)
====================================================================
Models a 5-year fixed deferred annuity with rate-sensitive (dynamic) lapse
behavior and compares effective duration under two assumptions:

  - **Dynamic duration**: Lapse rate at renewal responds to rate changes via AD
  - **Static duration**: Cashflows frozen at base rates; only discounting changes

Two-curve discounting:
  - Base curve = US Treasury (CMT)
  - Credit curve = flat 120bps spread over Treasury
  - Total discount at time t: D(t) = D_base(t) · D_credit(t)
  - Produces separate base (IR) and credit (CS) key rate durations

Product:
  - 5-year surrender charge period, 10-year projection horizon
  - Credited rate = 5yr CMT par + 60bps (fixed at issue, continuously compounded)
  - During SC period (years 1–5): 2% p.a. fixed lapse (SC suppresses optionality)
  - At year-5 renewal: dynamic lapse via logistic function of (market − credited),
    calibrated to ~40% at base rates
    (logistic form from https://modernfinancialmodeling.com/autodiff_alm)
  - **Credited rate resets at renewal** to new 5yr market rate + 60bps
    for policyholders who stay — the post-renewal liability reprices to market
  - Post-renewal (years 6–10): 2% p.a. fixed lapse
  - All remaining account value paid out at year 10

Run:  julia examples/dynamic_vs_static_duration.jl
=#

# ── Environment setup ────────────────────────────────────────────────────────
using Pkg
Pkg.activate(; temp=true)
Pkg.develop(path=dirname(@__DIR__))   # use local FinanceModels
Pkg.add("ActuaryUtilities")

using FinanceModels
using ActuaryUtilities
using Printf

# ── 1. CMT Yield Curve ──────────────────────────────────────────────────────

# Representative US Treasury CMT rates
cmt_rates = [5.25, 5.30, 5.35, 5.40, 5.50, 5.60, 5.65, 5.70, 5.75, 5.80] ./ 100
cmt_mats  = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0]

qs = CMTYield.(cmt_rates, cmt_mats)
bootstrapped = fit(Spline.Linear(), qs, Fit.Bootstrap())

# ZeroRateCurve for AD-based sensitivities
tenors = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
zrc_base = ZeroRateCurve(bootstrapped, tenors)

# Credit spread curve: flat 120bps over treasury (same tenors required)
credit_spread = 0.0120   # 120bps
zrc_credit = ZeroRateCurve(fill(credit_spread, length(tenors)), tenors)

println("── Base Curve (continuously compounded zero rates) ──")
for (t, r) in zip(zrc_base.tenors, zrc_base.rates)
    @printf("  %5.1fy: %.4f%%  (tsy) + %.0fbps (spread) = %.4f%%\n",
            t, r * 100, credit_spread * 10_000, (r + credit_spread) * 100)
end

# ── 2. Annuity Parameters ───────────────────────────────────────────────────

new_money_spread = 0.0060                    # 60bps new money spread over par
par_5yr = par(bootstrapped, 5.0)
credited_rate = rate(Continuous(par_5yr)) + new_money_spread  # continuously compounded

fixed_lapse_annual = 0.02                    # 2% p.a. during SC & post-renewal
fixed_lapse_monthly = -expm1(-fixed_lapse_annual / 12)  # exact monthly equivalent

renewal_month = 60                           # year 5
n_months = 120                               # 10-year projection
dt = 1 / 12

@printf("\n── Annuity Parameters ──\n")
@printf("  5yr par rate:    %.4f%% (semi-annual)\n", rate(par_5yr) * 100)
@printf("  Credited rate:   %.4f%% (cc) = par(cc) + %.0fbps\n",
        credited_rate * 100, new_money_spread * 10_000)
@printf("  Rate reset:      at renewal, new 5yr fwd rate + %.0fbps\n",
        new_money_spread * 10_000)
@printf("  Discount spread: %.0fbps over treasury\n", credit_spread * 10_000)
@printf("  Fixed lapse:     %.1f%% p.a. (SC period & post-renewal)\n",
        fixed_lapse_annual * 100)
@printf("  Projection:      %d years, monthly\n", n_months ÷ 12)

# ── 3. Dynamic Lapse Function (renewal only) ────────────────────────────────

# Logistic surrender rate from modernfinancialmodeling.com/autodiff_alm
#   sr(Δr) = 1 / (1 + exp(a − b·Δr))
# where Δr = market_rate − credited_rate.
#
# Calibrate midpoint `a` so sr(0) = 40% (base renewal lapse):
#   1/(1+exp(a)) = 0.40  →  a = log(1/0.40 − 1) = log(1.5)

target_renewal_lapse = 0.40
midpoint = log(1 / target_renewal_lapse - 1)
sensitivity = 60.0

surrender_rate(rate_diff) = 1 / (1 + exp(midpoint - sensitivity * rate_diff))

@printf("  Renewal lapse:   %.0f%% base (dynamic, logistic)\n", target_renewal_lapse * 100)
@printf("  Logistic params: midpoint=%.4f, sensitivity=%.1f\n", midpoint, sensitivity)
println()

# ── 4. Liability PV (two-curve, dynamic renewal lapse, rate reset) ──────────
#
# base_curve  = treasury (drives market rate for lapse + part of discounting)
# credit_curve = spread  (only affects discounting)
# total discount at t:  base(t) · credit(t)
#
# At renewal, the credited rate resets to the 5yr forward rate + new money
# spread.  This means the post-renewal AV growth tracks market rates,
# collapsing the duration of the tail.

function liability_pv(base_curve, credit_curve, original_credited_rate,
                      new_money_spread, fixed_lapse_monthly,
                      renewal_month, n_months, dt)
    # Initialize with correct numeric type for ForwardDiff compatibility
    df1 = discount(base_curve, dt) * discount(credit_curve, dt)
    av = one(df1)
    pv = zero(df1)

    # Credited rate starts at original; resets at renewal
    cr = original_credited_rate * one(df1)

    for m in 1:n_months
        t = m * dt

        # Credit interest at current credited rate
        av *= exp(cr * dt)

        # Total discount factor at time t
        df_total = discount(base_curve, t) * discount(credit_curve, t)

        if m == renewal_month
            # ── Dynamic lapse at renewal ──
            # Market rate: 5yr continuous forward from year 5 (treasury only)
            df_base_t = discount(base_curve, t)
            df_base_fwd = discount(base_curve, t + 5.0)
            mkt_rate = -log(df_base_fwd / df_base_t) / 5.0

            rate_diff = mkt_rate - original_credited_rate
            sr = surrender_rate(rate_diff)

            surrendered = av * sr
            av -= surrendered
            pv += surrendered * df_total

            # ── Reset credited rate for remaining policyholders ──
            cr = mkt_rate + new_money_spread
        else
            # ── Fixed lapse (SC suppresses optionality) ──
            surrendered = av * fixed_lapse_monthly
            av -= surrendered
            pv += surrendered * df_total
        end
    end

    # Terminal: all remaining AV paid out at year 10
    pv += av * discount(base_curve, n_months * dt) * discount(credit_curve, n_months * dt)
    return pv
end

# ── 5. Project Fixed Cashflows (for static duration) ────────────────────────

function project_cashflows(curve, original_credited_rate, new_money_spread,
                           fixed_lapse_monthly, renewal_month, n_months, dt)
    av = 1.0
    cfs = Float64[]
    times = Float64[]

    cr = original_credited_rate

    for m in 1:n_months
        t = m * dt
        av *= exp(cr * dt)

        if m == renewal_month
            df_t = discount(curve, t)
            df_fwd = discount(curve, t + 5.0)
            mkt_rate = -log(df_fwd / df_t) / 5.0

            rate_diff = mkt_rate - original_credited_rate
            sr = surrender_rate(rate_diff)

            surrendered = av * sr
            av -= surrendered

            # Reset credited rate
            cr = mkt_rate + new_money_spread
        else
            surrendered = av * fixed_lapse_monthly
            av -= surrendered
        end

        if surrendered > 1e-12
            push!(cfs, surrendered)
            push!(times, t)
        end
    end

    # Terminal payout
    push!(cfs, av)
    push!(times, n_months * dt)

    return cfs, times
end

# ── 6. Compute Sensitivities ────────────────────────────────────────────────

println("── Computing sensitivities... ──\n")

# Dynamic: AD flows through renewal lapse AND credited rate reset
dynamic = sensitivities(zrc_base, zrc_credit) do base_curve, credit_curve
    liability_pv(base_curve, credit_curve, credited_rate,
                 new_money_spread, fixed_lapse_monthly,
                 renewal_month, n_months, dt)
end

# Static: cashflows frozen at base curve — only discounting changes
static_cfs, static_times = project_cashflows(
    bootstrapped, credited_rate, new_money_spread,
    fixed_lapse_monthly, renewal_month, n_months, dt)
static = sensitivities(zrc_base, zrc_credit, static_cfs, static_times)

# ── 7. Results ───────────────────────────────────────────────────────────────

# Verify base-case renewal
base_df5 = discount(bootstrapped, 5.0)
base_df10 = discount(bootstrapped, 10.0)
base_mkt = -log(base_df10 / base_df5) / 5.0
base_sr = surrender_rate(base_mkt - credited_rate)
reset_rate = base_mkt + new_money_spread

println("══════════════════════════════════════════════════════════════════")
println("  DYNAMIC vs STATIC DURATION (Two-Curve, Rate Reset at Renewal)")
println("══════════════════════════════════════════════════════════════════")

@printf("\n  Base renewal market rate:  %.4f%% (5yr tsy fwd from yr 5)\n", base_mkt * 100)
@printf("  Base renewal lapse rate:  %.2f%%\n", base_sr * 100)
@printf("  Reset credited rate:      %.4f%% (cc) = fwd + %.0fbps\n",
        reset_rate * 100, new_money_spread * 10_000)

@printf("\n  Liability PV (dynamic): %.6f\n", dynamic.value)
@printf("  Liability PV (static):  %.6f\n", static.value)

# ── Base (Treasury) Key Rate Durations ──

println("\n  Base (Treasury) Key Rate Durations:")
println("  ─────────────────────────────────────────────────────")
@printf("  %8s  %12s  %12s  %12s\n", "Tenor", "Dynamic", "Static", "Δ")
println("  ─────────────────────────────────────────────────────")
for (i, t) in enumerate(zrc_base.tenors)
    d = dynamic.base_durations[i]
    s = static.base_durations[i]
    @printf("  %6.1fy    %10.4f    %10.4f    %10.4f\n", t, d, s, d - s)
end
println("  ─────────────────────────────────────────────────────")

total_base_dyn = sum(dynamic.base_durations)
total_base_sta = sum(static.base_durations)
@printf("  %8s  %10.4f    %10.4f    %10.4f\n",
        "TOTAL", total_base_dyn, total_base_sta, total_base_dyn - total_base_sta)

# ── Credit (Spread) Key Rate Durations ──

println("\n  Credit (Spread) Key Rate Durations:")
println("  ─────────────────────────────────────────────────────")
@printf("  %8s  %12s  %12s  %12s\n", "Tenor", "Dynamic", "Static", "Δ")
println("  ─────────────────────────────────────────────────────")
for (i, t) in enumerate(zrc_base.tenors)
    d = dynamic.credit_durations[i]
    s = static.credit_durations[i]
    @printf("  %6.1fy    %10.4f    %10.4f    %10.4f\n", t, d, s, d - s)
end
println("  ─────────────────────────────────────────────────────")

total_credit_dyn = sum(dynamic.credit_durations)
total_credit_sta = sum(static.credit_durations)
@printf("  %8s  %10.4f    %10.4f    %10.4f\n",
        "TOTAL", total_credit_dyn, total_credit_sta, total_credit_dyn - total_credit_sta)

# ── Summary ──

println("\n  ─────────────────────────────────────────────────────")
println("  Summary:                Dynamic     Static        Δ")
println("  ─────────────────────────────────────────────────────")
@printf("  Base (IR) duration:   %10.4f  %10.4f  %10.4f\n",
        total_base_dyn, total_base_sta, total_base_dyn - total_base_sta)
@printf("  Credit (CS) duration: %10.4f  %10.4f  %10.4f\n",
        total_credit_dyn, total_credit_sta, total_credit_dyn - total_credit_sta)
@printf("  Total eff. duration:  %10.4f  %10.4f  %10.4f\n",
        total_base_dyn + total_credit_dyn,
        total_base_sta + total_credit_sta,
        (total_base_dyn + total_credit_dyn) - (total_base_sta + total_credit_sta))
println("  ─────────────────────────────────────────────────────")

println("\n  Interpretation:")
@printf("  The credited rate resets at renewal → the post-renewal liability\n")
@printf("  reprices to market. The tail (years 6–10) has near-zero IR\n")
@printf("  duration because AV growth tracks the curve.\n\n")

base_gap = total_base_dyn - total_base_sta
@printf("  IR duration: dynamic=%.2f, static=%.2f, gap=%.4f\n",
        total_base_dyn, total_base_sta, base_gap)

if abs(base_gap) < 0.05
    println("  → Dynamic ≈ Static: the rate reset neutralizes post-renewal")
    println("    duration, validating a 5yr asset match + roll strategy.")
elseif base_gap < 0
    println("  → Dynamic < Static: renewal lapse optionality still shortens")
    println("    IR duration, but the magnitude is small vs the no-reset model.")
else
    println("  → Dynamic > Static: rate reset plus lapse optionality adds")
    println("    some IR sensitivity at the renewal point.")
end

@printf("\n  CS duration: %.2f (unaffected by lapse or rate reset —\n", total_credit_dyn)
println("    spread discounting applies uniformly to all cashflows).")

println()
println("  A 5yr asset match + reinvestment at renewal is well-supported")
println("  by this model: the effective IR duration is concentrated in")
println("  the first 5 years, with minimal tail risk after repricing.\n")
