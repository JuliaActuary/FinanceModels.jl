# ConstructionBase protocol for the knot curves, whose derived caches must follow their inputs:
#  * the properties are (spline, rates, tenors, extrapolation); `getproperties` excludes caches;
#  * `setproperties` rebuilds through `reconstruct`, whose keywords are exactly those properties:
#    a patch to a cache (or any unknown key) is a MethodError rather than silently ignored, and
#    `setproperties(c, getproperties(c)) == c`;
#  * `constructorof` discards the cache arguments and rebuilds, so
#    `constructorof(typeof(c))(getfields(c)...) == c`.
# Every `@set`/`setall`/`fit` reconstruction therefore revalidates and rebuilds the caches.
Accessors.ConstructionBase.getproperties(c::Yield.AbstractInterpolatedZeroCurve) =
    (spline = c.spline, rates = c.rates, tenors = c.tenors, extrapolation = c.extrapolation)
Accessors.ConstructionBase.setproperties(c::Yield.AbstractInterpolatedZeroCurve, patch::NamedTuple) =
    Yield.reconstruct(c; patch...)
__rebuild_from_fields(spline, rates, tenors, extrapolation, _...) =
    Yield.ZeroRateCurve(rates, tenors, spline; extrapolation)
Accessors.ConstructionBase.constructorof(::Type{<:Yield.Spline}) = __rebuild_from_fields
Accessors.ConstructionBase.constructorof(::Type{<:Yield.MonotoneConvex}) = __rebuild_from_fields

# `FlatForwardAt` stores its forward continuously compounded and has no bare-number method,
# so Accessors rebuilds it from the stored value through `Continuous`.
Accessors.ConstructionBase.constructorof(::Type{<:Yield.FlatForwardAt}) =
    f -> Yield.FlatForwardAt(Continuous(f))
