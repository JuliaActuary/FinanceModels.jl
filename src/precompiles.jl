
# created with the help of SnoopCompile.jl
@precompile_setup begin

    
    # 2021-03-31 rates from Treasury.gov
    rates =[0.01, 0.01, 0.03, 0.05, 0.07, 0.16, 0.35, 0.92, 1.40, 1.74, 2.31, 2.41] ./ 100
    tenors = [1/12, 2/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30]

    @precompile_all_calls begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        Yields.Par(rates,tenors)
        Yields.CMT(rates,tenors)
        Yields.Forward(rates,tenors)
        Yields.OIS(rates,tenors)
        Yields.Zero(NelsonSiegel(), rates,tenors)
        Yields.Zero(NelsonSiegelSvensson(), rates,tenors)
        Yields.Zero(rates,tenors)
        c = Yields.Zero(rates,tenors)
        Yields.zero(c,10)
        Yields.par(c,10)
        Yields.forward(c,5,6)
    end
        
end