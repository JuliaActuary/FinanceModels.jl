function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

    let 
        rates = [0.057, 0.0755, 0.0837, 0.078, 0.1084, 0.0702, 0.1167]
        tenors = [1,3,5,10,15,20,25]
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
    let 
        rates = [0.057, 0.0755, 0.0837, 0.078, 0.1084, 0.0702, 0.1167]
        tenors = [1.0,3,5,10,15,20,25]
        Yields.Par(rates,tenors)
        Yields.CMT(rates,tenors)
        Yields.Forward(rates,tenors)
        Yields.OIS(rates,tenors)
        Yields.Zero(rates,tenors)
        Yields.Zero(NelsonSiegel(), rates,tenors)
        Yields.Zero(NelsonSiegelSvensson(), rates,tenors)
                
    end
    
end