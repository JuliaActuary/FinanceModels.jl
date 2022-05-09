@testset "NelsonSiegelSvensson" begin
    
    @testset "NelsonSiegel" begin

        # EURAAA_20191111
        @testset "EURAAA_20191111" begin
            euraaa_yields = [-0.602, -0.6059, -0.6096, -0.613, -0.6215, -0.6279, -0.6341, -0.6327, -0.6106, -0.5694, -0.5161, -0.456, -0.3932, -0.3305, -0.2698, -0.2123, -0.1091, 0.0159, 0.0818, 0.1601, 0.2524, 0.315]
            euraaa_maturities = [0.25, 0.333, 0.417, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 17, 20, 25, 30]
            ns_param = NelsonSiegel(0.6440719852424944, -1.2135915867261, -1.8281745438894712, 0.1)
            ns = est_ns_params(euraaa_yields, euraaa_maturities)
            @test ns = ns_param
        end
    end

    @testset "NelsonSiegelSvensson" begin

        # EURAAA_20191111
        @testset "EURAAA_20191111" begin
            euraaa_yields = [-0.602, -0.6059, -0.6096, -0.613, -0.6215, -0.6279, -0.6341, -0.6327, -0.6106, -0.5694, -0.5161, -0.456, -0.3932, -0.3305, -0.2698, -0.2123, -0.1091, 0.0159, 0.0818, 0.1601, 0.2524, 0.315]
            euraaa_maturities = [0.25, 0.333, 0.417, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 17, 20, 25, 30]
            nss_param = NelsonSiegelSvensson(0.6343710125821183, -1.225225531496011, -2.107914746999533, 0.2892748213968784, 3.0, 1.5)
            nss = est_nss_params(euraaa_yields, euraaa_maturities)
            @test nss = nss_param
        end
    end

end
