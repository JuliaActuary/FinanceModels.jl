@testset "NelsonSiegelSvensson" begin

    # Per this technical note, the target/source rates should be interpreted as continuous rates: 
    # https://www.ecb.europa.eu/stats/financial_markets_and_interest_rates/euro_area_yield_curves/html/technical_notes.pdf
    
    @testset "NelsonSiegel" begin

        # EURAAA_20191111 at https://www.ecb.europa.eu/stats/financial_markets_and_interest_rates/euro_area_yield_curves/html/index.en.html
        @testset "EURAAA_20191111" begin
            euraaa_zeros = Continuous.([-0.602009,-0.612954,-0.621543,-0.627864,-0.632655,-0.610565,-0.569424,-0.516078,-0.455969,-0.39315,-0.33047,-0.269814,-0.21234,-0.158674,-0.109075,-0.063552,-0.021963,0.015929,0.050407,0.081771,0.110319,0.136335,0.160083,0.181804,0.201715,0.220009,0.23686,0.252419,0.26682,0.280182,0.292608,0.304191,0.31501] ./ 100)
            euraaa_pars = Continuous.([-0.601861,-0.612808,-0.621403,-0.627731,-0.63251,-0.610271,-0.568825,-0.515067,-0.454526,-0.391338,-0.328412,-0.26767,-0.210277,-0.156851,-0.107632,-0.062605,-0.021601,0.015642,0.049427,0.080073,0.107892,0.133178,0.156206,0.17722,0.196444,0.214073,0.230281,0.245221,0.259027,0.271818,0.283696,0.294753,0.305069] ./ 100)
            euraaa_maturities = vcat([0.25, 0.5, 0.75], 1:30)
            c_param = Yields.NelsonSiegelCurve(0.6202603126029489 /100, -1.1621281759833935 /100, -1.930016080035979 /100, 3.0)
            c = Yields.Zero(NelsonSiegel(),euraaa_zeros, euraaa_maturities)
            # @testset "parameter: $param" for param in [:β₀, :β₁, :β₂, :τ₁]
            #     @test getfield(c, param) ≈ getfield(c_param, param)
            # end

            @testset "zero rates: $t" for (t, r) in zip(euraaa_maturities, euraaa_zeros)
                @test Yields.zero(c, t) ≈ r   atol = 0.001
            end

            @testset "par rates: $t" for (t, r) in zip(euraaa_maturities, euraaa_pars)
                @test Yields.zero(c, t) ≈ r   atol = 0.001
            end

            @test discount(c,0) == 1.0
            @test discount(c,10) ≈ 1 / accumulation(c,10)
        end

        # Nelson-Siegel-Svensson package example at https://nelson-siegel-svensson.readthedocs.io/en/latest/usage.html
        @testset "pack" begin
            pack_yields = Continuous.([0.01, 0.011, 0.013, 0.016, 0.019, 0.021, 0.026, 0.03, 0.035, 0.037, 0.038, 0.04])
            pack_maturities = [10e-5, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
            c_param = Yields.NelsonSiegelCurve(0.04495841387198023, -0.03537510042719209, 0.0031561222355027227, 5.0)
            c = Yields.Zero(NelsonSiegel(),pack_yields, pack_maturities)

            @testset "zero rates: $t" for (t, r) in zip(pack_maturities, pack_yields)
                @test Yields.zero(c, t) ≈ r   atol = 0.003
            end
        end
    end

    @testset "NelsonSiegelSvensson" begin

        @testset "EURAAA_20191111" begin
            # EURAAA_20191111 at https://www.ecb.europa.eu/stats/financial_markets_and_interest_rates/euro_area_yield_curves/html/index.en.html

            euraaa_zeros = Continuous.([-0.602009,-0.612954,-0.621543,-0.627864,-0.632655,-0.610565,-0.569424,-0.516078,-0.455969,-0.39315,-0.33047,-0.269814,-0.21234,-0.158674,-0.109075,-0.063552,-0.021963,0.015929,0.050407,0.081771,0.110319,0.136335,0.160083,0.181804,0.201715,0.220009,0.23686,0.252419,0.26682,0.280182,0.292608,0.304191,0.31501] ./ 100)
            euraaa_pars = Continuous.([-0.601861,-0.612808,-0.621403,-0.627731,-0.63251,-0.610271,-0.568825,-0.515067,-0.454526,-0.391338,-0.328412,-0.26767,-0.210277,-0.156851,-0.107632,-0.062605,-0.021601,0.015642,0.049427,0.080073,0.107892,0.133178,0.156206,0.17722,0.196444,0.214073,0.230281,0.245221,0.259027,0.271818,0.283696,0.294753,0.305069] ./ 100)
            euraaa_maturities = vcat([0.25, 0.5, 0.75], 1:30)
            c_zero = Yields.Zero(NelsonSiegelSvensson(),euraaa_zeros, euraaa_maturities)
            c_par = Yields.Par(NelsonSiegelSvensson(),euraaa_pars, euraaa_maturities)

            @testset "par and zero constructors" for c in [c_zero,c_par]
                @testset "zero rates: $t" for (t, r) in zip(euraaa_maturities, euraaa_zeros)
                    @test Yields.zero(c, t) ≈ r   atol = 0.0001
                end

                @testset "par rates: $t" for (t, r) in zip(euraaa_maturities, euraaa_pars)
                    # are the target rates on the ECB site continuous rates or periodic/bond-equivalent?
                    @test Yields.par(c, t) ≈ r   atol = 0.0001
                end
            end

        end
        @testset "EURAAA_20191111 w parms given" begin
            # EURAAA_20191111 at https://www.ecb.europa.eu/stats/financial_markets_and_interest_rates/euro_area_yield_curves/html/index.en.html
            euraaa_zeros = Continuous.([-0.602009,-0.612954,-0.621543,-0.627864,-0.632655,-0.610565,-0.569424,-0.516078,-0.455969,-0.39315,-0.33047,-0.269814,-0.21234,-0.158674,-0.109075,-0.063552,-0.021963,0.015929,0.050407,0.081771,0.110319,0.136335,0.160083,0.181804,0.201715,0.220009,0.23686,0.252419,0.26682,0.280182,0.292608,0.304191,0.31501] ./ 100)
            euraaa_pars = Continuous.([-0.601861,-0.612808,-0.621403,-0.627731,-0.63251,-0.610271,-0.568825,-0.515067,-0.454526,-0.391338,-0.328412,-0.26767,-0.210277,-0.156851,-0.107632,-0.062605,-0.021601,0.015642,0.049427,0.080073,0.107892,0.133178,0.156206,0.17722,0.196444,0.214073,0.230281,0.245221,0.259027,0.271818,0.283696,0.294753,0.305069] ./ 100)
            euraaa_maturities = vcat([0.25, 0.5, 0.75], 1:30)
            c = Yields.NelsonSiegelSvenssonCurve(0.629440 / 100, -1.218082 /100, 	12.114098 /100, -14.181117 /100, 2.435976, 2.536963)

            @testset "zero rates: $t" for (t, r) in zip(euraaa_maturities, euraaa_zeros)
                @test Yields.zero(c, t) ≈ r   atol = 0.0001
            end

            @testset "par rates: $t" for (t, r) in zip(euraaa_maturities, euraaa_pars)
                # are the target rates on the ECB site continuous rates or periodic/bond-equivalent?
                @test Yields.par(c, t) ≈ r atol = 0.0001
            end

            @test Yields.zero(c,30) ≈ last(euraaa_zeros) atol = 0.0001 
            @test discount(c,0) == 1.0
            @test discount(c,10) ≈ 1 / accumulation(c,10)
        end
    end

end
