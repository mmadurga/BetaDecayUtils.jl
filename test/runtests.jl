using BetaDecayUtils
using Test

@testset "BetaDecayUtils.jl" begin
    @testset "square-well neutron width" begin
        E, V0, AM1, AM2 = 1.0, 50.0, 100.0, 1.0
        R = 1.4 * (AM1 + AM2)^(1/3)

        # Near threshold, neutron penetrability scales as (kR)^(2ℓ+1).
        for ℓ in 0:3
            P1 = nPenetrability(1e-8, AM1, AM2, ℓ; R=R)
            P2 = nPenetrability(4e-8, AM1, AM2, ℓ; R=R)
            @test P2/P1 ≈ 2^(2ℓ + 1) rtol=1e-5
        end

        μ = AM1*AM2/(AM1+AM2)*931.502
        for ℓ in 0:7
            P = nPenetrability(E, AM1, AM2, ℓ; R=R)
            factor = ℓ == 0 ? 1 : (2ℓ-1)/(2ℓ+1)
            expected = 2*P*197.326^2/(μ*R^2)*factor
            @test Gamma(E, AM1, AM2, ℓ; R=R) ≈ expected rtol=1e-12
        end

        # Regression values from the Nucracker calculator.
        AAg, RAg = 131.0, 1.4*132^(1/3)
        @test Gamma(1.0, AAg, 1.0, 0; R=RAg) ≈ 2.57517 rtol=5e-5
        @test Gamma(1.0, AAg, 1.0, 1; R=RAg) ≈ 0.606863 rtol=5e-5
        @test Gamma(1.0, AAg, 1.0, 2; R=RAg) ≈ 0.407737 rtol=5e-5
        ASn, RSn = 132.0, 1.4*133^(1/3)
        @test Gamma(1.0, ASn, 1.0, 1; R=RSn) ≈ 0.606223 rtol=5e-5

        # Ratios to the Wigner limit, Gamma_W = 3P*(ħc)^2/(μ R^2).
        for (ℓ, ratio) in ((1, 2/9), (2, 2/5), (3, 10/21))
            P = nPenetrability(E, AM1, AM2, ℓ; R=R)
            wigner = 3*P*197.326^2/(μ*R^2)
            @test Gamma(E, AM1, AM2, ℓ; R=R)/wigner ≈ ratio rtol=1e-12
        end

        @test Gamma(0.0, AM1, AM2, 0; R=R) == 0
        @test Gamma(0.0, AM1, AM2, 1; R=R) == 0
        @test_deprecated Gamma(E, V0, AM1, AM2, 1; R=R)
        @test_deprecated Gamma(E, V0, AM1, AM2, 1, R)
        @test_throws DomainError Gamma(-1.0, AM1, AM2, 1; R=R)
        @test_throws DomainError Gamma(E, 0.0, AM2, 1; R=R)
        @test_throws DomainError Gamma(E, AM1, AM2, 1; R=0.0)
        @test_throws ArgumentError Gamma(E, AM1, AM2, 0.5; R=R)
        @test_throws ArgumentError Gamma(E, AM1, AM2, 8; R=R)
    end
end
