@testset "Partial Trace" begin

    rng = MersenneTwister(42)

    ks = 2:6
    ds = 1:3

    for k in ks
        for d in ds
            for k_A in 1:(k-1)
                basis_full = LtrAscCutoffBasis(k,d)
                basis_A = LtrAscCutoffBasis(k_A, d)
                basis_B = LtrAscCutoffBasis(k-k_A, d)

                D = length(basis_full)

                ψ = randn(D)
                ψ ./= norm(ψ)

                D_A = length(basis_A)
                O_A_1 = zeros(D_A, D_A)
                D_B = length(basis_B)
                partialtr!(O_A_1, ψ, basis_full, basis_A, basis_B)

                ρ = ψ * ψ'
                O_A_2 = zeros(D_A, D_A)
                partialtr!(O_A_2, ρ, basis_full, basis_A, basis_B)
                @test O_A_1 ≈ O_A_2
            end
        end
    end

end
