@testset "Test operators" begin

    #############################  Explicit Tunnel ############################
    @testset "tunnel operators" begin
        @testset "2-site, N=1, 1 -> 2" begin
            N = 1
            k = 2
            basis = LtrAscBasis(k, N)
            tun_corr = spzeros(2,2)
            tun_corr[1,2] = 1.

            tun_calc = tunnel_spmatrix(1, 2, basis)
            @test tun_corr == tun_calc
        end

        @testset "2-site, N=1, 2 -> 1" begin
            N = 1
            k = 2
            basis = LtrAscBasis(k, N)
            tun_corr = spzeros(2,2)
            tun_corr[2,1] = 1.

            tun_calc = tunnel_spmatrix(2, 1, basis)
            @test tun_corr == tun_calc
        end

        @testset "2-site, N=2, 1 -> 2" begin
            N = 2
            k = 2
            basis = LtrAscBasis(k, N)
            tun_corr = spzeros(3,3)
            tun_corr[1,2] = sqrt(2)
            tun_corr[2,3] = sqrt(2)

            tun_calc = tunnel_spmatrix(1, 2, basis)
            @test isapprox(tun_corr, tun_calc;atol=1e-15)
        end

        @testset "2-site, N=2, 2 -> 1" begin
            N = 2
            k = 2
            basis = LtrAscBasis(k, N)
            tun_corr = spzeros(3,3)
            tun_corr[2,1] = sqrt(2)
            tun_corr[3,2] = sqrt(2)

            tun_calc = tunnel_spmatrix(2, 1, basis)
            @test isapprox(tun_corr, tun_calc;atol=1e-15)
        end

        @testset "3-site, N=1, 1 -> 2" begin
            N = 1
            k = 3
            basis = LtrAscBasis(k, N)
            tun_corr = spzeros(3,3)
            tun_corr[2,3] = 1.

            tun_calc = tunnel_spmatrix(1, 2, basis)
            @test isapprox(tun_corr, tun_calc;atol=1e-15)
        end

        @testset "3-site, N=1, 2 -> 1" begin
            N = 1
            k = 3
            basis = LtrAscBasis(k, N)
            tun_corr = spzeros(3,3)
            tun_corr[3,2] = 1.

            tun_calc = tunnel_spmatrix(2, 1, basis)
            @test isapprox(tun_corr, tun_calc;atol=1e-15)
        end

        @testset "3-site, N=1, 2 -> 3" begin
            N = 1
            k = 3
            basis = LtrAscBasis(k, N)
            tun_corr = spzeros(3,3)
            tun_corr[1,2] = 1.

            tun_calc = tunnel_spmatrix(2, 3, basis)
            @test isapprox(tun_corr, tun_calc;atol=1e-15)
        end

        @testset "3-site, N=1, 3 -> 2" begin
            N = 1
            k = 3
            basis = LtrAscBasis(k, N)
            tun_corr = spzeros(3,3)
            tun_corr[2,1] = 1.

            tun_calc = tunnel_spmatrix(3, 2, basis)
            @test isapprox(tun_corr, tun_calc;atol=1e-15)
        end

        @testset "3-site, N=1, 1 -> 2" begin
            N = 2
            k = 3
            basis = LtrAscBasis(k, N)
            tun_corr = spzeros(6,6)
            tun_corr[2,4] = 1.
            tun_corr[3,5] = sqrt(2)
            tun_corr[5,6] = sqrt(2)

            tun_calc = tunnel_spmatrix(1, 2, basis)
            @test isapprox(tun_corr, tun_calc;atol=1e-15)
        end
    end

#############################  Explicit Number ###############################
    @testset "number operators" begin
        @testset "2-site, N=1, n_1" begin
            N = 1
            k = 2
            basis = LtrAscBasis(k, N)
            corr = spzeros(2,2)
            corr[2,2] = 1.

            calc = number_spmatrix(1, basis)
            @test corr == calc
        end

        @testset "2-site, N=1, n_2" begin
            N = 1
            k = 2
            basis = LtrAscBasis(k, N)
            corr = spzeros(2,2)
            corr[1,1] = 1.

            calc = number_spmatrix(2, basis)
            @test corr == calc
        end

        @testset "2-site, N=2, n_1" begin
            N = 2
            k = 2
            basis = LtrAscBasis(k, N)
            corr = spzeros(3,3)
            corr[2,2] = 1.
            corr[3,3] = 2.

            calc = number_spmatrix(1, basis)
            @test corr == calc
        end

        @testset "2-site, N=2, n_2" begin
            N = 2
            k = 2
            basis = LtrAscBasis(k, N)
            corr = spzeros(3,3)
            corr[1,1] = 2.
            corr[2,2] = 1.

            calc = number_spmatrix(2, basis)
            @test corr == calc
        end

        @testset "3-site, N=1, n_1" begin
            N = 1
            k = 3
            basis = LtrAscBasis(k, N)
            corr = spzeros(3,3)
            corr[3,3] = 1.

            calc = number_spmatrix(1, basis)
            @test corr == calc
        end

        @testset "3-site, N=1, n_2" begin
            N = 1
            k = 3
            basis = LtrAscBasis(k, N)
            corr = spzeros(3,3)
            corr[2,2] = 1.

            calc = number_spmatrix(2, basis)
            @test corr == calc
        end

        @testset "3-site, N=1, n_3" begin
            N = 1
            k = 3
            basis = LtrAscBasis(k, N)
            corr = spzeros(3,3)
            corr[1,1] = 1.

            calc = number_spmatrix(3, basis)
            @test corr == calc
        end

        @testset "3-site, N=2, n_1" begin
            N = 2
            k = 3
            basis = LtrAscBasis(k, N)
            corr = spzeros(6,6)
            corr[4,4] = 1.
            corr[5,5] = 1.
            corr[6,6] = 2.

            calc = number_spmatrix(1, basis)
            @test corr == calc
        end

        @testset "3-site, N=2, n_2" begin
            N = 2
            k = 3
            basis = LtrAscBasis(k, N)
            corr = spzeros(6,6)
            corr[2,2] = 1.
            corr[3,3] = 2.
            corr[5,5] = 1.

            calc = number_spmatrix(2, basis)
            @test corr == calc
        end

        @testset "3-site, N=2, n_3" begin
            N = 2
            k = 3
            basis = LtrAscBasis(k, N)
            corr = spzeros(6,6)
            corr[1,1] = 2.
            corr[2,2] = 1.
            corr[4,4] = 1.

            calc = number_spmatrix(3, basis)
            @test corr == calc
        end

    end

    #=
    @testset "Test operators with LtrAsc0NBasis" begin

        basis = LtrAscBasis0N(3, 4)

        tunnel_spmatrix(1,2, basis)
        number_spmatrix(1, basis)


    end
    =#


#############################  creation/annihilation ###########################
    @testset "Creation/Annihilation operators" begin

        @testset "a_i^dagger a_i = n_i" begin

            ks = 2:5
            Ns = 1:5

            for k in ks
                for N in Ns
                    basis1 = LtrAscBasis(k, N)
                    basis2 = LtrAscBasis(k, N-1)

                    for i in 1:k
                        n_i = number_spmatrix(i, basis1)
                        a_i_dagger = creation_spmatrix(i, basis2, basis1)
                        a_i = annihilation_spmatrix(i, basis1, basis2)

                        @test isapprox(n_i, a_i_dagger*a_i ;atol=1e-14)
                    end
                end
            end
        end

        @testset "a_i^dagger a_j = tunnel_i_j" begin

            ks = 2:5
            Ns = 1:5

            for k in ks
                for N in Ns
                    basis1 = LtrAscBasis(k, N)
                    basis2 = LtrAscBasis(k, N-1)

                    for i in 1:k
                        for j in 1:k
                            i==j && continue

                            tunnel = tunnel_spmatrix(j, i, basis1)
                            a_i_dagger = creation_spmatrix(i, basis2, basis1)
                            a_j = annihilation_spmatrix(j, basis1, basis2)

                            @test isapprox(tunnel, a_i_dagger*a_j ;atol=1e-14)
                        end
                    end
                end
            end

        end





    end

end
