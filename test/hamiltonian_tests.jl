@testset "Hamiltonian creation" begin

    @testset "create tunnel matrices dimer N=1" begin
        N = 1
        k = 2

        graph = dimer_graph()

        basis = LtrAscBasis(k, N)

        tun_corr = spzeros(2,2)
        tun_corr[1,2] = 1.
        tun_corr[2,1] = 1.

        tun_calc = tunnel_spmatrices(graph, basis)[1]

        @test tun_corr == tun_calc
    end

    @testset "hamiltonia dimer N=1" begin

        N = 1
        k = 2

        # Define the correct Hamiltonian
        function f_H_correct(J,U)
            H = Matrix{Float64}(I,2,2)
            H[1,1] = 0
            H[1,2] = -J
            H[2,1] = -J
            H[2,2] = 0

            H = 0.5 * H
            return H
        end

        # Define the values to test
        J_arr = 1:1:10

        # Prepare the system
        basis = LtrAscBasis(k, N)
        graph = dimer_graph()
        bhh = BoseHubbardHamiltonian(graph, basis)

        for J in J_arr
            H_corr = f_H_correct(J, 0)
            H_calc = matrix(bhh, [J],[0,0])

            @test H_corr == H_calc
        end
    end


    @testset "hamiltonian dimer N=2" begin

        N = 2
        k = 2

        # Define the correct Hamiltonian
        function f_H_correct(J,U)
            H = Matrix{Float64}(I,3,3)
            H[1,1] = 2*U
            H[1,2] = -J*sqrt(2)
            H[1,3] = 0
            H[2,1] = -J*sqrt(2)
            H[2,2] = 0
            H[2,3] = -J*sqrt(2)
            H[3,1] = 0
            H[3,2] = -J*sqrt(2)
            H[3,3] = 2*U

            H = 0.5 * H
            return H
        end

        # Define values to test
        U_arr = 0:2:10
        J_arr = 0:2:10

        # Prepare the system
        basis = LtrAscBasis(k, N)
        graph = dimer_graph()
        bhh = BoseHubbardHamiltonian(graph, basis)

        for U in U_arr
            for J in J_arr
                H_corr = f_H_correct(J, U)
                H_calc = matrix(bhh, [J],[U,U])

                @test H_corr == H_calc
            end
        end
    end

    @testset "hamiltonian dimer N=3" begin

        N = 3
        k = 2

        function f_H_correct(J,U)
            H = Matrix{Float64}(I, 4, 4)
            H[1,1] = 3*U
            H[1,2] = -0.5*J*sqrt(3)
            H[1,3] = 0
            H[1,4] = 0
            H[2,1] = -0.5*J*sqrt(3)
            H[2,2] = U
            H[2,3] = -J
            H[2,4] = 0
            H[3,1] = 0
            H[3,2] = -J
            H[3,3] = U
            H[3,4] = -0.5*J*sqrt(3)
            H[4,1] = 0
            H[4,2] = 0
            H[4,3] = -0.5*J*sqrt(3)
            H[4,4] = 3*U

            return H
        end

        # Define test values
        U_arr = 0:2:10
        J_arr = 0:2:10

        # Prepare the system
        basis = LtrAscBasis(k, N)
        graph = dimer_graph()
        bhh = BoseHubbardHamiltonian(graph, basis)

        for U in U_arr
            for J in J_arr
                H_corr = f_H_correct(J, U)
                H_calc = matrix(bhh, [J],[U,U])

                @test H_corr == H_calc
            end
        end
    end




    @testset "hamiltonian trimer N=1" begin

        N = 1
        k = 3

        function f_H_correct(J,eps)
            H = Matrix{Float64}(I, 3, 3)
            H[1,1] = eps[3]
            H[1,2] = -J[2]
            H[1,3] = 0
            H[2,1] = -J[2]
            H[2,2] = eps[2]
            H[2,3] = -J[1]
            H[3,1] = 0
            H[3,2] = -J[1]
            H[3,3] = eps[1]

            H *= 0.5

            return H
        end

        # Create test values
        J_arr = -5:1:5
        J_arr = Iterators.product(J_arr, J_arr)
        eps_arr = -0.25:0.25:0.25
        eps_arr = Iterators.product(eps_arr, eps_arr, eps_arr)

        # Prepare the system
        basis = LtrAscBasis(k, N)
        graph = trimer_graph()
        bhh = BoseHubbardHamiltonian(graph, basis)

        for J in J_arr
            for eps in eps_arr
                H_corr = f_H_correct(J, eps)
                H_calc = matrix(bhh,J,eps,[0,0,0])

                @test H_corr == H_calc
            end
        end
    end


    @testset "hamiltonian trimer N=2" begin

        N = 2
        k = 3

        function f_H_correct(J,eps,U)
            H = Matrix{Float64}(I, 6, 6)
            H[1,1] = U[3] + eps[3]
            H[1,2] = -J[2]/sqrt(2)
            H[1,3] = 0
            H[1,4] = 0
            H[1,5] = 0
            H[1,6] = 0
            H[2,1] = -J[2]/sqrt(2)
            H[2,2] = (eps[2] + eps[3])/2
            H[2,3] = -J[2]/sqrt(2)
            H[2,4] = -J[1]*0.5
            H[2,5] = 0
            H[2,6] = 0
            H[3,1] = 0
            H[3,2] = -J[2]/sqrt(2)
            H[3,3] = U[2] + eps[2]
            H[3,4] = 0
            H[3,5] = -J[1]/sqrt(2)
            H[3,6] = 0
            H[4,1] = 0
            H[4,2] = -J[1]*0.5
            H[4,3] = 0
            H[4,4] = (eps[1] + eps[3])/2
            H[4,5] = -J[2]*0.5
            H[4,6] = 0
            H[5,1] = 0
            H[5,2] = 0
            H[5,3] = -J[1]/sqrt(2)
            H[5,4] = -J[2]*0.5
            H[5,5] = (eps[1] + eps[2])/2
            H[5,6] = -J[1]/sqrt(2)
            H[6,1] = 0
            H[6,2] = 0
            H[6,3] = 0
            H[6,4] = 0
            H[6,5] = -J[1]/sqrt(2)
            H[6,6] = U[1] + eps[1]

            return H
        end

        # Define test values
        U_arr = 0:2:5
        #U_arr = 1
        U_arr = Iterators.product(U_arr,U_arr,U_arr)
        J_arr = 0:2:5
        #J_arr = 2
        J_arr = Iterators.product(J_arr,J_arr)
        eps_arr = -0.25:0.25:0.25
        #eps_arr = 3
        eps_arr = Iterators.product(eps_arr, eps_arr, eps_arr)

        # Prepare the system
        basis = LtrAscBasis(k, N)
        graph = trimer_graph()
        bhh = BoseHubbardHamiltonian(graph, basis)

        for U in U_arr
            for J in J_arr
                for eps in eps_arr

                    H_corr = f_H_correct(J, eps, U)
                    H_calc = matrix(bhh, J, eps, U)

                    @test isapprox(H_corr, H_calc; atol=1e-14)
                end
            end
        end
    end


end
