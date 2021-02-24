"""
Created on 10 Sep 2019


TODO: 
    - add more tests for LtrAscBasis0N combined with operators
"""

using BoseHubbardQuantum
using Test
using SparseArrays
using LinearAlgebra

import BoseHubbardQuantum: tunnel_spmatrices

@testset "LtrAscBasis" begin

    @testset "length" begin
        for k in 3:5
        for Nmin in 2:10
        for Nmax in Nmin:10
            basis = LtrAscBasis(k, Nmin, Nmax)
            l = 0
            for _ in basis
                l += 1
            end
            @test l == length(basis)
        end
        end
        end

        Nmin = 0
        for k in 3:5
        for Nmax in Nmin:10
            basis = LtrAscBasis(k, Nmin, Nmax)
            l = 0
            for _ in basis
                l += 1
            end
            @test l == bose_hubbard_hilbert_space_size(k, Nmin, Nmax)
        end
        end
    end


    @testset "N=1, k=1" begin
        # Define the correct basis list
        basis_corr = [[0,1],[1,0]]
        basis_calc = LtrAscBasis(2,1).basis

        @test basis_corr == basis_calc
    end

    @testset "N=2, k=3" begin
        # Define the correct basis list
        basis_corr = [[0,0,2],
                      [0,1,1],
                      [0,2,0],
                      [1,0,1],
                      [1,1,0],
                      [2,0,0]]
        basis_calc = LtrAscBasis(3,2).basis
        @test basis_corr == basis_calc
    end

    @testset "k=2, Nmin=0, Nmax=1" begin
        # Define the correct basis list
        basis_corr = [[0,0],[0,1],[1,0]]
        basis_calc = LtrAscBasis(2,0,1).basis

        @test basis_corr == basis_calc
    end

    @testset "k=3, Nmin=0, Nmax=2" begin
        # Define the correct basis list
        basis_corr = [[0,0,0],
                      [0,0,1],
                      [0,1,0],
                      [1,0,0],
                      [0,0,2],
                      [0,1,1],
                      [0,2,0],
                      [1,0,1],
                      [1,1,0],
                      [2,0,0]]
        basis_calc = LtrAscBasis(3,0,2).basis
        @test basis_corr == basis_calc
    end

    @testset "k=3, Nmin=1, Nmax=2" begin
        # Define the correct basis list
        basis_corr = [[0,0,1],
                      [0,1,0],
                      [1,0,0],
                      [0,0,2],
                      [0,1,1],
                      [0,2,0],
                      [1,0,1],
                      [1,1,0],
                      [2,0,0]]
        basis_calc = LtrAscBasis(3,1,2).basis
        @test basis_corr == basis_calc
    end


end


@testset "PonomarevBasis" begin

    @testset "length" begin
    for k in 3:5
        for N in 2:10
            basis = PonomarevBasis(k, N)
            l = 0
            for _ in basis
                l += 1
            end
            @test l == length(basis)
        end
        end
    end

    @testset "getposition" begin
        @testset "k=3, N=2" begin
            k = 3
            N = 2
            basis = PonomarevBasis(k,N)
            @test getposition(basis, [0,0,2]) == 1
            @test getposition(basis, [0,1,1]) == 2
            @test getposition(basis, [0,2,0]) == 3
            @test getposition(basis, [1,0,1]) == 4
            @test getposition(basis, [1,1,0]) == 5
            @test getposition(basis, [2,0,0]) == 6

        end


    end

    @testset "getstate" begin
        @testset "k=3, N=2" begin
            k = 3
            N = 2
            basis = PonomarevBasis(k,N)
            @test all(getstate(basis, 1) .== [0,0,2])
            @test all(getstate(basis, 2) .== [0,1,1])
            @test all(getstate(basis, 3) .== [0,2,0])
            @test all(getstate(basis, 4) .== [1,0,1])
            @test all(getstate(basis, 5) .== [1,1,0])
            @test all(getstate(basis, 6) .== [2,0,0])
        end
    end


    @testset "getposition ∘ getstate == id" begin
        for k in 3:6
        for N in 2:10
            basis = PonomarevBasis(k,N)
            for i in 1:bose_hubbard_hilbert_space_size(k,N)
                @test getposition(basis, getstate(basis, i)) == i
            end
        end
        end
    end


    @testset "index map is surjective" begin
    for k in 3:6
        for N in 2:10
            basis = PonomarevBasis(k,N)
            indices = [getposition(basis, state) for state in basis]
            @test all(indices .== sort(indices))
            @test all(indices .== 1:bose_hubbard_hilbert_space_size(k,N))
        end
        end
    end

    @testset "getstate! with views" begin
        k = 6
        N = 10
        D = bose_hubbard_hilbert_space_size(k,N)
        basis = PonomarevBasis(k,N)
        basis_states = zeros(k,D)
        for i in 1:D
            @views getstate!(basis_states[:,i], basis, i)
        end

    end



end

#=
@testset "LtrAsc0NBasis" begin

    @testset "length" begin
        for k in 3:5
        for N in 2:10
            basis = LtrAscBasis0N(k, N)
            l = 0
            for _ in basis
                l += 1
            end
            @test l == length(basis)
        end
        end
    end

    @testset "getposition" begin

        @testset "k=3, N=1" begin
            k = 3
            N = 1
            basis = LtrAscBasis0N(k, N)
            @test getposition(basis, [0,0,0]) == 1
            @test getposition(basis, [0,0,1]) == 2
            @test getposition(basis, [0,1,0]) == 3
            @test getposition(basis, [0,1,1]) == 4
            @test getposition(basis, [1,0,0]) == 5
            @test getposition(basis, [1,0,1]) == 6
            @test getposition(basis, [1,1,0]) == 7
            @test getposition(basis, [1,1,1]) == 8
        end

        @testset "k=3, N=2" begin
            k = 3
            N = 2
            basis = LtrAscBasis0N(k, N)
            @test getposition(basis, [0,0,0]) == 1
            @test getposition(basis, [0,0,1]) == 2
            @test getposition(basis, [0,0,2]) == 3
            @test getposition(basis, [0,1,0]) == 4
            @test getposition(basis, [0,1,1]) == 5
            @test getposition(basis, [0,1,2]) == 6
            @test getposition(basis, [0,2,0]) == 7
            @test getposition(basis, [0,2,1]) == 8
            @test getposition(basis, [0,2,2]) == 9
            @test getposition(basis, [1,0,0]) == 10
            @test getposition(basis, [1,0,1]) == 11
            @test getposition(basis, [1,0,2]) == 12
            @test getposition(basis, [1,1,0]) == 13
            @test getposition(basis, [1,1,1]) == 14
            @test getposition(basis, [1,1,2]) == 15
            @test getposition(basis, [1,2,0]) == 16
            @test getposition(basis, [1,2,1]) == 17
            @test getposition(basis, [1,2,2]) == 18
            @test getposition(basis, [2,0,0]) == 19
            @test getposition(basis, [2,0,1]) == 20
            @test getposition(basis, [2,0,2]) == 21
            @test getposition(basis, [2,1,0]) == 22
            @test getposition(basis, [2,1,1]) == 23
            @test getposition(basis, [2,1,2]) == 24
            @test getposition(basis, [2,2,0]) == 25
            @test getposition(basis, [2,2,1]) == 26
            @test getposition(basis, [2,2,2]) == 27
        end

    end




end
=#

@testset "Test operators" begin
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

end


@testset "Test Hamiltonian creation" begin



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
