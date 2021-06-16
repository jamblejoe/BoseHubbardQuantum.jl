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

    @testset "getposition" begin
        @testset "k=3, N=2" begin
            k = 3
            N = 2
            basis = LtrAscBasis(k,N)
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
            basis = LtrAscBasis(k,N)
            @test all(getstate(basis, 1) .== [0,0,2])
            @test all(getstate(basis, 2) .== [0,1,1])
            @test all(getstate(basis, 3) .== [0,2,0])
            @test all(getstate(basis, 4) .== [1,0,1])
            @test all(getstate(basis, 5) .== [1,1,0])
            @test all(getstate(basis, 6) .== [2,0,0])
        end
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

@testset "LtrAscCutoffBasis" begin
    @testset "getposition" begin

        @testset "k=3, d=1" begin
            k = 3
            d = 1
            basis = LtrAscCutoffBasis(k, d)
            @test getposition(basis, [0,0,0]) == 1
            @test getposition(basis, [0,0,1]) == 2
            @test getposition(basis, [0,1,0]) == 3
            @test getposition(basis, [0,1,1]) == 4
            @test getposition(basis, [1,0,0]) == 5
            @test getposition(basis, [1,0,1]) == 6
            @test getposition(basis, [1,1,0]) == 7
            @test getposition(basis, [1,1,1]) == 8
        end

        @testset "k=3, d=2" begin
            k = 3
            d = 2
            basis = LtrAscCutoffBasis(k, d)
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


    @testset "getstate" begin

        @testset "k=3,d=1" begin
            k = 3
            d = 1
            basis = LtrAscCutoffBasis(k, d)
            @test getstate(basis, 1) == [0,0,0]
            @test getstate(basis, 2) == [0,0,1]
            @test getstate(basis, 3) == [0,1,0]
            @test getstate(basis, 4) == [0,1,1]
            @test getstate(basis, 5) == [1,0,0]
            @test getstate(basis, 6) == [1,0,1]
            @test getstate(basis, 7) == [1,1,0]
            @test getstate(basis, 8) == [1,1,1]

        end

        @testset "k=3,d=2" begin
            k = 3
            d = 2
            basis = LtrAscCutoffBasis(k, d)
            @test getstate(basis, 1) == [0,0,0]
            @test getstate(basis, 2) == [0,0,1]
            @test getstate(basis, 3) == [0,0,2]
            @test getstate(basis, 4) == [0,1,0]
            @test getstate(basis, 5) == [0,1,1]
            @test getstate(basis, 6) == [0,1,2]
            @test getstate(basis, 7) == [0,2,0]
            @test getstate(basis, 8) == [0,2,1]
            @test getstate(basis, 9) == [0,2,2]
            @test getstate(basis, 10) == [1,0,0]
            @test getstate(basis, 11) == [1,0,1]
            @test getstate(basis, 12) == [1,0,2]
            @test getstate(basis, 13) == [1,1,0]
            @test getstate(basis, 14) == [1,1,1]
            @test getstate(basis, 15) == [1,1,2]
            @test getstate(basis, 16) == [1,2,0]
            @test getstate(basis, 17) == [1,2,1]
            @test getstate(basis, 18) == [1,2,2]
            @test getstate(basis, 19) == [2,0,0]
            @test getstate(basis, 20) == [2,0,1]
            @test getstate(basis, 21) == [2,0,2]
            @test getstate(basis, 22) == [2,1,0]
            @test getstate(basis, 23) == [2,1,1]
            @test getstate(basis, 24) == [2,1,2]
            @test getstate(basis, 25) == [2,2,0]
            @test getstate(basis, 26) == [2,2,1]
            @test getstate(basis, 27) == [2,2,2]

        end

    end

    @testset "getstate and getposition" begin

        ks = 2:4
        ds = 1:3
        for k in ks
            for d in ds
                basis = LtrAscCutoffBasis(k,d)
                for j in 1:length(basis)
                    @test getposition(basis, getstate(basis, j)) == j
                end
            end
        end
    end



    @testset "iterate" begin

        @testset "k=3, d=1" begin

            k = 3
            d = 1
            basis = LtrAscCutoffBasis(k,d)
            for (i,state) in enumerate(basis)
                @test state == getstate(basis,i)
            end

        end

    end

    @testset "is in" begin
        k = 3
        d = 2
        basis = LtrAscCutoffBasis(k,d)
        @test !(1 in basis)
        @test !([0,0] in basis)
        @test [0,0,0] in basis
        @test [0,0,1] in basis
        @test !([0,0,-1] in basis)
        @test !([3,0,0] in basis)
    end




end
