##############################################################################
#
# partial traces
#
##############################################################################

#
# The caching structure should be thought more about
#

struct PartialtrCache
	basis_full
	bases_A
end

function PartialtrCache(k, N, k_A)
	basis_full=LtrAscBasis(k, N)
	bases_A=[LtrAscBasis(k_A, n) for n in 0:N]
	PartialtrCache(basis_full, bases_A)
end

#=
"""
ONLY WORKS FOR 1D CHAINS
NOT OPTIMIZED

traces out subsystem B and returns O acting on A

Operator O
"""
=#
#=
function partialtr(k::Integer, N::Integer, k_B::Integer, O::AbstractMatrix)

	D = bose_hubbard_hilbert_space_size(k, N)

	@assert size(O)[1] == size(O)[2] == D
	@assert 1 <= k_B <= k

	k_A = k - k_B

	basis_full = LtrAscBasis(k, N)
	basis_full_index = create_basis_index(basis_full)
	basis_A = LtrAscBasis(k_A, 0, N)
	D_A = length(basis_A)

	O_A = zeros(D_A,D_A)
	state_row_cache = zeros(k)
	state_col_cache = zeros(k)


	for (i, state_A_col) in enumerate(basis_A)
        N_i = sum(state_A_col)
		state_col_cache[1:k_A] .= state_A_col

		for (j, state_A_row) in enumerate(basis_A)
			N_j = sum(state_A_row)
			state_row_cache[1:k_A] .= state_A_row

            if N_i == N_j
                for state_B in LtrAscBasis(k_B, N-N_i)
                    state_col_cache[k_A+1:end] .= state_B
                    state_row_cache[k_A+1:end] .= state_B

                    O_ind_col = basis_full_index[state_col_cache]
                    O_ind_row = basis_full_index[state_row_cache]
                    O_A[j,i] += O[O_ind_row, O_ind_col]

                    #for state_B_row in LtrAscBasis(k_B, N-N_j)
                    #	state_row_cache[k_A+1:end] .= state_B_row
                    #
                    #	O_ind_row = basis_full_index[state_row_cache]
                    #
                    #	O_A[j,i] += O[O_ind_row, O_ind_col]
                    #end
                end
            end

		end
	end

	return O_A
end
=#


"""
Partial trace for operators. There is a faster version for partial traces of
density matrices of a single pure state.
"""
function partialtr!(O_A::AbstractMatrix, O::AbstractMatrix,
		basis_full, basis_A, basis_B;
		#cache::PartialtrCache=PartialtrCache(k, N, k_A)
    )



	D = length(basis_full)
	D_A = length(basis_A)
	D_B = length(basis_B)
	D != D_A * D_B && error("D != D_A * D_B")

	# make sure that O_A is initialized to 0
    O_A .= 0

	k = basis_full.k
	k_A = basis_A.k
	k_B = basis_B.k

	state_row_cache = zeros(Int, k)
	state_col_cache = zeros(Int, k)

	for (i, state_A_col) in enumerate(basis_A)
		N_i = sum(state_A_col)
		state_col_cache[1:k_A] .= state_A_col

		for (j, state_A_row) in enumerate(basis_A)
			N_j = sum(state_A_row)
			state_row_cache[1:k_A] .= state_A_row

			for state_B in basis_B
				state_col_cache[k_A+1:end] .= state_B
				state_row_cache[k_A+1:end] .= state_B

				O_ind_col = getposition(basis_full, state_col_cache)
				O_ind_row = getposition(basis_full, state_row_cache)
				O_A[j,i] += O[O_ind_row, O_ind_col]
			end


		end
	end


	return O_A
end

#=
function partialtr(k::Integer, N::Integer, k_A::Integer, O::AbstractMatrix)
    D_A = bose_hubbard_hilbert_space_size(k_A, 0, N)
    O_A = zeros(D_A,D_A)
    return partialtr!(O_A, k, N, k_A, O)
end
=#
#=
function partialtr!(O_A::AbstractMatrix, k::Integer, N::Integer, k_A::Integer, O::AbstractMatrix;
		cache::PartialtrCache=PartialtrCache(k, N, k_A)
    )

	basis_full = cache.basis_full
	bases_A = cache.bases_A

    O_A .= 0

	D = bose_hubbard_hilbert_space_size(k, N)

	@assert size(O)[1] == size(O)[2] == D
	@assert 1 <= k_A <= k

	k_B = k - k_A

	#basis_full = LtrAscBasis(k, N)
    #basis_full_index = create_basis_index(basis_full)

	#basis_A = [LtrAscBasis(k_A, n) for n in 0:N]
	#D_A = bose_hubbard_hilbert_space_size(k_A, 0, N)

	#O_A = zeros(D_A,D_A)
	state_row_cache = zeros(k)
	state_col_cache = zeros(k)

    offset = 0
    for basis_A in bases_A

        for (i, state_A_col) in enumerate(basis_A)
            N_i = sum(state_A_col)
            state_col_cache[1:k_A] .= state_A_col

            for (j, state_A_row) in enumerate(basis_A)
                N_j = sum(state_A_row)
                state_row_cache[1:k_A] .= state_A_row

                for state_B in LtrAscBasis(k_B, N-N_i)
                    state_col_cache[k_A+1:end] .= state_B
                    state_row_cache[k_A+1:end] .= state_B

                    #O_ind_col = basis_full_index[state_col_cache]
                    #O_ind_row = basis_full_index[state_row_cache]
					O_ind_col = getposition(basis_full, state_col_cache)
                    O_ind_row = getposition(basis_full, state_row_cache)
                    O_A[offset+j,offset+i] += O[O_ind_row, O_ind_col]

                    #for state_B_row in LtrAscBasis(k_B, N-N_j)
                    #	state_row_cache[k_A+1:end] .= state_B_row
                    #
                    #	O_ind_row = basis_full_index[state_row_cache]
                    #
                    #	O_A[j,i] += O[O_ind_row, O_ind_col]
                    #end
                end


            end
        end

        offset += length(basis_A)
    end

	return O_A
end
=#


function partialtr(ψ::AbstractVector,
	basis_full::AbstractBasis, basis_A::AbstractBasis, basis_B::AbstractBasis;
		sqrtρ = zeros(length(basis_A), length(basis_B)))
	D_A = length(basis_A)
	O_A = zeros(D_A, D_A)
	partialtr!(O_A, ψ, basis_full, basis_A, basis_B; sqrtρ=sqrtρ)
end



function partialtr!(O_A::AbstractMatrix, ψ::AbstractVector,
        basis_full::AbstractBasis, basis_A::AbstractBasis, basis_B::AbstractBasis;
		sqrtρ = zeros(length(basis_A), length(basis_B))
    )

	D = length(basis_full)
	D_A = length(basis_A)
	D_B = length(basis_B)
	D != D_A * D_B && error("D != D_A * D_B")

	# this is for safety purposes to ensure that sqrtρ really is initialized
	# with zeros
	sqrtρ .= 0

	k_A = basis_A.k

	for (i, state) in enumerate(basis_full)
		@views j = getposition(basis_A, state[1:k_A])
		@views l = getposition(basis_B, state[k_A+1:end])
		sqrtρ[j,l] = ψ[i]
	end

	mul!(O_A, sqrtρ, transpose(sqrtρ))
end
