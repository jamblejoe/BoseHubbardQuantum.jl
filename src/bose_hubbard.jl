##############################################################################
#
# Basis
#
##############################################################################
import Base: length, iterate, isequal, getindex

"""
For given number of particles N and site count k returns the size of the Hilbert space
The Hilbert space size is binom(N+k-1, k-1).
"""
function bose_hubbard_hilbert_space_size(k::Integer, N::Integer)
  return binomial(N+k-1,k-1)
end

function bose_hubbard_hilbert_space_size(parameters::Dict)
    k = parameters["k"]
    N = parameters["N"]
  return bose_hubbard_hilbert_space_size(k, N)
end



# Abstract Basis is not implemented yet and is just a placeholder
# one should think of AbstractBasis needing to be abstract type
# or needing to have fields
abstract type AbstractBasis end

##############################################################################


struct LtrAscBasis <: AbstractBasis
	k::Int
	N::Int
	basis::Vector{Vector{Int}}
	index::Dict{Vector{Int}, Int}
end


"""
Creates a list of basis states for
N particles
k sites
The states are ordered in ascending order from left to right, e.g.
for 2 particles and 3 sites
(0,0,2) -> (0,1,1) -> (0,2,0) -> (1,0,1) -> (1,1,0) -> (2,0,0).
The length of the list is binomial(N+k-1,k-1).
"""
function LtrAscBasis(parameters::Dict)
	k = parameters["k"]
	N = parameters["N"]
  LtrAscBasis(k,N)
end

function LtrAscBasis(k::Integer, N::Integer)

	@assert k>0
    @assert N>0
    D = bose_hubbard_hilbert_space_size(k, N)

	# create the basis
	state = zeros(typeof(N), k)
    basis = typeof(state)[]

    ltr_asc_loop!(basis, state, N, 1)
	# calculate the index
	index = create_basis_index(basis)

    return LtrAscBasis(k, N, basis, index)
end

function ltr_asc_loop!(basis, state, n::Integer, pos::Integer)
    k = length(state)
    if pos < k
        for i in 0:n
            state[pos] = i
            ltr_asc_loop!(basis, state, n-i, pos+1)
        end
    else
        state[pos] = n
        push!(basis, copy(state))
    end
end

length(basis::LtrAscBasis) = length(basis.basis)
sites(basis::LtrAscBasis) = basis.k
iterate(basis::LtrAscBasis) = iterate(basis.basis)
iterate(basis::LtrAscBasis, state) = iterate(basis.basis, state)
isequal(b1::LtrAscBasis, b2::LtrAscBasis) = b1.k == b2.k && b1.N == b2.N

getposition(basis::LtrAscBasis, state::AbstractVector) = basis.index[state]


##############################################################################

"""
https://arxiv.org/pdf/1410.7280.pdf
"""
struct PonomarevBasis <: AbstractBasis
    k::Int
    N::Int
end

PonomarevBasis(parameters::Dict) = PonomarevBasis(parameters["k"], parameters["N"])

length(basis::PonomarevBasis) = bose_hubbard_hilbert_space_size(basis.k,basis.N)
sites(basis::PonomarevBasis) = basis.k
isequal(b1::PonomarevBasis, b2::PonomarevBasis) = b1.k == b2.k && b1.N == b2.N

#=
Think of caching, precomputing the binomials. Pascals triangle?
=#
"""
Returns the position of a state in the basis.
"""
function getposition(basis::PonomarevBasis, state::AbstractVector)
    k = basis.k
    N = basis.N

    @assert length(state) == k
    @assert all(0 .<= state .<= N)
    @assert sum(state) == N
    
    index = 1
    l = N
    for (m_l,n_i) in enumerate(state)
        for _ in 1:n_i
            index += bose_hubbard_hilbert_space_size(k-m_l, l)
            l -= 1
        end
    end
    index
end

"""
Returns state of the basis with position index.
"""
function getstate(basis::PonomarevBasis, index::Integer)
    state = zeros(Int, basis.k)
    getstate!(state, basis, index)
end

function getstate!(state::AbstractVector{<:Integer}, basis::PonomarevBasis, index::Integer)
    k = basis.k
    N = basis.N

    @assert 1 <= index <= length(basis) "index $index out of bounds [1,$(length(basis))]"

    state .= 0

    for i in N:-1:1
        m_i = 0
        while bose_hubbard_hilbert_space_size(m_i+1,i) < index
            m_i += 1
        end
        state[k-m_i] += 1
        index -= bose_hubbard_hilbert_space_size(m_i,i)
    end
    state
end


"""
Returns the basis element at position index.
"""
getindex(basis::PonomarevBasis, index::Integer) = getstate(basis::PonomarevBasis, index::Integer)


function iterate(basis::PonomarevBasis, state=1)
    count = state
    count > length(basis) ? nothing : ( getstate(basis, count), count+1 )
end


##############################################################################

struct LtrAscBasis0N <: AbstractBasis
	k::Int
	N::Int
end

LtrAscBasis0N(parameters::Dict) = LtrAscBasis0N(parameters["k"], parameters["N"])

length(basis::LtrAscBasis0N) = (basis.N+1)^basis.k
sites(basis::LtrAscBasis0N) = basis.k
#=
This probably can be made non-allocating
=#
function iterate(basis::LtrAscBasis0N, state=(zeros(Int, basis.k), 1))
	k = basis.k
	N = basis.N
	basis_state, count = state


	if count > (N+1)^k
		return nothing
	end

	(basis_state, ( reverse(digits(count, base = N+1, pad = k)) , count+1 ) )
end
isequal(b1::LtrAscBasis0N, b2::LtrAscBasis0N) = b1.k == b2.k && b1.N == b2.N

function getposition(basis::LtrAscBasis0N, state::AbstractVector)
	k = basis.k
	N = basis.N

	@assert length(state) == k
	@assert all(0 .<= state .<= N)
	index = 1
	for i in eachindex(state)
		index += (N+1)^(k-i) * state[i]
	end
	index
end


##############################################################################




STD_BASIS = LtrAscBasis

# the return type of this function is not defined at compile time!!
function get_or_create_basis(parameters::Dict)
	haskey(parameters, "basis") && return parameters["basis"]
	STD_BASIS(parameters)
end


"""
A helper function which creates a dictionary assigning basis elements
to their position/index in the basis for quick lookup of the index.
"""
function create_basis_index(basis)
    basis_element_type = eltype(basis)
	basis_element_element_type = eltype(basis_element_type)
    basis_to_index = Dict{basis_element_type, basis_element_element_type}()
    for (i,basis_element) in enumerate(basis)
        push!(basis_to_index, basis_element => i)
    end
    return basis_to_index
end




##############################################################################
#
# Predefined graphs, e.g. representing lattices
#
##############################################################################

"""
Creates a dimer graph
resulting in a Hamiltonian
H = -J/2 (a_dag_1 a_2 + a_dag_2 a_1) + U/2 (n_1(n_1-1) + n_2(n_2-1))
  + delta/2 (n_1-n_2)
"""
function dimer_graph()
    return path_graph(2)
end

"""
Creates a trimer graph
"""
function trimer_graph()
    return path_graph(3)
end

"""
Creates a chain graph of length k
"""
function k_chain_graph(k::Integer)
    return path_graph(k)
end






##############################################################################
#
# Various operators
#
##############################################################################

"""
Creates a tunnel operator from source to destination in the standard basis.
"""
function tunnel_operator(source::Integer, destination::Integer, parameters::Dict)
    basis = get_or_create_basis(parameters)
    tunnel = tunnel_spmatrix(source,destination,basis)
    return tunnel
end

"""
Creates a sparse matrix representing a tunnel operator from source to destination
in the given basis.
"""
function tunnel_spmatrix(source::Integer, destination::Integer,
	basis::AbstractBasis)

	# ask the discord forum about the following line
	#@assert typeof(source) == typeof(destination)

	k = basis.k
	N = basis.N

    basis_length = length(basis)
    basis_element_length = k

    @assert 1 <= source <= k && 1 <= destination <= k

    # create the tunneling vector
    tunnel_operator = zeros(Int, basis_element_length)
    tunnel_operator[source] = -1
    tunnel_operator[destination] = 1


    rows = Int[]
    cols = Int[]
    values = Float64[]	# is Float64 really necessary here?

    # create a helper dictionary to quickly look up indices of basis elements
    #basis_to_index = basis.index


    for (i,basis_element) in enumerate(basis)
        # Create tunnel states
        # let the tunneling operator act on the basis_element
        neighbouring_state = tunnel_operator .+ basis_element

        # check if the "tunneled"/neighbouring state is a valid state
        if all(x -> 0 <= x <= N, neighbouring_state)

            # get the index of the neighbouring state
            neighbour_index = getposition(basis, neighbouring_state)

            # store the indeces
            push!(rows, neighbour_index)
            push!(cols, i)

            # calculate and store the value
            value = sqrt((basis_element[destination]+1)*basis_element[source])
            push!(values, value)

        end
    end

    # create the sparse matrices from the rows/columns/values
    tunnel = sparse(rows, cols, values, basis_length, basis_length)

    return tunnel
end

"""
Creates a symmetric tunnel operator between site_1 and site_2 in the standard
basis.
"""
function tunnel_symm_operator(site_1::Integer, site_2::Integer, parameters::Dict)
    k = parameters["k"]
    N = parameters["N"]

    basis = get_or_create_basis(parameters)

    tunnel = tunnel_symm_spmatrix(site_1,site_2,basis)
    return tunnel
end


"""
Create a sparse matrix representing a symmetric tunnel operator between
site_1 and site_2 in the given basis.
"""
function tunnel_symm_spmatrix(site_1::Integer, site_2::Integer, basis::AbstractBasis)
    t = tunnel_spmatrix(site_1, site_2, basis)
    t += transpose(t)
    return t
end

"""
Creates the kinetic operator in the standard basis.
"""
function kinetic_operator(parameters::Dict)
    basis = get_or_create_basis(parameters)
    op = kinetic_spmatrix(basis)
    return op
end

"""
Creates a sparse matrix representing the kinetic operator in the given basis.
"""
function kinetic_spmatrix(basis::AbstractBasis)
	k = basis.k
	N = basis.N

	@assert k > 1

	op = tunnel_symm_spmatrix(1,2,basis)
	for i in 2:(k-1)
		op += tunnel_symm_spmatrix(i,i+1,basis)
	end
	op
end



"""
Create a number operator on the given site in the standard basis.
"""
function number_operator(site::Integer, parameters::Dict)
    k = parameters["k"]
    N = parameters["N"]

    basis = LtrAscBasis(k, N)

    tunnel = number_spmatrix(site, basis)
    return tunnel
end


"""
Create a sparse matrix representing the number operator on the
given site in the given basis.
"""
function number_spmatrix(site::Integer, basis::AbstractBasis)
		k = basis.k
		basis_length = length(basis)

    @assert 1 <= site <= k

    # create arrays to store the indices and values of the sparse matrix
    # matrix will only have diagonal non-zero entries
    diags = Int[]
    values = Float64[]  # again, really necessary?

    for (i, basis_element) in enumerate(basis)
        value = basis_element[site]
        if value != 0
            push!(diags, i)
            push!(values, value)
        end
    end

    number = sparse(diags, diags, values, basis_length, basis_length)

    return number
end


##############################################################################
#
# Optimized Bose Hubbard Hamiltonian creation
#
##############################################################################

"""
Representation of a Bose-Hubbard Hamiltonian with
N particles and
k sites.
"""
struct BoseHubbardHamiltonian
    N 			# the number of particles
    k 			# the number of sites
    tunnel 		# stores the tunneling matrices defined by the edges of given graphs
    #interaction # stores the interaction matrices defined by the nodes of given graph
	potential   # stores the potential matrices defined by the nodes of the given graph
end

function BoseHubbardHamiltonian(graph, basis::AbstractBasis)

    k = basis.k
	N = basis.N

	@assert k>0
  	@assert N>0
	@assert k == nv(graph)


    tunnel = tunnel_spmatrices(graph, basis)
    #interaction = interaction_spmatrices(graph, basis)
	potential = potential_spmatrices(graph, basis)

    #return BoseHubbardHamiltonian(N, k, tunnel, interaction)
	return BoseHubbardHamiltonian(N, k, tunnel, potential)

end

function BoseHubbardHamiltonianChain(parameters::Dict)

	k = parameters["k"]
	N = parameters["N"]

	@assert k>0
	@assert N>0

	graph = k_chain_graph(k::Integer)
	basis = get_or_create_basis(parameters)

	return BoseHubbardHamiltonian(graph, basis)
end


"""
Create a Vector of sparse matrices holding the nearest neighbourhood tunneling
matrices assigned to every edge of the given graph/lattice with respect to
the given basis.
"""
function tunnel_spmatrices(graph, basis::AbstractBasis)

	k = basis.k

	@assert k>0
	@assert k == nv(graph)

    #basis_length = length(basis)
    sites_count = k
    tunnel_operators_count = ne(graph)

    tunnel = Vector{SparseMatrixCSC{Float64,Int}}(undef, tunnel_operators_count)

    for (i,edge) in enumerate(edges(graph))
        src_index = src(edge)
        dst_index = dst(edge)

        tunnel_spmatrix = tunnel_symm_spmatrix(src_index, dst_index, basis)
        tunnel[i] = tunnel_spmatrix
    end

    return tunnel
end


"""
Create a Vector of sparse matrices holding the onsite potential
matrices assigned to every vertex of the given graph/lattice with respect to
the given basis.
"""
function potential_spmatrices(graph, basis::AbstractBasis)

	k = basis.k

	@assert k>0
	@assert k == nv(graph)

	potential = Vector{SparseMatrixCSC{Float64,Int}}(undef, k)

	for i in 1:k
		potential[i] = number_spmatrix(i, basis)
	end

	potential
end

#"""
#Create tunneling operators represented by vectors which act on
#basis elements by addition, e.g. the tunneling vector assigned
#to tunneling from site 1 to site 2 in a 3-site model would be
#(-1, 1, 0).
#This tunneling operator then acts e.g. on the state (2,4,1) by
#(-1, 1, 0) + (2, 4, 1) = (1, 5, 1),
#pushing one particle from site 1 to site 2.
#
#Only one-at-a-time particle tunneling operators are created.
#
#Returned are
#a vector of source sites
#a vector of destination sites
#a vector of the tunneling operators/vectors
#"""
#function _create_tunnel_operators(graph)
#
#    # Number of sites
#    sites_count = nv(graph)
#    tunnel_count = ne(graph)
#
#    src_sites = Int64[]
#    dst_sites = Int64[]
#    tunnel_operators = Vector{Int64}[]
#
#    for edge in edges(graph)
#
#        #assert isinstance(edge[0], int) and isinstance(edge[1], int), "edges must connect integer values"
#        src_index = src(edge)
#        dst_index = dst(edge)
#        push!(src_sites, src_index)
#        push!(dst_sites, dst_index)
#
#        tunnel_operator = zeros(sites_count)
#        tunnel_operator[src_index] = -1
#        tunnel_operator[dst_index] = 1
#        push!(tunnel_operators, tunnel_operator)
#    end

#    return (src_sites, dst_sites, tunnel_operators)
#
#end

"""
Create a vector of sparse matrices representing the on-site interaction
of particles.
"""
function interaction_spmatrices(graph, basis::AbstractBasis)

    basis_length = length(basis)
    sites_count = nv(graph)

    interaction = Vector{SparseMatrixCSC{Float64,Integer}}(undef, sites_count)


    for i in 1:sites_count

        # create arrays to store the indices and values of the sparse matrices
        # matrix will only have diagonal non-zero entries
        diags = Int[]
        values = Float64[]

        for (j, basis_state) in enumerate(basis)
            value = basis_state[i]*(basis_state[i]-1)
            if value != 0
                push!(diags, j)
                push!(values, value)
            end
        end

        interaction[i] = sparse(diags, diags, values, basis_length, basis_length)
    end

    return interaction
end


"""
Returns the dense matrix representation of the given Bose-Hubbard
Hamiltonian
"""
function matrix(bhh::BoseHubbardHamiltonian, J, U)
    return Array(spmatrix(bhh, J, U))
end

function matrix(bhh::BoseHubbardHamiltonian, J, eps, U)
    return Array(spmatrix(bhh, J, eps, U))
end

function matrix(bhh::BoseHubbardHamiltonian, parameters::Dict)
	k = parameters["k"]
	J = parameters["J"]
	U = haskey(parameters, "U") ? parameters["U"] : zeros(k)
	U = length(U) == 1 ? [U for _ in 1:k] : U
	eps = haskey(parameters, "eps") ? parameters["eps"] : zeros(k)

	Array(spmatrix(bhh, J, eps, U))
end

"""
Returns the sparse matrix representation of the given Bose-Hubbard
Hamiltonian
"""
function spmatrix(bhh::BoseHubbardHamiltonian, J, eps, U)

    N = bhh.N
    k = bhh.k
    tunnel = bhh.tunnel
    potential = bhh.potential

    D = bose_hubbard_hilbert_space_size(k, N)

    @assert length(J) == length(tunnel) "Number of tunnel parameters must be $(length(tunnel))"
	@assert length(eps) == length(potential) "Number of on-site potential parameters must be $(length(potential))"
    @assert length(U) == length(potential) "Number of on-site interaction parameters must be $(length(potential))"

    H = spzeros(D, D)
	# tunnel part
    for i in 1:length(tunnel)
        H += -J[i] * tunnel[i]
    end

	# potential part
	for i in 1:length(potential)
        H += eps[i] * potential[i]
    end

	# interaction part
	#for i in 1:interaction_count
    #    H += U[i] * interaction[i]
    #end
	for i in 1:length(potential)
        H += U[i] * potential[i] * (potential[i] - I(D))
    end

    H *= 0.5

    return H
end

# TODO
# remove this function at some point
function spmatrix(bhh::BoseHubbardHamiltonian, J, U)
	k = bhh.k
	eps = zeros(k)
	return spmatrix(bhh, J, eps, U)
end





##############################################################################
#
# partial traces
#
##############################################################################

"""
ONLY WORKS FOR 1D CHAINS
NOT OPTIMIZED

traces out subsystem B and returns O acting on A

Operator O
"""
function partialtr(k::Integer, N::Integer, k_B::Integer, O::AbstractMatrix)

	D = bose_hubbard_hilbert_space_size(k, N)

	@assert size(O)[1] == size(O)[2] == D
	@assert 1 <= k_B <= k

	k_A = k - k_B

	basis_full = LtrAscBasis(k, N)
	basis_full_index = create_basis_index(basis_full)
	basis_A = LtrAscBasis0N(k_A, N)
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


			for state_B_col in LtrAscBasis(k_b, N-N_i)
				state_col_cache[k_A+1:end] .= state_B_col

				O_ind_col = basis_full_index[state_col_cache]

				for state_B_row in LtrAscBasis(k_b, N-N_j)
					state_row_cache[k_A+1:end] .= state_B_row

					O_ind_row = basis_full_index[state_row_cache]

					O_A[j,i] += O[O_ind_row, O_ind_col]
				end
			end

		end
	end

	return O_A
end
