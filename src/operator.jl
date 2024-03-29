##############################################################################
#
# Various operators
#
##############################################################################
"""
Creates a sparse matrix representing the creation operator on site index with
respect to a given basis.
"""
creation_spmatrix(index::Integer, basis::AbstractBasis) = creation_spmatrix(index, basis, basis)

"""
Creates a sparse matrix representing the creation operator on site index with
respect to given bases. The matrix will be mapping from basis1 to basis2, i.e.
the domain is basis1 and the codomain is basis2.
"""
function creation_spmatrix(index::Integer, basis1::AbstractBasis, basis2::AbstractBasis)

	L = sites(basis1)
	@assert L == sites(basis2)

    @assert 1 <= index <= L

    # operator acting on basis state by elementwise addition
    op = zeros(Int, L)
    op[index] = 1

    rows = Int[]
    cols = Int[]
    values = Float64[]	# is Float64 really necessary here?


    for (i,state) in enumerate(basis1)

        dst_state = op .+ state

        if dst_state in basis2
            # get the index of the destination state
            dst_index = getposition(basis2, dst_state)

            # store the indices
            push!(rows, dst_index)
            push!(cols, i)

            # calculate and store the value
            value = sqrt((state[index]+1))
            push!(values, value)

        end
    end

    # create the sparse matrix from the rows/columns/values
    return sparse(rows, cols, values, length(basis2), length(basis1))

end


annihilation_spmatrix(index::Integer, basis::AbstractBasis) = annihilation_spmatrix(index, basis, basis)
"""
Creates a sparse matrix representing the creation operator on site index with
respect to the given basis.
"""
function annihilation_spmatrix(index::Integer, basis1::AbstractBasis, basis2::AbstractBasis)

	L = sites(basis1)
	@assert L == sites(basis2)
    @assert 1 <= index <= L

    # create the tunneling vector
    op = zeros(Int, L)
    op[index] = -1

    rows = Int[]
    cols = Int[]
    values = Float64[]	# is Float64 really necessary here?


    for (i,state) in enumerate(basis1)

        dst_state = op .+ state

        if dst_state in basis2
            dst_index = getposition(basis2, dst_state)

            # store the indices
            push!(rows, dst_index)
            push!(cols, i)

            # calculate and store the value
            value = sqrt((state[index]))
            push!(values, value)

        end
    end

    # create the sparse matrix from the rows/columns/values
    return sparse(rows, cols, values, length(basis2), length(basis1))

end




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

	L = sites(basis)
	#N = basis.N

    basis_length = length(basis)
    basis_element_length = L

    @assert 1 <= source <= L && 1 <= destination <= L

    # create the tunneling vector
    tunnel_operator = zeros(Int, basis_element_length)
    tunnel_operator[source] = -1
    tunnel_operator[destination] = 1

    neighbouring_state = similar(tunnel_operator)
    basis_element = similar(tunnel_operator)

    rows = Int[]
    cols = Int[]
    values = Float64[]	# is Float64 really necessary here?

    # create a helper dictionary to quickly look up indices of basis elements
    #basis_to_index = basis.index


    #for (i,basis_element) in enumerate(basis)
    for i in eachindex(basis)
        getstate!(basis_element, basis, i)
        # Create tunnel states
        # let the tunneling operator act on the basis_element
        neighbouring_state .= tunnel_operator .+ basis_element


        if neighbouring_state in basis

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
    L = _get_L(parameters)
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
	L = sites(basis)

	@assert L > 1

	op = tunnel_symm_spmatrix(1,2,basis)
	for i in 2:(L-1)
		op += tunnel_symm_spmatrix(i,i+1,basis)
	end
	op
end



"""
Create a number operator on the given site in the standard basis.
"""
function number_operator(site::Integer, parameters::Dict)
    L = _get_L(parameters)
    N = parameters["N"]

    basis = LtrAscBasis(L, N)

    tunnel = number_spmatrix(site, basis)
    return tunnel
end


"""
Create a sparse matrix representing the number operator on the
given site in the given basis.
"""
function number_spmatrix(site::Integer, basis::AbstractBasis)
    L = sites(basis)

    @assert 1 <= site <= L

    # create arrays to store the indices and values of the sparse matrix
    # matrix will only have diagonal non-zero entries
    diags = Int[]
    values = Float64[]  # again, really necessary?

    basis_element = zeros(L)

    #for (i, basis_element) in enumerate(basis)
    for i in eachindex(basis)
        getstate!(basis_element, basis, i)
        value = basis_element[site]
        if value != 0
            push!(diags, i)
            push!(values, value)
        end
    end

    number = sparse(diags, diags, values, length(basis), length(basis))

    return number
end

function total_number_spmatrix(basis::AbstractBasis)
    L = sites(basis)
    basis_length = length(basis)

    # create arrays to store the indices and values of the sparse matrix
    # matrix will only have diagonal non-zero entries
    diags = Int[]
    values = Float64[]  # again, really necessary?

    for (i, basis_element) in enumerate(basis)
        value = sum(basis_element)
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
    #N 			# the number of particles
    L 			# the number of sites
    #D           # Dimension of the system, i.e. the length of the underlying basis
    basis       # the underlying basis
    tunnel 		# stores the tunneling matrices defined by the edges of given graphs
    #interaction # stores the interaction matrices defined by the nodes of given graph
	potential   # stores the potential matrices defined by the nodes of the given graph
end

function BoseHubbardHamiltonian(graph, basis::AbstractBasis)

    L = sites(basis)
    #N = basis.N
    #D = length(basis)

	@assert L>0
  	#@assert N>0
	@assert L == nv(graph)


    tunnel = tunnel_spmatrices(graph, basis)
    #interaction = interaction_spmatrices(graph, basis)
	potential = potential_spmatrices(graph, basis)

    #return BoseHubbardHamiltonian(N, k, tunnel, interaction)
	return BoseHubbardHamiltonian(L, basis, tunnel, potential)

end

function BoseHubbardHamiltonianChain(parameters::Dict)

	L = _get_L(parameters)
	N = parameters["N"]

	@assert L>0
	@assert N>0

	graph = chain_graph(L::Integer)
	basis = get_or_create_basis(parameters)

	return BoseHubbardHamiltonian(graph, basis)
end


"""
Create a Vector of sparse matrices holding the nearest neighbourhood tunneling
matrices assigned to every edge of the given graph/lattice with respect to
the given basis.
"""
function tunnel_spmatrices(graph, basis::AbstractBasis)

	L = sites(basis)

	@assert L>0
	@assert L == nv(graph)

    #basis_length = length(basis)
    sites_count = L
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

	L = sites(basis)

	@assert L>0
	@assert L == nv(graph)

	potential = Vector{SparseMatrixCSC{Float64,Int}}(undef, L)

	for i in 1:L
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
# This function is not used by the standard Hamiltonian creation anymore
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
matrix(bhh::BoseHubbardHamiltonian, parameters::Dict) = Array(spmatrix(bhh, parameters))
matrix(bhh::BoseHubbardHamiltonian, J, U) = Array(spmatrix(bhh, J, U))
matrix(bhh::BoseHubbardHamiltonian, J, eps, U) = Array(spmatrix(bhh, J, eps, U))


"""
Returns the sparse matrix representation of the given Bose-Hubbard
Hamiltonian
"""
function spmatrix(bhh::BoseHubbardHamiltonian, parameters::Dict)
	L = _get_L(parameters)
	J = parameters["J"]
	U = haskey(parameters, "U") ? parameters["U"] : zeros(L)
	U = length(U) == 1 ? [U for _ in 1:L] : U
	eps = haskey(parameters, "eps") ? parameters["eps"] : zeros(L)

	spmatrix(bhh, J, eps, U)
end

# TODO
# remove this function at some point
function spmatrix(bhh::BoseHubbardHamiltonian, J, U)
	L = bhh.L
	eps = zeros(L)
	return spmatrix(bhh, J, eps, U)
end

function spmatrix(bhh::BoseHubbardHamiltonian, J, eps, U)

    #N = bhh.N
    #k = bhh.k
    basis = bhh.basis
    tunnel = bhh.tunnel
    potential = bhh.potential

    #D = bose_hubbard_hilbert_space_size(k, N)
    D = length(basis)

    @assert length(J) == length(tunnel) "Number of tunnel parameters must be $(length(tunnel))"
	@assert length(eps) == length(potential) "Number of on-site potential parameters must be $(length(potential))"
    @assert length(U) == length(potential) "Number of on-site interaction parameters must be $(length(potential))"

    #H = similar(tunnel[1])
    #H .= 0
    H = spzeros(D,D)
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
