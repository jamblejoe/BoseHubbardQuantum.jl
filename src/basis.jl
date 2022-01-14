##############################################################################
#
# Basis
#
##############################################################################

"""
For given number of particles N and site count k returns the size of the Hilbert space
The Hilbert space size is binom(N+k-1, k-1).
"""
function bose_hubbard_hilbert_space_size(k::Integer, N::Integer)
  return binomial(N+k-1,k-1)
end

function bose_hubbard_hilbert_space_size(k::Integer, Nmin::Integer, Nmax::Integer)
    Nmin != 0 && throw(DomainError(Nmin, "Only implemented for Nmin=0!"))
    return binomial(Nmax+k,k)
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
    N::Int      # keep this for convinience and to not break old code
    Nmin::Int
    Nmax::Int
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

LtrAscBasis(k::Integer, N::Integer) = LtrAscBasis(k, N, N)
function LtrAscBasis(k::Integer, Nmin::Integer, Nmax::Integer)

    @assert k>0
    @assert Nmax >= Nmin >=0
    @assert typeof(Nmin) == typeof(Nmax)

    # create the basis
    state = zeros(typeof(Nmax), k)
    basis = typeof(state)[]

    for N in Nmin:Nmax
        ltr_asc_loop!(basis, state, N, 1)
    end
    # calculate the index
    index = create_basis_index(basis)

    return LtrAscBasis(k, Nmax, Nmin, Nmax, basis, index)
end

# can this be improved by making it non-recursive?
# see:
# https://github.com/georglind/bosehubbard/blob/master/bosehubbard.py
# http://iopscience.iop.org/article/10.1088/0143-0807/31/3/016
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

Base.length(basis::LtrAscBasis) = length(basis.basis)
sites(basis::LtrAscBasis) = basis.k
Base.iterate(basis::LtrAscBasis) = iterate(basis.basis)
Base.iterate(basis::LtrAscBasis, state) = iterate(basis.basis, state)
Base.isequal(b1::LtrAscBasis, b2::LtrAscBasis) = b1.k == b2.k && b1.N == b2.N

getposition(basis::LtrAscBasis, state::AbstractVector) = basis.index[state]
getstate(basis::LtrAscBasis, index::Integer) = basis.basis[index]
getstate!(state::AbstractVector, basis::LtrAscBasis, index::Integer) = state .= basis.basis[index]
Base.in(state::AbstractVector, basis::LtrAscBasis) = haskey(basis.index, state)
Base.eachindex(basis::LtrAscBasis) = 1:length(basis)


##############################################################################

"""
https://arxiv.org/pdf/1410.7280.pdf
"""
struct PonomarevBasis <: AbstractBasis
    k::Int
    N::Int
    D::Int # Hilbert space size
    Ds::Matrix{Int}
end

function PonomarevBasis(k::Integer, N::Integer)

    Ds = Matrix{Int}(undef, k+1, N+1)
    for i in 0:k
        for j in 0:N
            Ds[i+1,j+1] = bose_hubbard_hilbert_space_size(i,j)
        end
    end

    PonomarevBasis(k, N, bose_hubbard_hilbert_space_size(k,N), Ds)

end

PonomarevBasis(parameters::Dict) = PonomarevBasis(parameters["k"], parameters["N"])

@inline _D(basis::PonomarevBasis, k::Integer, N::Integer) = basis.Ds[k+1, N+1]
@inline _D_unsafe(basis::PonomarevBasis, k::Integer, N::Integer) = @inbounds basis.Ds[k+1, N+1]

Base.length(basis::PonomarevBasis) = basis.D
sites(basis::PonomarevBasis) = basis.k
Base.isequal(b1::PonomarevBasis, b2::PonomarevBasis) = b1.k == b2.k && b1.N == b2.N

"""
Returns the position of a state in the basis.
"""
function getposition(basis::PonomarevBasis, state::AbstractVector)
    k = basis.k
    N = basis.N

    length(state) == k || throw(DomainError(k, "length of state must be $k"))
    #all(0 .<= state .<= N) || throw(DomainError(N, "components of state must be between 0 and $N"))
    for s in state
        0 <= s <= N || throw(DomainError(N, "components of state must be between 0 and $N"))
    end
    sum(state) == N || throw(DomainError(N, "components of state must sum to $N"))

    index = 1
    l = N
    for (m_l,n_i) in enumerate(state)
        for _ in 1:n_i
            #index += bose_hubbard_hilbert_space_size(k-m_l, l)
            index += _D_unsafe(basis, k-m_l, l)
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

function getstate!(state::AbstractVector{T}, basis::PonomarevBasis, index::Integer) where T
    k = basis.k
    N = basis.N
    #Ds = basis.Ds

    one(index) <= index <= length(basis) || error("index $index out of bounds [1,$(length(basis))]")

    state .= zero(T)

    for i in N:-1:1
        m_i = 0
        #while bose_hubbard_hilbert_space_size(m_i+1,i) < index
        while _D_unsafe(basis, m_i+1, i) < index
            m_i += one(index)
        end
        state[k-m_i] += one(T)
        #index -= bose_hubbard_hilbert_space_size(m_i,i)
        index -= _D_unsafe(basis, m_i, i)
    
    end
    state
end


"""
Returns the basis element at position index.
"""
Base.getindex(basis::PonomarevBasis, index::Integer) = getstate(basis, index)


function Base.iterate(basis::PonomarevBasis, state=1)
    count = state
    count > length(basis) ? nothing : ( getstate(basis, count), count+1 )
end

Base.eachindex(basis::PonomarevBasis) = 1:length(basis)

function Base.in(state::AbstractVector, basis::PonomarevBasis)
    k = basis.k
    N = basis.N

    length(state) == k || return false
    #all(0 .<= state) || return false
    for s in state
        0 <= s <= N || return false
    end
    sum(state) == N || return false
    return true
end

##############################################################################


struct LtrAscCutoffBasis <: AbstractBasis
	k::Int
	d::Int
end
"""
Basis which has no restriction on number of particles but to make it finite
there is a cutoff d on each site, i.e. on each site are maximal d particles.
The cutoff makes sense for positive U, as a lot of particles on one site are
energetically punished. I think this makes the most sense when some product of
cutoff d and U is big enough.
"""
LtrAscCutoffBasis(parameters::Dict) = LtrAscCutoffBasis(parameters["k"], parameters["d"])

Base.length(basis::LtrAscCutoffBasis) = (basis.d+1)^basis.k
sites(basis::LtrAscCutoffBasis) = basis.k
#=
This probably can be made non-allocating
=#
function Base.iterate(basis::LtrAscCutoffBasis, state=1)
	count = state
	count > length(basis) ? nothing : ( getstate(basis, count), count+1 )
end
Base.isequal(b1::LtrAscCutoffBasis, b2::LtrAscCutoffBasis) = b1.k == b2.k && b1.d == b2.d

function getposition(basis::LtrAscCutoffBasis, state::AbstractVector{<:Integer})
	k = basis.k
	d = basis.d

	#@assert length(state) == k
	#@assert all(0 .<= state .<= d)
	state in basis || error("state $state not in basis")
	index = 1
	for i in eachindex(state)
		index += (d+1)^(k-i) * state[i]
	end
	index
end

"""
Returns the basis element at position index.
"""
Base.getindex(basis::LtrAscCutoffBasis, index::Integer) = getstate(basis, index)


"""
Returns state of the basis with position index.
"""
function getstate(basis::LtrAscCutoffBasis, index::Integer)
    state = zeros(Int, basis.k)
    getstate!(state, basis, index)
end

function getstate!(state::AbstractVector, basis::LtrAscCutoffBasis, index::Integer)
    k = basis.k
	d = basis.d

	# do we really need the following condition?
	length(state) != k && error("length(state) != k. Got $(length(state)) and $k")

    !(1 <= index <= length(basis)) && error("Index $index out of bounds [1,$(length(basis))]")

	digits!(state, index-1, base=d+1)
	reverse!(state)
	state
end

Base.in(state::AbstractVector{<:Integer}, basis::LtrAscCutoffBasis) = length(state)==basis.k && all(0 .<= state .<= basis.d)
Base.eachindex(basis::LtrAscCutoffBasis) = 1:length(basis)


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
