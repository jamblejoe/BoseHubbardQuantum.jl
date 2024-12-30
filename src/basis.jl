##############################################################################
#
# Basis
#
##############################################################################

"""
For given number of particles N and site count L returns the size of the Hilbert space
The Hilbert space size is binom(N+L-1, L-1).
"""
function bose_hubbard_hilbert_space_size(L::Integer, N::Integer)
  return binomial(N+L-1,L-1)
end

function bose_hubbard_hilbert_space_size(L::Integer, Nmin::Integer, Nmax::Integer)
    Nmin != 0 && throw(DomainError(Nmin, "Only implemented for Nmin=0!"))
    return binomial(Nmax+L,L)
end


# Abstract Basis is not implemented yet and is just a placeholder
# one should think of AbstractBasis needing to be abstract type
# or needing to have fields
abstract type AbstractBasis end

##############################################################################


struct LtrAscBasis <: AbstractBasis
    L::Int
    Nmin::Int
    Nmax::Int
    states::Vector{Vector{Int}}
    index::Dict{Vector{Int}, Int}
end


"""
    Creates a vector of basis states for N particles and L sites.
    The states are ordered in ascending order from left to right, e.g.
    for 2 particles and 3 sites
    (0,0,2) -> (0,1,1) -> (0,2,0) -> (1,0,1) -> (1,1,0) -> (2,0,0).
    The total number of basis states is binomial(N+L-1,L-1).
"""
function LtrAscBasis end
LtrAscBasis(L::Integer, N::Integer) = LtrAscBasis(L, N, N)
function LtrAscBasis(L::Integer, Nmin::Integer, Nmax::Integer)

    if L <= 0
        throw(ArgumentError("L must be greater than 0"))
    end
    if Nmax < Nmin || Nmin < 0
        throw(ArgumentError("Nmax must be greater than or equal to Nmin and Nmin must be non-negative"))
    end

    # create the basis
    tmp_state = zeros(typeof(Nmax), L)
    states = typeof(tmp_state)[]

    for N in Nmin:Nmax
        ltr_asc_loop!(states, tmp_state, N, 1)
    end
    # calculate the index
    index = create_basis_index(states)

    return LtrAscBasis(L, Nmin, Nmax, states, index)
end

# can this be improved by making it non-recursive?
# see:
# https://github.com/georglind/bosehubbard/blob/master/bosehubbard.py
# http://iopscience.iop.org/article/10.1088/0143-0807/31/3/016
function ltr_asc_loop!(states, state, n::Integer, pos::Integer)
    L = length(state)
    if pos < L
        for i in 0:n
            state[pos] = i
            ltr_asc_loop!(states, state, n-i, pos+1)
        end
    else
        state[pos] = n
        push!(states, copy(state))
    end
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

Base.length(basis::LtrAscBasis) = length(basis.states)
sites(basis::LtrAscBasis) = basis.L
Base.iterate(basis::LtrAscBasis) = iterate(basis.states)
Base.iterate(basis::LtrAscBasis, state) = iterate(basis.states, state)
Base.isequal(b1::LtrAscBasis, b2::LtrAscBasis) = b1.L == b2.L && b1.N == b2.N

getposition(basis::LtrAscBasis, state::AbstractVector) = basis.index[state]
getstate(basis::LtrAscBasis, index::Integer) = basis.states[index]
getstate!(state::AbstractVector, basis::LtrAscBasis, index::Integer) = state .= basis.states[index]
Base.in(state::AbstractVector, basis::LtrAscBasis) = haskey(basis.index, state)
Base.eachindex(basis::LtrAscBasis) = 1:length(basis)


##############################################################################

"""
    Ponomarev basis for the Bose-Hubbard model. 
    This is a non-allocating version of the LtrAscBasis for N=Nmin=Nmax.
    Based on the paper by Raventós et al. (2014) https://arxiv.org/pdf/1410.7280.pdf
"""
struct PonomarevBasis <: AbstractBasis
    L::Int
    N::Int
    D::Int # Hilbert space size
    Ds::Matrix{Int}
end

function PonomarevBasis(L::Integer, N::Integer)

    Ds = Matrix{Int}(undef, L+1, N+1)
    for i in 0:L
        for j in 0:N
            Ds[i+1,j+1] = binomial(j+i-1,i-1)
        end
    end

    PonomarevBasis(L, N, binomial(N+L-1,L-1), Ds)

end


@inline _D(basis::PonomarevBasis, L::Integer, N::Integer) = basis.Ds[L+1, N+1]
@inline _D_unsafe(basis::PonomarevBasis, L::Integer, N::Integer) = @inbounds basis.Ds[L+1, N+1]

Base.length(basis::PonomarevBasis) = basis.D
sites(basis::PonomarevBasis) = basis.L
particles(basis::PonomarevBasis) = basis.N
Base.isequal(b1::PonomarevBasis, b2::PonomarevBasis) = b1.L == b2.L && b1.N == b2.N

"""
Returns the position of a state in the basis.
"""
function getposition(basis::PonomarevBasis, state::AbstractVector)
    L = basis.L
    N = basis.N

    length(state) == L || throw(DomainError(L, "length of state must be $L"))
    for s in state
        0 <= s <= N || throw(DomainError(N, "components of state must be between 0 and $N"))
    end
    sum(state) == N || throw(DomainError(N, "components of state must sum to $N"))

    index = 1
    l = N
    for (m_l,n_i) in enumerate(state)
        for _ in 1:n_i
            index += _D_unsafe(basis, L-m_l, l)
            l -= 1
        end
    end
    index
end

"""
Returns state of the basis with position index.
"""
function getstate(basis::PonomarevBasis, index::Integer)
    state = zeros(Int, basis.L)
    getstate!(state, basis, index)
end

function getstate!(state::AbstractVector{T}, basis::PonomarevBasis, index::Integer) where T
    L = basis.L
    N = basis.N

    one(index) <= index <= length(basis) || error("index $index out of bounds [1,$(length(basis))]")

    state .= zero(T)

    for i in N:-1:1
        m_i = 0
        #while bose_hubbard_hilbert_space_size(m_i+1,i) < index
        while _D_unsafe(basis, m_i+1, i) < index
            m_i += one(index)
        end
        state[L-m_i] += one(T)
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
    L = basis.L
    N = basis.N

    length(state) == L || return false
    for s in state
        0 <= s <= N || return false
    end
    sum(state) == N || return false
    return true
end

##############################################################################

"""
Basis which has no restriction on number of particles but to make it finite
there is a cutoff d on each site, i.e. on each site are maximal d particles.
The cutoff makes sense for positive U, as a lot of particles on one site are
energetically punished. I think this makes the most sense when some product of
cutoff d and U is big enough.
"""
struct LtrAscCutoffBasis <: AbstractBasis
	L::Int
	d::Int
end



Base.length(basis::LtrAscCutoffBasis) = (basis.d+1)^basis.L
sites(basis::LtrAscCutoffBasis) = basis.L
#=
This probably can be made non-allocating
=#
function Base.iterate(basis::LtrAscCutoffBasis, state=1)
	count = state
	count > length(basis) ? nothing : ( getstate(basis, count), count+1 )
end
Base.isequal(b1::LtrAscCutoffBasis, b2::LtrAscCutoffBasis) = b1.L == b2.L && b1.d == b2.d

function getposition(basis::LtrAscCutoffBasis, state::AbstractVector{<:Integer})
	L = basis.L
	d = basis.d

	#@assert length(state) == L
	#@assert all(0 .<= state .<= d)
	state in basis || error("state $state not in basis")
	index = 1
	for i in eachindex(state)
		index += (d+1)^(L-i) * state[i]
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
    state = zeros(Int, basis.L)
    getstate!(state, basis, index)
end

function getstate!(state::AbstractVector, basis::LtrAscCutoffBasis, index::Integer)
    L = basis.L
	d = basis.d

	# do we really need the following condition?
	length(state) != L && error("length(state) != L. Got $(length(state)) and $L")

    (1 <= index <= length(basis)) || error("Index $index out of bounds [1,$(length(basis))]")

	digits!(state, index-1, base=d+1)
	reverse!(state)
	state
end

Base.in(state::AbstractVector{<:Integer}, basis::LtrAscCutoffBasis) = length(state)==basis.L && all(0 .<= state .<= basis.d)
Base.eachindex(basis::LtrAscCutoffBasis) = 1:length(basis)


##############################################################################



