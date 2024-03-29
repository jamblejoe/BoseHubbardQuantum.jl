module BoseHubbardQuantum

#=
    TODO:
        - implement on the fly Hamiltonian creation for big systems
        - make basis creation multithreaded (req. Julia 1.3+)

=#
using Reexport
@reexport using LinearAlgebra
@reexport using Graphs
@reexport using SparseArrays


# basis
export bose_hubbard_hilbert_space_size
export AbstractBasis, get_or_create_basis
export LtrAscBasis, PonomarevBasis, LtrAscBasis0N, LtrAscCutoffBasis
export getposition, getstate, getstate!

# graphs
export dimer_graph, trimer_graph, chain_graph

# operators
export tunnel_operator, tunnel_symm_operator, number_operator, kinetic_operator
export tunnel_spmatrix, tunnel_symm_spmatrix, number_spmatrix, kinetic_spmatrix
export creation_spmatrix, annihilation_spmatrix

# bosehubbard hamiltonian
export BoseHubbardHamiltonian, BoseHubbardHamiltonianChain
export matrix, spmatrix

# partial traces
export partialtr, partialtr!

# depricated
export k_chain_graph

_get_L(parameters::Dict) = haskey(parameters, "k") ? parameters["k"] : parameters["L"]



include("bose_hubbard.jl")
include("basis.jl")
include("operator.jl")
include("partial_trace.jl")
include("deprecated.jl")

end
