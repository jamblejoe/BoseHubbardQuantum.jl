module BoseHubbardQuantum


using LightGraphs
using SparseArrays


# basis
export bose_hubbard_hilbert_space_size
export AbstractBasis, LtrAscBasis, get_or_create_basis

# graphs
export dimer_graph, trimer_graph, k_chain_graph

# operators
export tunnel_operator, tunnel_symm_operator, number_operator, kinetic_operator
export tunnel_spmatrix, tunnel_symm_spmatrix, number_spmatrix, kinetic_spmatrix

# bosehubbard hamiltonian
export BoseHubbardHamiltonian, BoseHubbardHamiltonianChain
export matrix, spmatrix


include("bose_hubbard.jl")

end
