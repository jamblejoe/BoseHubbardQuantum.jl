"""
Created on 10 Sep 2019


TODO:
    - add more tests for LtrAscBasis0N combined with operators
"""

using BoseHubbardQuantum
using Test
using SparseArrays
using LinearAlgebra

# partialtr
using Random

import BoseHubbardQuantum: tunnel_spmatrices

include("basis_tests.jl")
include("operator_tests.jl")
include("hamiltonian_tests.jl")
include("partialtr_tests.jl")
