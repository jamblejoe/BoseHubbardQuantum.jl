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
# exponentials of Bose-Hubbard Hamiltonian
#
##############################################################################

exp_scaled(c::Number, F::Eigen) = exp_scaled(c, F.values, F.vectors)

function exp_scaled(c::Number, eigenvalues::AbstractVector, eigenstates::AbstractMatrix)
	#exp!(similar(eigenstates), c, eigenvalues, eigenstates)
	eigenstates * Diagonal(exp.(c .* eigenvalues)) * eigenstates'
end



#function myexp!(res::AbstractMatrix, c::Number, eigenvalues::AbstractVector, eigenstates::AbstractMatrix)
	#res .= eigenstates'
	#res .*= c .* eigenvalues
	#mul!(res, )
	#mul!(res, eigenstates, (exp.(c .* eigenvalues) .* eigenstates'))
#	res .= eigenstates * Diagonal(exp.(c .* eigenvalues)) * eigenstates'
#	res
#end
