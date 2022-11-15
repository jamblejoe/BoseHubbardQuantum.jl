"""
Creates a chain graph of length k
"""
function k_chain_graph(k::Integer)
    return path_graph(k)
end