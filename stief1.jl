using LinearAlgebra
function stief1(X, Z)
#
# Syntax: Q = stief1(X, Z)
#
# Inputs: X is a K x N matrix and Z is a M x N matrix, where M >= K
#
# Outputs: Q is an M x K matrix with orthonormal columns,
# i.e., a matrix in the Stiefel manifold V_K(C^M),
# that minimizes over Q \sum_n |Q X[:,n] - Z[:,n]|_2^2
# This is a generalized version of the Procrustes problem.
    U,s,V = svd(Z*X');
    return U*V'
end
