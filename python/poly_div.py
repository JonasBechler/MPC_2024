# function [ R, S ] = poly_div( Q, P, d )
# % Calculates the polynomials R and S solving
# %
# % Q / P = R + q^{-d} S / P
# %
# % by solving the diophantic equation
# %
# %  Q = P*R + q^{-d} S
# %
# % using a set of linear equations to perform a comparison 
# % of the coefficients
# 
# if d == 0 
#     R = 0;
#     S = Q;
#     return;
# end
# 
# 
# % degrees of q and p
# nq = length(Q)-1; 
# np = length(P)-1;
# 
# %First the degrees of the polynomials are determined
# %
# %  the degree of R is 1 less the number of division operations
# %  that needs to be performed which is d
# 
# nr = d - 1;
# 
# % The degree of the remainder depends on the degrees of Q and P
# % and the number of divisions that needs to be performed according to d.
# % If the degree of divident Q is higher than the degree of the divisor
# % P, then it takes nq-np of division steps until the degrees
# % of the Divident and the remainder are equal. If d is even higher than
# % that number, the degree of the remainder increases with each step. 
# 
# ns = max(nq,np+d-1);
# 
# % The number of coefficients to be determined is
# nc = nr+1 + ns-d+1; % in S th efirst d coefficients are 0 -> ns-d+1 parameters
# 
# %in case of nq <= np+d-1 the number of coefficients in S
# % remains constant equal to np
# 
# % The linear equations are put in matrix form
# % q = A*[R;S] 
# % with 
# A = zeros(nc,nc); % reserve memory
# q = [Q,zeros(1,nc-nq-1)]'; % -1 because nq is degree and nc is number of coef.
# % each element in q represents a power of q
# % and so does every line in A
# % The comparison of coefficients leads to the band structure of A
# for i = 1:nr+1
#     A(i:i+np,i) = P';
# end
# % In this case the coefficients for S are all multiplied by one
# % and the polynom starts with q^{-d}, therefore the coefficients 
# % 1 start at line d+1 and create a partion of an identity matrix
# A(d+1:ns+1,nr+2:nc) = eye(ns-d+1);
# x = A \ q; % better than inv
# R = x(1:nr+1)';
# S = x(nr+2:end)';
# 
# %
# %Test
# 
# % PR = conv(P,R); 
# % q_dS = [zeros(1,d),S]; % zero padding for degree q^{-d}
# % 
# % if np+nr >= d+ns
# %    TEST_Q = PR;
# %    TEST_Q(1:d+ns) =  TEST_Q(1:d+ns) + q_dS
# % else
# %    TEST_Q = q_dS;
# %    TEST_Q(1:np+nr+1) =  TEST_Q(1:np+nr+1) + PR
# % end  
# 
# end


import numpy as np
import control as ct
 
def  poly_div( Q, P, d ):
    if d == 0:
        R = 0
        S = Q
        return R, S
    n_q = len(Q) - 1
    n_p = len(P) - 1
    n_r = d - 1
    n_s = max(n_q, n_p + d - 1)
    n_c = n_r + 1 + n_s - d + 1
    A = np.zeros((n_c, n_c))
    q = np.concatenate((Q, np.zeros(n_c - n_q - 1)))
    for i in range(n_r + 1):
        A[i:i + n_p + 1, i] = P
    A[d:d + n_s + 1, n_r + 1:n_c] = np.eye(n_s - d + 1)
    x = np.linalg.solve(A, q)
    R = x[0:n_r + 1]
    S = x[n_r + 1:]
    return R, S