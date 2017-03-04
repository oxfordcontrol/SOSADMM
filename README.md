# SOSADMM
An open source first-order MATLAB solver for conic programs with row sparsity.

## Description<a name="Description"></a>

SOSADMM solves in the standard SDP in the primal vectorized form

		minimize 	c'x						
	(1)	subject to	Ax = b,					
					x \in K							

SOSADMM supports cartesian products of the following cones:

* R^n (free variables)
* Non-negative orthant
* Second-order cone
* Positive semidefinite cone

SOSADMM is called with the syntax

	>> [x,y,z,info] = sosadmm(At,b,c,K,options);
	
where `At` is the transpose of the matrix `A` in problem (1) above. 
Note that the inputs and outputs are in the same format used by SeDuMi. 

