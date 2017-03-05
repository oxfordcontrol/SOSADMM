# SOSADMM
An open source first-order MATLAB solver for conic programs with row sparsity. SOSADMM implements the alternating direction method of multipliers (ADMM) described in our paper 
* Exploiting Sparsity of Coefficient Matching Conditions in Sum-of-Squares Programs using ADMM (included in the `doc/` folder)

In particular, SOSADMM exploits the sparsity of the coefficient matching conditions when SOS programs are formulated in the usual monomial basis to reduce the computational costs of the ADMM algorithm (see the detailed discriptions in our paper).

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
