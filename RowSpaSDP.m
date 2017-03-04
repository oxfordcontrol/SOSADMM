function [At,b,c,K] = RowSpaSDP(m,n,density)

% Setup problem for SDPs with row sparsity. Inputs:
% m: number of equality constraints
% n: size of SDP cones
% density: density of each row constraint

%% Data 
At = [];
for i = 1:m
    Amat = full(sprandsym(n,density));   %% random symmetric sparse matrices,
    At = [At, sparse(vec( (Amat+Amat')/2 ))];
end


%% stict primal feasible point
Temp = rand(n);%.*Spa;
Temp = (Temp+Temp')/2;
X = Temp + (-min(eig(Temp))+1)*eye(n);  
b = At'*vec(X);    

%% strict dual feasible point 
y = rand(m,1);
Temp = rand(n);%*Spa;
Temp = (Temp+Temp')/2;
S = Temp + (-min(eig(Temp))+1)*eye(n);
c = vec(S) + At*y;

%% cone
K.f = 0;
K.l = 0;
K.q = 0;
K.s = n; 





