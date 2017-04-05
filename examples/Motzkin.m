
clc;clear
%% Motzkin polynomial 
% f(X,Y) = X^4Y^2 + X^2Y^4 ? 3*X^2Y^2 + 1
% This function is non-negative but not an SOS
% The performance of our algorithm is not good

Maxiter = 1e3;
Tol     = 1e-4;

mpol x 2;   % create a 2x1 vector x of variables x(1),x(2)   
p = x(1)^4*x(2)^2 + x(1)^2*x(2)^4 - 3 * x(1)^2*x(2)^2 + 1; % Motzkin polynomial 
    
P = msdp(min(p));         
[A,b,c,K] = msedumi(P);
[m,n] = size(A);
Density = [m,n,sum(sum(spones(A)))/m/n];

%% solution by sedumi
[x,y,info] = sedumi(A,b,c,K);
 
%% by sosadmm -- exploiting row sparsity
opts.Max_iter = Maxiter;
opts.eps      = Tol;
[x1,y1,z1,info1] = sosadmm(A',b,c,K,opts);
    
%% by cdcs - primal
opts.relTol = Tol;
opts.solver = 'primal';
opts.maxIter = Maxiter;
[x2,y2,z2,info2] = cdcs(A',b,c,K,opts);

Cost = [c'*x, c'*x1, c'*x2]
Time = [info.wallsec,info1.time.total,info2.time.total]

Error = (Cost - Cost(1))./Cost(1)