

%% The classical problem of minimizing globally the two-dimensional
%  six-hump camel back function
%  f(x1,x2) = x_1^2(4 - 2.1 x_1^2 + x_1^4/3) + x_1x_2+x_2^2(-4+4x_2^2)

%  To run this example, the following packages are required
%  1. GloptiPoly (Genertate the problem)
%  2. SeDuMi     (Interior-point method, for comparsion)
%  3. CDCS       (First-order method, for comparsion)

clc;clear
Maxiter = 1e3;
Tol     = 1e-4;

mpol x 2;   % create a 2x1 vector x of variables x(1),x(2)   
p = x(1)^2*(4 - 2.1*x(1)^2 + x(1)^4/3) + x(1)*x(2)+x(2)^2*(-4+4*x(2)^2); % hump camel function
    
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

%% statistics
Cost = [c'*x, c'*x1, c'*x2]
Time = [info.wallsec,info1.time.total,info2.time.total]
Error = (Cost - Cost(1))./Cost(1)

    