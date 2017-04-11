
%% Example 1: Testing random unconstrained POP f(x)

%  To run this example, the following packages are required
%  1. GloptiPoly (Genertate the problem, essentail)
%  2. SeDuMi     (Interior-point method, for comparsion)
%  3. CDCS       (First-order method, for comparsion)
%  4. SCS        (First-order method, for comparsion, optional)

clear;

N = 2:2:10;   % number of variables
d = 2;        % half degree of the polynomial

TimePro = zeros(length(N),1);     % time for problem generation
TimeTotal = zeros(length(N),3);   % sedumi, sosadmm, cdcs(primal), scs(direct)
TimeSetup = zeros(length(N),2);   % sosadmm, cdcs(primal), scs(direct)
TimeADMM = zeros(length(N),2);
TimeAver = zeros(length(N),2);
Cost = zeros(length(N),3);
Iter = zeros(length(N),3);
Density = zeros(length(N),3);

data = cell(length(N),1);
%%
Maxiter = 1e3;
Tol     = 1e-4;

for i = 1:length(N)
    fprintf('testing : N = %i', N(i))
    %% using GloptiPloy 3 to construct random dense polynomials
    tic
    mpol('v',N(i),1);   
    b  = mmon(v,0,2*d-1);    % degree up tp 2d - 1;
    p0 = randn(1,length(b));  p0 = p0/norm(p0);
    p = p0*b + sum(mmon(v,d).^2);
    
    %% Problem data in Sedumi form
    P = msdp(min(p));                
    [A,b,c,K] = msedumi(P);
    [m,n] = size(A);
    Density(i,:) = [m,n,sum(sum(spones(A)))/m/n];
    
    TimePro(i) = toc;
    data{i}.A = A;
    data{i}.b = b;
    data{i}.c = c;
    data{i}.K = K;
    %% solutions using different method
    x = zeros(length(c),1);
    [x,y,info] = sedumi(A,b,c,K);

    
    % by sosadmm -- exploiting row sparsity
    opts.Max_iter = Maxiter;
    opts.eps      = Tol;
    [x1,y1,z1,info1] = sosadmm(A',b,c,K,opts);
    
    % by cdcs - primal
    opts.relTol = Tol;
    opts.solver = 'primal';
    opts.maxIter = Maxiter;
    [x2,y2,z2,info2] = cdcs(A',b,c,K,opts);
    
%     % by SCS
%     params.max_iters = Maxiter;
%     params.eps = Tol;
%     [x3,y3,cscs,info3] = solveWithSCSdirect(A',full(b),full(c),K,params);
    
    %% statistics
    TimeTotal(i,:) = [info.wallsec,info1.time.total,info2.time.total];   
    TimeSetup(i,:) = [info1.time.data,info2.time.setup]; 
    TimeADMM(i,:) = [info1.time.admm,info2.time.admm]; 
    Cost(i,:) = [c'*x,c'*x1,c'*x2];
    Iter(i,:) = [info.iter,info1.iter,info2.iter];
    TimeAver(i,:)  = TimeADMM(i,:)./Iter(i,2:end);
end
