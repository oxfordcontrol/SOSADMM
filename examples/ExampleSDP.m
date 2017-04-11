
%% Example 1: Testing random SDP

clear;

n = 300;
density = 0.0001;
m = 2000:2000:20000;


TimeTotal = zeros(length(m),3);   % sedumi, sosadmm, cdcs(primal), scs(direct)
TimeSetup = zeros(length(m),3);   % sosadmm, cdcs(primal), scs(direct)
TimeADMM = zeros(length(m),3);
TimeAver = zeros(length(m),3);
Cost = zeros(length(m),3);
Iter = zeros(length(m),3);

%%
Maxiter = 2e3;
Tol     = 1e-4;

for i = 1:length(m)
    
    %% construct random row sparse SDP
        tic
        [At,b,c,K] = RowSpaSDP(m(i),n,density);
        A = At';
        %% solutions using different method
        x = zeros(length(c),1);
        if m(i) < 6000     %% only test sedumi if it took < 1800 for last instancem
            [x,y,info] = sedumi(A,b,c,K);
        end

        % by sosadmm -- exploiting row sparsity
        opts.Max_iter = Maxiter;
        opts.eps      = Tol;
        [x1,y1,z1,info1] = sosadmm(A',b,c,K,opts);

        % by cdcs - primal
        opts.relTol = Tol;
        opts.solver = 'primal';
        opts.maxIter = Maxiter;
        [x2,y2,z2,info2] = cdcs(A',b,c,K,opts);

        % by SCS
        params.max_iters = Maxiter;
        params.eps = Tol;
        [x3,y3,cscs,info3] = solveWithSCSdirect(A',full(b),full(c),K,params);
        
        %% statistics
        TimeTotal(i,:) = [info1.time.total,info2.time.total,(info3.solveTime+info3.setupTime)/1e3];   
        TimeSetup(i,:) = [info1.time.data,info2.time.setup,info3.setupTime/1e3]; 
        TimeADMM(i,:)  = [info1.time.admm,info2.time.admm,info3.solveTime/1e3]; 
        Cost(i,:)  = [c'*x1,c'*x2,cscs];
        Iter(i,:)  = [info1.iter,info2.iter,info3.iter];
        TimeAver(i,:)  = TimeADMM(i,:)./Iter(i,:);

end




