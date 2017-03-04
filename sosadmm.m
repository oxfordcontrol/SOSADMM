function [x,y,cost,info] = sosadmm(At,b,c,K,options)
%  Solving the SDP formulation for SOS polynomials using ADMM approach
%  min  c^Tx                        (primal SDP)
%  s.t. Ax = b
%       x \in K
%  [x,y,cost,info] = sosadmm(At,b,c,K,options)


%% Default parameters
opts.Max_iter  = 2e3;          % maximum steps
opts.eps       = 1.e-4;        % relative tolerances for the stopping criterion of ADMM
opts.rho       = 1;            % penalty 
opts.adaptive  = true;         % adaptive penalty factor?
opts.tau       = 2;            % increase factor for adaptive penalty scheme (must be > 1)
opts.rhoMin    = 1e-2;         % min value for rho
opts.rhoMax    = 1e3;          % max value for rho
opts.rhoIt     = 50;           % if the ratio of pinf and dinf is on the same side of eta for rhoIt iterations, adjust rho by tau
opts.mu        = 10;           % ratio of residuals for adaptive penalty scheme
opts.dispIter  = 50;           % for output information printing
opts.verbose   = 1;

% Set user options
if nargin > 4
    opts = setstructfields(opts,options);
end

%% Parameters and initial checks
[At,b,c,K,opts]  = checkInputs(At,b,c,K,opts);
%[At,b,c,K,opts] = rescaleData2(At,b,c,K,opts);
[At,c] = svecData(At,c,K);
[n,m]  = size(At);

%% =======================================================================
%    Main Function: data preprocessing + iterations of ADMM
%  =======================================================================
myline = repmat('=',1,70);
myline1 = repmat('-',1,70);
fprintf([myline,'\n'])
fprintf('Solving row-sparse SDPs via ADMM v1.0\n')
fprintf([myline,'\n'])

%%   Data Preprocessing 
tdata = tic;
[Hind,invX,nnzR,nzA,sbA,normAi,nzAind]  = preData(At);
time.data = toc(tdata);
fprintf(['Input data processed in ' num2str(time.data) ' seconds.\n'])
fprintf('Free variables         : %i                \n',K.f);
fprintf('Non-negative variables : %i                \n',K.l);
fprintf('Second-order cones     : %i (max. size: %i)\n',length(K.q),max(K.q));
fprintf('Semidefinite cones     : %i (max. size: %i)\n',length(K.s),max(K.s));
fprintf('Affine constraints     : %i                \n',m);
fprintf('Non-zero elements      : %i                \n',sum(nnzR));
fprintf('Non-zero density       :%10.3e            \n',sum(nnzR)/m/n);
fprintf([myline1,'\n'])

%% Iterations of ADMM
% initialize vectorized versions of all variables
x  = zeros(n,1);
z  = zeros(n,1);
zt = zeros(sum(nnzR),1); %% in a vector
mu = zeros(sum(nnzR),1);
xi = zeros(n,1);

% initialize blockified versions of all variables
[X,Xi] = makeConeVariables(K);

dresi = zeros(opts.Max_iter,1);       %% dual residual
presi = zeros(opts.Max_iter,1);       %% primal residual
pcost = zeros(opts.Max_iter,1);       %% primal cost
dcost = zeros(opts.Max_iter,1);       %% dual cost
gap   = zeros(opts.Max_iter,1);       %% duality gap

%% output information
linetitle = ' iter |   presi   |   dresi   |   cost(c*x) |    rho    |   time (s) \n';
fprintf(linetitle)

admmtime = tic;
for iter = 1 : opts.Max_iter 
    
    %% X-minimization step
    [x_n,X_n] = Xblockmin(z,zt,xi,mu);      
    
    %% Y-minimization step
    [z_n,zt_n,Omega,xtmp] = Yblockmin(X_n,x_n,Xi,xi,mu);
    
    %% dual update step
    [mu_n,xi_n,Xi_n] = Dupdate(mu,xi,Xi,z_n,zt_n,x_n,xtmp);   
       
    %% calculate  residual
    [isConverged,pres,dres] = ConverCheck(z_n,x_n,zt_n,z,zt,xtmp,mu_n,xi_n);
   
    %dcost(iter) = -b'*Omega;     %% dual cost
    pcost(iter) = c'*x_n;        %% primal cost
    presi(iter) = pres;
    dresi(iter) = dres;
    %gap(iter)   = abs(pcost(iter) - dcost(iter))/(1 + abs(pcost(iter)) + abs(dcost(iter)));
       
    %% whether to stop && output information 
    if opts.verbose && (iter == 1 || ~mod(iter,opts.dispIter) || isConverged)
        fprintf('%5d |  %8.2e |  %8.2e |  %10.3e |  %8.2e |   %8.2e\n',...
            iter,pres,dres,pcost(iter),opts.rho, toc(admmtime))
    end    
    if(isConverged)
        break;
    end
    %% iteration of next step
    x  = x_n;  X  = X_n;
    z  = z_n;  zt = zt_n;
    mu = mu_n; xi = xi_n; Xi = Xi_n;
end

Xvec = cellfun(@(x)vec(x),X_n,'UniformOutput',false);
x = vertcat(Xvec{:});
y = Omega;
time.admm = toc(admmtime);
time.total= time.admm + time.data;
cost      = pcost(iter);
info.iter = iter;
info.dual = dresi(1:iter);
info.time = time;
info.pri  = presi(1:iter);

%% Solution summary
fprintf([myline1,'\n'])
fprintf('SOLUTION SUMMARY\n')
fprintf('Number of iterations: %11.d\n',info.iter)
fprintf('Primal residual     : %11.4e\n',info.pri(iter))
fprintf('Dual residual       : %11.4e\n',info.dual(iter))
fprintf('Optimal value       : %11.4e\n',cost)
fprintf('Data processing time: %11.4e\n',info.time.data)
fprintf('Admm time           : %11.4e\n',info.time.admm)
fprintf('Total time          : %11.4e\n',info.time.total)
fprintf([myline,'\n'])


    %% ======================================================================= %
    %  NESTED FUNCTIONS
    %  ======================================================================= %

    %% data pre-processing
    function [Hind,invX,nnzR,nzA,sbA,normAi,nzAind] = preData(At)
        %% At \in n * m; m : the number of constraints 
        %[n,m] = size(At);
        nnzR = full(sum(spones(At)));  %% the number of non-zeros of each constraint in At
        [Hind,cols] = find(At);
        H = accumarray(Hind,1);
        if length(H) < n
            H(length(H)+1:n) = 0;      %% in case of inconsistent dimension
        end
        invX = 1./(ones(n,1) + H);
        
        tmp = At(:);
        nzA = full(tmp(find(At)));     %% the non-zero elements stacked in a vector
        
        %% sparse block diagnal form       
        inzA = cols; jnzA = [1:length(nzA)]';
        sbA  = sparse(inzA,jnzA,nzA,m,length(nzA));
        
        %% norm of each row
        normAi = zeros(m,1);  
        nzAind = [0,cumsum(nnzR)];
        for i = 1:m
            tmp = nzA(nzAind(i)+1:nzAind(i+1));
            normAi(i) = tmp'*tmp;
        end    
    end

    %% Minimization over Block X 
    function [x_n,X_n] = Xblockmin(z,zt,xi,mu)
        rho = opts.rho;
        tmp = accumarray(Hind,zt+mu/rho);
        if length(tmp) < n   %% in case all zero row in At
            tmp(length(tmp)+1:n) = 0;
        end
        x_n = invX.*(z+xi/rho-c/rho+tmp);   
        
        % block form
        X_n = blockify(X,x_n,K);      %% block form
    end

    %% Minimization over Block Y 
    function [z_n,zt_n,Omega,xtmp] = Yblockmin(X,x,Xi,xi,mu)
       %% projection in parrallel
       rho = opts.rho; 
       S   = cellfun(@(X,Xi)(X - Xi./rho),X,Xi,'UniformOutput',false);
       Z_n = projectK(S,K,0);
       z_n = flatten(zeros(size(x)),Z_n);  
       
       %% update zt in parrallel
       xtmp  = x(Hind);
       Omega = (-b + sbA*(xtmp-mu/rho))./normAi;
             
       tmp = repval(nzA,Omega,nzAind); %% using a mex function
%        tmp = zeros(size(nzA));
%        for i = 1:m
%            tmp(nzAind(i)+1:nzAind(i+1)) = nzA(nzAind(i)+1:nzAind(i+1))*Omega(i);
%        end     
       zt_n  = xtmp - mu/rho - tmp;     
    end

    %% update scaled multipliers, and calculate dual residua
    function [mu_n,xi_n,Xi_n] = Dupdate(mu,xi,Xi,z,zt,x,xtmp)
        mu_n = mu + opts.rho*(zt-xtmp);
        xi_n = xi + opts.rho*(z-x);
        Xi_n = blockify(Xi,xi_n,K);
    end

    function [isConverged,pres,dres] = ConverCheck(z_n,x_n,zt_n,z,zt,xtmp,mu,xi)
        persistent itPinf itDinf
        %Use the basic convergence test in the Boyd survey paper

        rho = opts.rho;
        %primal residual
        r     = (norm(z_n -x_n,'fro')^2 + norm(zt_n - xtmp,'fro')^2)^(1/2);
        pres  = r./max([norm(z_n,'fro'),norm(x_n,'fro'),norm(zt,'fro'),norm(xtmp,'fro')]);

        %dual residual
        s = rho.*(norm(z_n -z,'fro')^2 + norm(zt_n - zt,'fro')^2)^(1/2);
        dres  = s./((norm(mu,'fro')^2 + norm(xi,'fro')^2)^(1/2));

        %stopping criteria
        if(max(pres,dres)<opts.eps)
            isConverged = true;
        else
            isConverged = false;
        end

        % Update penalty parameter
%         if opts.adaptive
%             resRat = pres/dres;
%             if resRat > opts.mu
%                 opts.rho = opts.rho*opts.tau;
%             elseif 1/resRat > opts.mu
%                 opts.rho = opts.rho/opts.tau;
%             end
%         end
        if opts.adaptive
            resRat = pres/dres;
            if resRat > opts.mu
                itPinf = itPinf+1;
                itDinf = 0;
                if itPinf >= opts.rhoIt
                    % resRat remained small for long => rescale rho
                    itPinf = 0;
                    opts.rho = min(opts.rhoMax,opts.rho*opts.tau);
                end
            else
                itDinf = itDinf+1;
                itPinf = 0;
                if itDinf >= opts.rhoIt
                    % resRat remained small for long => rescale rho
                    itDinf = 0;
                    opts.rho = max(opts.rhoMin,opts.rho/opts.tau);
                end
            end
        end

    end


end


