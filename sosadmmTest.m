function sosadmmTest

% SOSADMMTEST
%
% Run some examples to test sosadmm.
%

% Preliminaries
clc;
opts.maxIter = 1e+3;
opts.relTol  = 1e-3;
here = pwd;
cd examples

% ---------------------------------------------------------------------------- %
%           SDP with row sparsity
% ---------------------------------------------------------------------------- %
% Parameters
m = 500;                     % # constraints
n = 200;                     % 
density = 1.e-3;

% Setup
fprintf('\nSetting up row sparse SDP, m=%i...',m);
tsetup = tic;
[At,b,c,K] = RowSpaSDP(m,n,density);
tsetup = toc(tsetup);
fprintf('done in %.2f seconds. \n',tsetup);

% solution by admm
sosadmm(At,b,c,K,opts);

%%
% Parameters
m = 1000;                    % # constraints
n = 300;                     % 
density = 1.e-3;

% Setup
fprintf('\nSetting up row sparse SDP, m=%i...',m);
tsetup = tic;
[At,b,c,K] = RowSpaSDP(m,n,density);
tsetup = toc(tsetup);
fprintf('done in %.2f seconds. \n',tsetup);

% solution by admm
sosadmm(At,b,c,K,opts);


% ---------------------------------------------------------------------------- %
%                              SDPs in arising in SOS programs
% ---------------------------------------------------------------------------- %
% exSOS
fprintf('Testing an SOS feasibility problem\n');
pause(1)
load(['examples',filesep,'exSOS.mat'])
sosadmm(At,b,c,K,opts);
 
 
% exLyapunov
fprintf('\nTesting an example of Finding Lyapunov functions \n');
pause(1)
load(['examples',filesep,'exLyapunov.mat'])
sosadmm(At,b,c,K,opts);


% ---------------------------------------------------------------------------- %
%                                   END
% ---------------------------------------------------------------------------- %
cd(here);
fprintf('\n\nSOSADMM was successfully tested.\n\n')
end



% ============================================================================ %
%                               NESTED FUNCTIONS                               %
% ============================================================================ %

% -------------------
% SDP with row sparsity
% -------------------
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


end


% ============================================================================ %
%                        END OF NESTED FUNCTIONS                               %
% ============================================================================ %
