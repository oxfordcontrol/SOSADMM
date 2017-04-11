function opts = setSOSADMMopts

% SETSOSADMMOPTS
%
% Set defalt options for sosadmm
%
% See also SOSDAMM

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

end