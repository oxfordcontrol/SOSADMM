
%% Non-orthongnality of the constraints in an SOS program
% You need to install SOSTOOLS, to run this example
clc;clear

pvar a b x

prog = sosprogram(x);       % initialization

prog = sosdecvar(prog,a);   
prog = sosdecvar(prog,b);

p1 = a*x^4 + b*x^2 + x + 1;
p2 = b*x^4 + a*x^2 + x + 1;
prog = sosineq(prog, p1);
prog = sosineq(prog, p2);
prog = sossolve(prog);
