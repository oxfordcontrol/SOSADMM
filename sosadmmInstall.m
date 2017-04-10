function sosadmmInstall

% SOSADMMINSTALL
%
% Install SOSADMM and compile required binaries.
%
% See also SOSADMM



% Compile some mex files from this package
here = pwd;
cd('private')
if (~isempty (strfind (computer, '64')))
    mexcmd = 'mex -largeArrayDims' ;
else
    mexcmd = 'mex' ;
end
eval([mexcmd, ' svec.c']);
eval([mexcmd, ' smat.c']);
eval([mexcmd, ' repval.c']);
cd(here)

% Finally add to path and save
addpath(here);
savepath

fprintf('\nCompilation completed successfully.\n');

end
