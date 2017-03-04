function [varargout] = makeConeVariables(K)

nCones = length(K.s);
varargout = repmat({{}},[1 nargout]);

shift = 0;
%populate the free variables
if(isfield(K,'f') && K.f > 0)
    for j = 1:nargout
        varargout{j}{1+shift} = zeros(K.f,1);
    end
    shift = 1;
end

%populate the zero variables
if(isfield(K,'z') && K.z > 0)
    for j = 1:nargout
        varargout{j}{1+shift} = zeros(K.z,1);
    end
    shift = shift+1;
end

%populate the linear cone variables
if(isfield(K,'l') && K.l > 0)
    for j = 1:nargout
        varargout{j}{1+shift} = zeros(K.l,1);
    end
    shift = shift+1;
end

if isfield(K,'q') && sum(K.q) > 0
    error(sprintf('Second-order cone not yet implemented in %s!',mfilename));
end

%populate the cones
for i = (1:nCones)
    for j = 1:nargout
        varargout{j}{end+1} = zeros(K.s(i),K.s(i));
    end
end
