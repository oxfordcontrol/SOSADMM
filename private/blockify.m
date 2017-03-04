function X = blockify(X,x,K)

%Turn portions of x into the matrices X_i 
%that must satisfy X_i \succeq 0

%NB: X is both the input and output, since this
%enables a sort-of call by reference in matlab,
%and no copies are made (seemingly)

if(isfield(K,'f') && K.f > 0)
    X{1}(:) = x(1:K.f);
    shift = 1; 
    count = K.f;
else
    shift = 0; 
    count = 0; 
end

if(isfield(K,'z') && K.z > 0)
    X{1+shift}(:) = x(count+(1:K.z));
    shift = shift+1; 
    count = count + K.z;
end

if(isfield(K,'l') && K.l > 0)
    X{1+shift}(:) = x(count+(1:K.l));
    shift = shift+1; 
    count = count + K.l;
end

if isfield(K,'q') && sum(K.q) > 0
    error('Second-order cone not yet implemented in %s!',mfilename);
end


for ii = 1:length(K.s)
    d = K.s(ii);    %dimension of this SDP cone
    if d>0
%         X{ii+shift}(:) = x((count+1):(count+d^2),1);
%         count = count + d^2;
        coneDim = d*(d+1)/2;
        X{ii+shift} = smat(x((count+1):(count+coneDim),1));
        count = count + coneDim;
    end
end

end

