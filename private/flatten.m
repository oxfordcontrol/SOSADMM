function z = flatten(z,Z)

%flatten a cell array in a long vector

%Zvec = cellfun(@(x)x(:),Z,'UniformOutput',false);
Zvec = cellfun(@(x)svec(x),Z,'UniformOutput',false);
z = vertcat(Zvec{:});

end
