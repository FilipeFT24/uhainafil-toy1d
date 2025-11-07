function [g] = mymatlabFunction(f, vars)
if ~has(f, vars(1, 1))
    c = matlabFunction(f, 'Vars', vars);
    g = @(x, varargin) c(varargin{:}).*ones(size(x));
else
    g = matlabFunction(f, 'Vars', vars);
end
end