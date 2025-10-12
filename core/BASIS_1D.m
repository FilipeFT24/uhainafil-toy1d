function [d, f, xp] = BASIS_1D(P)
N = P+1;
d = cell(1, N);
f = cell(1, N);
switch P
    case 0
        xp = [1./2]; %#ok<NBRAK2> 
    case 1
        xp = [0; 1];
    case 2
        xp = [0; 1; 1./2];
    case 3
        xp = [0; 1; 1./3; 2./3];
    case 4
        xp = [0; 1; 1./4; 1./2; 3./4];
    otherwise
        return;
end
%{
if P == 0
    xp = 0.5;
else
    xp = linspace(0, 1, P+1)';
end
%}
syms x;
for i = 1:N
    fi = sym(1);
    for j = 1:N
        if j ~= i
            fi = fi.*(x-xp(j, 1))./(xp(i, 1)-xp(j, 1));
        end
    end
    di = simplify(diff(fi, x));
    if di == 0
        d{1, i} = @(x) repmat(0         , size(x)); %#ok<REPMAT> 
        f{1, i} = @(x) repmat(double(fi), size(x));
    else
        if simplify(diff(di, x)) == 0
            d{1, i} = @(x) repmat(double(di), size(x));
        else
            d{1, i} = matlabFunction(di, 'Vars', x);
        end
        f{1, i} = matlabFunction(simplify(fi), 'Vars', x);
    end
end
end