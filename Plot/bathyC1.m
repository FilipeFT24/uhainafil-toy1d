function [f] = bathyC1(slope, ds, flag)
d1 =-ds;
d2 = ds;
A  = [...
    d1.^[3, 2, 1, 0].*[1, 1, 1, 1];
    d1.^[2, 1, 0, 0].*[3, 2, 1, 0];
    d2.^[3, 2, 1, 0].*[1, 1, 1, 1];
    d2.^[2, 1, 0, 0].*[3, 2, 1, 0]];
switch flag
    case 1
        b = [0; 0; slope.*ds; +slope]; % _/
    case 2
        b = [slope.*ds; -slope; 0; 0]; % \_
    otherwise
        return
end
x  = A\b;
f  = @(y) ...
    x(1, 1).*y.^3+...
    x(2, 1).*y.^2+...
    x(3, 1).*y.^1+...
    x(4, 1);
end