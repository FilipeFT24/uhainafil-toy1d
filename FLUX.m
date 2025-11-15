function [g] = FLUX(g, penParam, t)
%--------------------------------------------------------------------------
g     = SWEc1(g, t, penParam);
g     = SWEf1(g, t);
g.res = g.Fluxc+g.Fluxf;
g.res =-g.res;
%--------------------------------------------------------------------------
% P0 correction:
%{
if g.p > 0
    res_aux = meanval(g, g.res);
    x       = g.fix;
    for v = 1:2
        g.res(x, :, v) = repmat(res_aux(x, v), 1, g.N);
    end
end
%}
%--------------------------------------------------------------------------
end