function [g, obj] = ADVt(g, obj, penParam)
%--------------------------------------------------------------------------
K    = g.numE;
N    = g.N;
D    = 1;
V    = 1+D;
ns   = size(g.rkc, 1);
t0   = g.t;
xold = g.x;
xtmp = zeros(K, N, V);
ks   = zeros(K, N, V, ns);
%--------------------------------------------------------------------------
tk   = g.data.tk;
nk   = size(tk, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g              = FLUX(g, penParam, t0);
ks(:, :, :, 1) = g.res;
dt             = g.CFL*min(g.vollambda, [], 'all');
for i = 1:nk
    if t0 < tk(1, i) && t0+dt > tk(1, i)
        dt = tk(1, i)-t0;
        break;
    end
end
xtmp = xtmp+g.rkb(1, 1).*ks(:, :, :, 1);
%--------------------------------------------------------------------------
if ns > 1
    for i = 2:ns
        ts  = t0+g.rkc(i, 1).*dt;
        xs  = xold;
        for j = 1:i-1
            xs = xs+dt.*g.rka(i, j).*ks(:, :, :, j);
        end
        g.x = xs;
        %{
        if g.data.wetdry
            g = limN(g);
        end
        %}
        g              = FLUX(g, penParam, ts);
        ks(:, :, :, i) = g.res;
        xtmp           = xtmp+g.rkb(1, i).*ks(:, :, :, i);
    end
end
%--------------------------------------------------------------------------
g.x         = xold+dt.*xtmp;
%{
if g.data.wetdry
    g = limN(g);
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% IMPLICIT MANNING/FRICTION:
M            = manning(g);
g.x(:, :, 2) = g.x(:, :, 2)./(1+dt.*M);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g.t   = g.t+dt;
g.nit = g.nit+1;
resn2 = zeros(1, 2);
%     = squeeze(max(max(-g.res, [], 1), [], 2))';
for i = 1:2
    res_quad    = g.res(:, :, i)*g.bf';
    resc        = res_quad.^2*g.W'.*g.detJ0T(:, 1);
    resn2(1, i) = sqrt(sum(resc, 1));
end
fprintf("%5d | %.3e %.3e | %.3e, %.3e\n", g.nit, g.t-dt, dt, resn2(1, 1), resn2(1, 2));
%--------------------------------------------------------------------------
obj = PlotB(g, obj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end