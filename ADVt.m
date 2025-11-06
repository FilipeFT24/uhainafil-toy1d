function [g, obj] = ADVt(g, obj, penParam)
%--------------------------------------------------------------------------
ns   = size(g.tlvls, 2);
t0   = g.t;
dt   = 0;
tk   = g.data.tk;
nk   = numel(tk);
test = g.test;
%--------------------------------------------------------------------------
for j = 1:ns
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if test == 17
        %------------------------------------------------------------------
        K            = g.numE;
        N            = g.N;
        V            = 2;
        xv           = g.coordV(:, 1);
        xbl          = xv(  1, 1);
        xbr          = xv(K+1, 1);
        X            = zeros(K, N, V);
        X(:, :, 1)   = inittype(g.itype, @(x) g.data.Ur1(x, g.t), g.xydc, g.xyqc, g.fi_aux);
        X(:, :, 2)   = inittype(g.itype, @(x) g.data.Ur2(x, g.t), g.xydc, g.xyqc, g.fi_aux); % em p2 isto nao deve ser g.t
        %------------------------------------------------------------------
        La           = 20;
        Lr           = 10;
        xra          = (g.xydc-(g.xydc(end)-La))./La;
        xrg          = (g.xydc-(g.xydc(1)))./Lr;
        xra(xra < 0) = 0;
        xrg(xrg > 1) = 1;
        Cra          = sqrt(1-xra.^2);
        Crg          = (1-xrg).^1;

        sol1_old = g.x(:, :, 1);
        sol2_old = g.x(:, :, 2);

        for i = 1:2
            g.x(:, :, i) = (1-Crg).*g.x(:, :, i)+Crg.*X(:, :, i);
            switch i
                case 1
                    sol1_new = g.x(:, :, 1);
                case 2
                    sol2_new = g.x(:, :, 2);
                otherwise
                    return
            end
        end
        for i = 1:V
            g.x(:, :, i) = Cra.*g.x(:, :, i);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = t0+dt.*g.tlvls(1, j);
    g = FLUX(g, penParam, t);
    if j == 1
        x0 = g.x;
        dt = g.CFL*min(g.vollambda, [], 'all');
        for i = 1:nk
            if g.t < tk(1, i) && g.t+dt > tk(1, i)
                dt = tk(1, i)-g.t;
                break;
            end
        end
    end
    %----------------------------------------------------------------------
    x   = g.x;
    x   = g.omega(j, 1).*x0+g.omega(j, 2).*(x+dt.*g.res);
    g.x = x;
    if g.data.wetdry
        g = limN(g);
    end
    %----------------------------------------------------------------------
end
g.nit = g.nit+1;
g.t   = g.t+dt;
%--------------------------------------------------------------------------
%{
% IMPLICIT MANNING/FRICTION:
M            = manning(g);
g.x(:, :, 2) = g.x(:, :, 2)./(1+dt.*M);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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