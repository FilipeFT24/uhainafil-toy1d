function [g, obj] = ADVt(g, obj, penParam)
%--------------------------------------------------------------------------
ns = size(g.tlvls, 2);
t0 = g.t;
dt = 0;
tk = g.data.tk;
nk = numel(tk);
%--------------------------------------------------------------------------
for j = 1:ns
    %----------------------------------------------------------------------
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
    x    = g.x;
    x    = g.omega(j, 1).*x0+g.omega(j, 2).*(x+dt.*g.res);
    g.x  = x;

    Z_dof            = g.x(:, :, 1);
    Hudof            = g.x(:, :, 2);
    Zbdof            = g.zb;
    H_dof            = Z_dof-Zbdof;
    H_m              = meanval(g, H_dof);


    xx = 1;

    if g.nit == 8242
        stop = 1;
    end

    
    if g.data.wetdry
        g = limX(g);
    end


    Z_dof1            = g.x(:, :, 1);
    Hudof1            = g.x(:, :, 2);
    Zbdof1            = g.zb;
    H_dof1            = Z_dof1-Zbdof1;
    H_m1               = meanval(g, H_dof1);

    xx = 1;



    %----------------------------------------------------------------------
end
%--------------------------------------------------------------------------
% IMPLICIT MANNING/FRICTION:
% M            = manning(g);
% g.x(:, :, 2) = g.x(:, :, 2)./(1+dt.*M);
%--------------------------------------------------------------------------
g.nit = g.nit+1;
g.t   = g.t+dt;
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
end