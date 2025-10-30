function [g] = limN(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS:
%--------------------------------------------------------------------------
K      = g.numE;
Kf     = K+1;
fl     = g.bfl;
fr     = g.bfr;
drytol = g.drytol;
g.fix  = false(K, 1);
g.zb   = g.zbinit;
%--------------------------------------------------------------------------
g      = lim0(g);
%--------------------------------------------------------------------------
Z_dof  = g.x(:, :, 1);
Hudof  = g.x(:, :, 2);
Zbdof  = g.zb;
H_dof  = Z_dof-Zbdof;
H_m    = meanval(g, H_dof);
Hum    = meanval(g, Hudof);
Z_m    = meanval(g, Z_dof);
Zbm    = meanval(g, Zbdof);
H_l    = H_dof*fl';
H_r    = H_dof*fr';
Z_l    = Z_dof*fl';
Z_r    = Z_dof*fr';
Zbl    = Zbdof*fl';
Zbr    = Zbdof*fr';


% if g.nit >= 6245 || any(abs(H_dof(12, :)) > eps, 'all')
%     xx = 1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through "problematic" faces:
%--------------------------------------------------------------------------
log2l           = H_l(1:K, 1) < drytol;
log2r           = H_r(1:K, 1) < drytol;
log2            = false(Kf, 1);
log2(1:Kf-1, 1) = log2(1:Kf-1, 1) | log2l;
log2(2:Kf  , 1) = log2(2:Kf  , 1) | log2r;
log3            = true(Kf, 1);
for i = 2:Kf-1
    if log2(i-1, 1) && log2(i+1, 1)
        log3(i, 1) = false;
    end
end
log2( 1, 1)     = false;
log2(Kf, 1)     = false;
log2            = log2 & log3;
f               = find(log2);
n               = size(f, 1);
%--------------------------------------------------------------------------


flag = 0;
if g.nit > 6145
    flag = 0;
    xx = 1;
end

if n > 1
    yy = 1;
end

%--------------------------------------------------------------------------
for i = 1:n
    %----------------------------------------------------------------------
    r    = f(i, 1);
    l    = r-1;
    L    = l-1;
    R    = r+1;
    Z_L  = Z_r(L, 1);
    Zrl  = Z_r(l, 1);
    Zlr  = Z_l(r, 1);
    Z_R  = Z_l(R, 1);
    dryl = H_l(l, 1) < drytol & H_r(l, 1) < drytol;
    dryr = H_l(r, 1) < drytol & H_r(r, 1) < drytol;
    %----------------------------------------------------------------------
    %
    if dryr
        if Z_l(r, 1) > Z_L%Zbl(r, 1) > Z_L-drytol
            g.x  (l, :, 1) = Z_L;
            g.x  (l, :, 2) = Hum(l, 1);
            g.zb (l, :)    =-H_m(l, 1)+Z_L;
            g.fix(l, 1)    = true;
        end
        if Zlr < Zrl
            g.x  (r, :, 1) = Zrl;
            g.x  (r, :, 2) = 0;
            g.zb (r, :)    = Zrl;
            g.fix(r, 1)    = true;
        end
        PLOT(0, g, l, r, Z_dof, Zbdof, drytol);
    end
    %
    if dryl
        if Z_r(l, 1) > Z_R%Zbr(l, 1) > Z_R-drytol
            g.x  (r, :, 1) = Z_R;
            g.x  (r, :, 2) = Hum(r, 1);
            g.zb (r, :)    =-H_m(r, 1)+Z_R;
            g.fix(r, 1)    = true;
        end
        if Zrl < Zlr
            g.x  (l, :, 1) = Zlr;
            g.x  (l, :, 2) = 0;
            g.zb (l, :)    = Zlr;
            g.fix(l, 1)    = true;
        end
        PLOT(0, g, l, r, Z_dof, Zbdof, drytol);
    end
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSITIVITY-PRESERVING LIMITER:
%--------------------------------------------------------------------------
if g.p > 0
    g = lim1(g);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end








function [] = PLOT(flag, g, l, r, Z_dof, Zbdof, drytol)
%
if flag
    switch g.p
        case 0
            aux = 1;
        case 1
            aux = [1, 2];
        case 2
            aux = [1, 3, 2];
        case 3
            aux = [1, 3, 4, 2];
        otherwise
            return
    end
    %
    f1 = figure;
    subplot(1, 2, 1);
    hold on;
    for j = l-1:r+1
        plot(g.xydc(j, aux), Z_dof(j, aux), '--ob');
        plot(g.xydc(j, aux), Zbdof(j, aux)+drytol,  ':*k');
    end
    subplot(1, 2, 2);
    hold on;
    for j = l-1:r+1
        plot(g.xydc(j, aux), g.x(j, aux, 1), '--ob');
        plot(g.xydc(j, aux), g.zb(j, aux)+drytol, ':*k');
    end
    close(f1);
end
end