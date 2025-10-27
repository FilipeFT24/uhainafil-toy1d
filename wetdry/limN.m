function [g] = limN(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS #1:
%--------------------------------------------------------------------------
K      = g.numE;
Kf     = K+1;
fc     = g.bf;
fl     = g.bfl;
fr     = g.bfr;
drytol = g.drytol;
g.fix  = false(K, 1);
g.zb   = g.zbinit;
%--------------------------------------------------------------------------
Z_dof  = g.x(:, :, 1);
Hudof  = g.x(:, :, 2);
Zbdof  = g.zbinit;
H_dof  = Z_dof-Zbdof;
H_m    = meanval(g, H_dof);
Hum    = meanval(g, Hudof);
Z_l    = Z_dof*fl';
Z_r    = Z_dof*fr';
Zbl    = Zbdof*fl';
Zbr    = Zbdof*fr';
H_l    = Z_l-Zbl;
H_r    = Z_r-Zbr;

% Zbmax  = max(Zbl(2:K, 1), Zbr(1:K-1, 1));
% H_l    = Z_l-[Zbl(1, 1); Zbmax];
% H_r    = Z_r-[Zbmax; Zbr(K, 1)];%max(Zbl, Zbr);%
% gg = H_l(2:K)-H_r(1:K-1);
% xx = 1;
%     if g.nit == 0
%     log2l(41, 1) = true;
%     log2r(40, 1) = true;
%     end
%% have to do an initial limiting for projection!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through "problematic" cells:
%--------------------------------------------------------------------------
log1 = H_m(:, 1) < drytol | (H_l(:, 1) < drytol & H_r(:, 1) < drytol);
c    = find(log1);
m    = size(c, 1);
%--------------------------------------------------------------------------
for i = 1:m
    %----------------------------------------------------------------------
    o              = c    (i, 1);
    Zbo            = Zbdof(o, :);
    H_o            = drytol-eps.*10;
    %----------------------------------------------------------------------
    H_dof(o, :)    = H_o;
    Hudof(o, :)    = 0;
    Z_dof(o, :)    = H_o+Zbo;
    %----------------------------------------------------------------------
    H_m  (o, 1)    = H_o;
    Hum  (o, 1)    = 0;
    H_l  (o, 1)    = H_o;
    H_r  (o, 1)    = H_o;
    Z_l  (o, 1)    = Z_dof(o, :)*fl';
    Z_r  (o, 1)    = Z_dof(o, :)*fr';
    Zbl  (o, 1)    = Z_l  (o, 1);
    Zbr  (o, 1)    = Z_r  (o, 1);
    %----------------------------------------------------------------------
    g.x  (o, :, 1) = Z_dof(o, :);
    g.x  (o, :, 2) = Hudof(o, :);
    g.zb (o, :)    = Z_dof(o, :);
    %----------------------------------------------------------------------
end
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

H1    = g.x(:, :, 1)-g.zb-drytol;
Haux1 = [H1*fl', H1*fr'];

flag = 0;
if n > 2 || g.nit > 6753
    flag = 0;
end

if n > 2
    xx = 1;
end

%--------------------------------------------------------------------------
for i = 1:n
    %----------------------------------------------------------------------
    r   = f(i, 1);
    l   = r-1;
    L   = l-1;
    R   = r+1;
    Z_L = Z_r(L, 1);
    Z_R = Z_l(R, 1);
    %----------------------------------------------------------------------
    if Z_L <= Z_R
        %------------------------------------------------------------------
        g.x  (l, :, 1) = Z_L;
        g.x  (l, :, 2) = Hum(l, 1);
        g.zb (l, :)    =-H_m(l, 1)+Z_L;
        g.fix(l, 1)    = true;
        if (Z_r(l, 1) > Z_l(r, 1) && Z_l(r, 1) < Z_L) || (Z_r(l, 1) < Z_l(r, 1) && Z_l(r, 1) > Z_L) % run-up/down
            g.x  (r, :, 1) = Z_L;
            g.x  (r, :, 2) = Hum(r, 1);
            g.zb (r, :)    =-H_m(r, 1)+Z_L;
            g.fix(r, 1)    = true;
        end
        if (Z_r(l, 1) < Z_l(r, 1) && Z_l(r, 1) > Z_L)
            xx = 1;
        end


        PLOT(0, g, l, r, Z_dof, Zbdof, drytol);
        %------------------------------------------------------------------
    else
        %------------------------------------------------------------------
        g.x  (r, :, 1) = Z_R;
        g.x  (r, :, 2) = Hum(r, 1);
        g.zb (r, :)    =-H_m(r, 1)+Z_R;
        g.fix(r, 1)    = true;
        if (Z_l(r, 1) > Z_r(l, 1) && Z_r(l, 1) < Z_R) || (Z_l(r, 1) < Z_r(l, 1) && Z_r(l, 1) > Z_R) % run-up/down
            g.x  (l, :, 1) = Z_R;
            g.x  (l, :, 2) = Hum(l, 1);
            g.zb (l, :)    =-H_m(l, 1)+Z_R;
            g.fix(l, 1)    = true;
        end
        PLOT(flag, g, l, r, Z_dof, Zbdof, drytol);
        


        %------------------------------------------------------------------
    end
    %----------------------------------------------------------------------
end

H2    = g.x(:, :, 1)-g.zb-drytol;
Haux2 = [H2*fl', H2*fr'];
a1f   = find(Haux1(:, 2) < 0, 1, 'first');
a1l   = find(Haux1(:, 1) < 0, 1, 'last');

if ~isempty(a1f)
    if Haux2(a1f-1, 2) < 0% && Haux1(a1l, 1) < 0
        xx = 1;
    end
end
if ~isempty(a1l)
    if Haux2(a1l+1, 1) < Haux2(a1l, 2)
        xx = 1;
    end
end

if Haux2(323, 1) < Haux2(322, 2)
    xx = 1;
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
        plot(g.xydc(j, aux), Zbdof(j, aux),  ':*k');
    end
    subplot(1, 2, 2);
    hold on;
    for j = l-1:r+1
        plot(g.xydc(j, aux), g.x  (j, aux, 1), '--ob');
        plot(g.xydc(j, aux), Zbdof(j, aux), ':*k');
    end
    close(f1);
end
end