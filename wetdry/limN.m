function [g] = limN(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS #1:
%--------------------------------------------------------------------------
K      = g.numE;
Kf     = K+1;
fc     = g.bf;
fl     = g.bfl;
fr     = g.bfr;
rtol   = eps;
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
log1 = H_m(:, 1) < 0 | (H_l(:, 1) < drytol & H_r(:, 1) < drytol);
c    = find(log1);
m    = size(c, 1);
%--------------------------------------------------------------------------
for i = 1:m
    o              = c    (i, 1);
    Zbo            = Zbdof(o, :);
    H_dof(o, :)    = 0;
    Hudof(o, :)    = 0;
    Z_dof(o, :)    = Zbo;
    H_m  (o, 1)    = 0;
    Hum  (o, 1)    = 0;
    g.x  (o, :, 1) = Zbo;
    g.x  (o, :, 2) = 0;
    H_l  (o, 1)    = 0;
    H_r  (o, 1)    = 0;
    Z_l  (o, 1)    = Zbl(o, 1);
    Z_r  (o, 1)    = Zbr(o, 1);
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

H1    = g.x(:, :, 1)-g.zbinit-drytol;
Haux1 = [H1*fl', H1*fr'];
flag = 0;
if n > 2
    stop = 1;
    flag = 0;
end

if g.nit > 6500
    flag = 0;
    x = 1;
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
    H_L = H_r(L, 1);
    H_R = H_l(r, 1);
    %----------------------------------------------------------------------
    if log2l(r, 1) && (~log2r(l, 1) || log2r(l, 1) && H_R > H_L)
        %------------------------------------------------------------------
        % WET/DRY
        %------------------------------------------------------------------
        Z_X            = min(Z_L, Z_R);
        g.x  (l, :, 1) = Z_X;
        g.x  (l, :, 2) = Hum(l, 1);
        g.zb (l, :)    =-H_m(l, 1)+Z_X;
        g.fix(l, 1)    = true;
        if (Z_L < Z_R && Z_l(r, 1) < Z_X) || (Z_L > Z_R && Z_l(r, 1) < Z_X) % up/down %%% ta aqui qq coisa mal mas as condiçoes parecem ser estas
            g.x  (r, :, 1) = Z_X;
            g.x  (r, :, 2) = Hum(r, 1);
            g.zb (r, :)    =-H_m(r, 1)+Z_X;
            g.fix(r, 1)    = true;

            flag = 0;
        end
        %------------------------------------------------------------------
        %flag = 0;
        PLOT(flag, g, l, r, Z_dof, Zbdof, drytol);
    else
        %------------------------------------------------------------------
        % DRY/WET
        %------------------------------------------------------------------
        Z_X            = min(Z_L, Z_R);
        g.x  (r, :, 1) = Z_X;
        g.x  (r, :, 2) = Hum(r, 1);
        g.zb (r, :)    =-H_m(r, 1)+Z_X;
        g.fix(r, 1)    = true;
        if (Z_R < Z_L && Z_r(l, 1) < Z_X) || (Z_R > Z_L && Z_r(l, 1) < Z_X) % up/down
            g.x  (l, :, 1) = Z_X;
            g.x  (l, :, 2) = Hum(l, 1);
            g.zb (l, :)    =-H_m(l, 1)+Z_X;
            g.fix(l, 1)    = true;

            flag = 0;
        end
        %------------------------------------------------------------------
        %flag = 0;
        PLOT(flag, g, l, r, Z_dof, Zbdof, drytol);
    end
    %----------------------------------------------------------------------
end

H2    = g.x(:, :, 1)-g.zbinit-drytol;
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
        plot(g.xydc(j, aux), g.x (j, aux, 1), '--ob');
        plot(g.xydc(j, aux), g.zb(j, aux)+drytol, ':*k');
    end
    close(f1);
end
end