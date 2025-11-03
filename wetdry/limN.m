function [g] = limN(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS:
%--------------------------------------------------------------------------
K     = g.numE;
Kf    = K+1;
fl    = g.bfl;
fr    = g.bfr;
tol   = g.drytol;
g.fix = false(K, 1);
g.zb  = g.zbinit;
%--------------------------------------------------------------------------
g     = lim0(g);
%--------------------------------------------------------------------------
Z_dof = g.x(:, :, 1);
Hudof = g.x(:, :, 2);
B_dof = g.zb;
H_dof = Z_dof-B_dof;
H_m   = meanval(g, H_dof);
Hum   = meanval(g, Hudof);
B_l   = B_dof*fl';
B_r   = B_dof*fr';
H_l   = H_dof*fl';
H_r   = H_dof*fr';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through "problematic" faces:
%--------------------------------------------------------------------------
log2l           = H_l(1:K, 1) < tol;
log2r           = H_r(1:K, 1) < tol;
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

H_aux = H_dof-tol;
if n > 2
    yy = 1;
end

%--------------------------------------------------------------------------
for i = 1:n
    %----------------------------------------------------------------------
    r    = f(i, 1);
    l    = r-1;
    L    = l-1;
    R    = r+1;
    Brl  = g.zb(l, :   )*fr';
    Blr  = g.zb(r, :   )*fl';
    Z_L  = g.x (L, :, 1)*fr';
    Zrl  = g.x (l, :, 1)*fr';
    Zlr  = g.x (r, :, 1)*fl';
    Z_R  = g.x (R, :, 1)*fl';
    dryl = log2l(l, 1) && log2r(l, 1);
    dryr = log2l(r, 1) && log2r(r, 1);
    %{
    Z_L = Z_r(L, 1);
    Zrl = Z_r(l, 1);
    Zlr = Z_l(r, 1);
    Z_R = Z_l(R, 1);
    %}

    % falta o tratamento de quando tens uma face entre wets

    %----------------------------------------------------------------------
    if ~dryl
        if ~dryr
            %--------------------------------------------------------------
            % wet/wet
            xx = 1;
            %--------------------------------------------------------------
        else
            %--------------------------------------------------------------
            % wet/dry
            if Z_L < Brl+tol % compare downstream water level to max. bathymetry level: if not enough to flood, set to Z_L
                Zrl            = Z_L;
                Brl            =-H_m(l, 1)+Zrl;
                g.x  (l, :, 1) = Zrl;
                g.x  (l, :, 2) = Hum(l, 1);
                g.zb (l, :)    = Brl;
                g.fix(l, 1)    = true;
            end
            %if Zrl < Zlr && B_r(l, 1) < B_l(r, 1)
                Zlr            = Zrl;%Brl+tol;%Zrl;%Brl+tol;
                g.x  (r, :, 1) = Zlr;
                g.x  (r, :, 2) = 0;
                g.zb (r, :)    = Zlr;
                g.fix(r, 1)    = true;
                
                PLOT(0, g, l, r, Z_dof, B_dof, tol);
            %end

        

            %--------------------------------------------------------------
        end
    else
        %------------------------------------------------------------------
        % dry/wet
        if Z_R < Brl+tol
            Zlr            = Z_R;
            Blr            =-H_m(r, 1)+Zlr; 
            g.x  (r, :, 1) = Zlr;
            g.x  (r, :, 2) = Hum(r, 1);
            g.zb (r, :)    = Blr;
            g.fix(r, 1)    = true;
        end
        Zrl            = Zlr; % Z_L
        g.x  (l, :, 1) = Zrl;
        g.x  (l, :, 2) = 0;
        g.zb (l, :)    = Zrl;
        g.fix(l, 1)    = true;
        %------------------------------------------------------------------
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


function [] = PLOT(flag, g, l, r, Z_dof, B_dof, tol)
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
        plot(g.xydc(j, aux), B_dof(j, aux),  ':*k');
    end
    subplot(1, 2, 2);
    hold on;
    for j = l-1:r+1
        plot(g.xydc(j, aux), g.x(j, aux, 1), '--ob');
        plot(g.xydc(j, aux), g.zb(j, aux), ':*k');
    end
    close(f1);
end
end