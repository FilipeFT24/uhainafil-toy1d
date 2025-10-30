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
rtol   = eps;
g.fix  = false(K, 1);
g.zb   = g.zbinit;
%--------------------------------------------------------------------------
Z_dof  = g.x(:, :, 1);
Hudof  = g.x(:, :, 2);
Zbdof  = g.zbinit;
H_dof  = Z_dof-Zbdof;
H_m    = meanval(g, H_dof);
Hum    = meanval(g, Hudof);
Z_m    = meanval(g, Z_dof);
Zbm    = meanval(g, Zbdof);
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
log1 = (H_l(:, 1) < 0 & H_r(:, 1) < 0); %H_m(:, 1) < 0 | 
c    = find(log1);
m    = size(c, 1);
%--------------------------------------------------------------------------
for i = 1:m
    %----------------------------------------------------------------------
    o              = c    (i, 1);
    Zbo            = Zbdof(o, :);
    H_o            = 0;
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

H0    = g.x(:, :, 1)-g.zb;
H1    = g.x(:, :, 1)-g.zb-drytol;
Haux0 = [H0*fl', H0*fr'];
Haux1 = [H1*fl', H1*fr'];

flag = 0;
if n > 2
    xx = 1;
end

%--------------------------------------------------------------------------
for i = 1:min(2, n)
    %----------------------------------------------------------------------
    r    = f(i, 1);
    l    = r-1;
    L    = l-1;
    R    = r+1;
    Z_L  = Z_r(L, 1);
    Z_R  = Z_l(R, 1);
    dryl = H_m(l, 1) < drytol;
    dryr = H_m(r, 1) < drytol;
    %----------------------------------------------------------------------
%     if dryr
%         if Zbl(r, 1) > Z_L%-drytol
%             g.x  (l, :, 1) = Z_L;
%             g.x  (l, :, 2) = Hum(l, 1);
%             g.zb (l, :)    =-H_m(l, 1)+Z_L;
%             g.fix(l, 1)    = true;
%             if Z_l(r, 1) < Z_L
%                 g.x  (r, :, 1) = Z_L;
%                 g.x  (r, :, 2) = Hum(r, 1);
%                 g.zb (r, :)    =-H_m(r, 1)+Z_L;
%                 g.fix(r, 1)    = true;
%             end
%         end
%         %PLOT(1, g, l, r, Z_dof, Zbdof, drytol);
%     end
    if dryl
        if Z_r(l, 1) > Z_R%-drytol
            g.x  (r, :, 1) = Z_R;
            g.x  (r, :, 2) = Hum(r, 1);
            g.zb (r, :)    =-H_m(r, 1)+Z_R;
            g.fix(r, 1)    = true;
            if Z_r(l, 1) < Z_R
                g.x  (l, :, 1) = Z_R;
                g.x  (l, :, 2) = Hum(l, 1);
                g.zb (l, :)    =-H_m(l, 1)+Z_R;
                g.fix(l, 1)    = true;
            end
            %PLOT(1, g, l, r, Z_dof, Zbdof, drytol);
            %disp(1);
        end
    end
end




% %--------------------------------------------------------------------------
% for i = 1:min(2, n)
%     %----------------------------------------------------------------------
%     r      = f(i, 1);
%     l      = r-1;
%     L      = l-1;
%     R      = r+1;
%     Z_L    = Z_r(L, 1);
%     Z_R    = Z_l(R, 1);
%     Z_X    = min(Z_L, Z_R);
%     %----------------------------------------------------------------------
%     if Z_r(l, 1) > Z_l(r, 1)
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %------------------------------------------------------------------
% 
%         flag11         = Z_L < Z_R && Z_l(r, 1) < Z_L; % run-up   (l->r)
%         flag12         = Z_L > Z_R && Z_r(l, 1) > Z_R; % run-down (l->r)
%         flag13         = Z_L < Z_R && Z_l(r, 1) > Z_L; % trans. 1 (l->r) : modified run-up   (l->r)
%         flag14         = Z_L > Z_R && Z_r(l, 1) < Z_R; % trans. 2 (l->r) : modified run-down (l->r)
%         %
%         %if flag11 || flag12
%             g.x  (l, :, 1) = Z_X;
%             g.x  (l, :, 2) = Hum(l, 1);
%             g.zb (l, :)    =-H_m(l, 1)+Z_X;
%             g.fix(l, 1)    = true;
% 
%             if flag11 || flag12
%                 g.x  (r, :, 1) = Z_X;
%                 g.x  (r, :, 2) = Hum(r, 1);
%                 g.zb (r, :)    =-H_m(r, 1)+Z_X;
%                 g.fix(r, 1)    = true;
%             end
%         %end
%         PLOT(0, g, l, r, Z_dof, Zbdof, drytol);
%         %------------------------------------------------------------------
%     else
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %------------------------------------------------------------------
% 
%         flag21         = Z_L > Z_R && Z_r(l, 1) < Z_R; % run-up   (r->l)
%         flag22         = Z_L < Z_R && Z_l(r, 1) > Z_L; % run-down (r->l)
%         flag23         = Z_L > Z_R && Z_r(l, 1) > Z_R; % trans. 1 (r->l) : modified run-up   (r->l)
%         flag24         = Z_L < Z_R && Z_l(r, 1) < Z_L; % trans. 2 (r->l) : modified run-down (r->l)
%         %
%         %if flag21 || flag22 || flag23 || flag24
%             g.x  (r, :, 1) = Z_X;
%             g.x  (r, :, 2) = Hum(r, 1);
%             g.zb (r, :)    =-H_m(r, 1)+Z_X;
%             g.fix(r, 1)    = true;
%             if flag21 || flag22
%             g.x  (l, :, 1) = Z_X;
%             g.x  (l, :, 2) = Hum(l, 1);
%             g.zb (l, :)    =-H_m(l, 1)+Z_X;
%             g.fix(l, 1)    = true;
%             end
%         %end
%         PLOT(0, g, l, r, Z_dof, Zbdof, drytol);
%         %------------------------------------------------------------------
%     end
%     %----------------------------------------------------------------------
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSITIVITY-PRESERVING LIMITER
%--------------------------------------------------------------------------
% Z_dof            = g.x(:, :, 1);
% Hudof            = g.x(:, :, 2);
% Zbdof            = g.zb;
% H_dof            = Z_dof-Zbdof;
% H_m              = meanval(g, H_dof);
% log_w            =~g.fix & H_m > drytol;
% H_dofw           = H_dof(log_w, :);
% Hudofw           = Hudof(log_w, :);
% Zbdofw           = Zbdof(log_w, :);
% H_mw             = H_m  (log_w, 1);
% Humw             = Hum  (log_w, 1);
% %--------------------------------------------------------------------------
% H_quad           = H_dofw*fc';
% min_H            = min(H_quad, [], 2);
% theta            = min(1, (H_mw-drytol)./(H_mw-min_H)); % uhaina uses 0 instead of drytol
% theta            = max(0, theta);
% %--------------------------------------------------------------------------
% g.x(log_w, :, 1) = H_mw+theta.*(H_dofw-H_mw)+Zbdofw;
% g.x(log_w, :, 2) = Humw+theta.*(Hudofw-Humw);
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