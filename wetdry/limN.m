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
Z_m    = meanval(g, Z_dof);
Zbm    = meanval(g, Zbdof);
Z_l    = Z_dof*fl';
Z_r    = Z_dof*fr';
Zbl    = Zbdof*fl';
Zbr    = Zbdof*fr';
H_l    = Z_l-Zbl;
H_r    = Z_r-Zbr;

%Zbmax  = max(Zbl(2:K, 1), Zbr(1:K-1, 1));
% H_l    = Z_l-[Zbl(1, 1); Zbmax];
% H_r    = Z_r-[Zbmax; Zbr(K, 1)];%max(Zbl, Zbr);%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through "problematic" cells:
%--------------------------------------------------------------------------
% log1 = H_l(:, 1) <= drytol & H_r(:, 1) <= drytol & H_m(:, 1) <= drytol;
% c    = find(log1);
% m    = size(c, 1);
% %--------------------------------------------------------------------------
% for i = 1:m
%     %----------------------------------------------------------------------
%     % VARS:
%     o              = c  (i, 1);
%     Zbo            = Zbm(o, 1);
%     %----------------------------------------------------------------------
%     % VARS MODIFIED:
%     H_o            = drytol;
%     Z_o            = H_o+Zbo;
%     Huo            = Hum(o, 1).*(H_m(o, 1)./H_o);
%     %----------------------------------------------------------------------
%     % SOLUTION:
%     H_dof(o, :)    = H_o;
%     Hudof(o, :)    = Huo;
%     Z_dof(o, :)    = Z_o;
%     Zbdof(o, :)    = Zbo;
%     H_m  (o, 1)    = H_o;
%     Hum  (o, 1)    = Huo;
%     Z_m  (o, 1)    = Z_o;
%     Zbm  (o, 1)    = Zbo;
%     %
%     g.x  (o, :, 1) = Z_o;
%     g.x  (o, :, 2) = Huo;
%     g.zb (o, :)    = Zbo;
%     H_l  (o, 1)    = H_o;
%     H_r  (o, 1)    = H_o;
%     Z_l  (o, 1)    = Z_o;
%     Z_r  (o, 1)    = Z_o;
%     Zbl  (o, 1)    = Zbo;
%     Zbr  (o, 1)    = Zbo;
%     g.fix(o, 1)    = true;
%     %----------------------------------------------------------------------
% end

z_dof = g.x(:, :, 1);
zbdof = g.zb;
h_dof = Z_dof-zbdof;

if any(h_dof < 0, 'all')
    xx = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through "problematic" faces:
%--------------------------------------------------------------------------
n = 0;
if g.p > 0
    %----------------------------------------------------------------------
    dry             = H_l(1:K, 1) < drytol & H_r(1:K, 1) < drytol;
    log2l           = H_l(1:K, 1) < drytol;
    log2r           = H_r(1:K, 1) < drytol;
    log2            = false(Kf, 1);
    log2(1:Kf-1, 1) = log2(1:Kf-1, 1) | log2l;
    log2(2:Kf  , 1) = log2(2:Kf  , 1) | log2r;
    for i = 2:Kf-1
        if dry(i-1, 1) && dry(i, 1)
            log2(i, 1) = false;
        end
    end
    log2( 1, 1)     = false;
    log2(Kf, 1)     = false;
    %----------------------------------------------------------------------
    f               = find(log2);
    n               = size(f, 1);
    %----------------------------------------------------------------------
end

H_x = [H_l, H_r];
H_aux = [H_l, H_r]-drytol;
switch g.p
    case 0
        aux = 1;
    case 1
        aux = [1, 2];
    case 2
        aux = [1, 3, 2];
    case 3
        aux = [1, 3, 4, 1];
    otherwise
        return
end
plot_ = 1;

%--------------------------------------------------------------------------
for i = 1:n
    %----------------------------------------------------------------------
    r = f(i, 1);
    l = r-1;
    L = l-1;
    R = r+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if log2l(r, 1)
        if ~log2r(l, 1) || (log2r(l, 1) && Zbr(l, 1) > Z_r(L, 1))
            %--------------------------------------------------------------
            % WET/DRY
            % DRY/DRY w/ Zbr(l, 1) > Z_r(L, 1) => negative water depths
            %--------------------------------------------------------------
%             if Z_m(l, 1) > Z_m(r, 1)
%                 Z_1 = Z_m(l, 1);
%             else
                Z_1 = Z_r(L, 1);
                if Z_l(r, 1) < Z_1
                    Hu2            = Hum(r, 1);
                    Zb2            =-H_m(r, 1)+Z_1;
                    g.x  (r, :, 1) = Z_1;
                    g.x  (r, :, 2) = Hu2;
                    g.zb (r, :)    = Zb2;
                    g.fix(r, 1)    = true;
                end
%             end
            Hu1            = Hum(l, 1);
            Zb1            =-H_m(l, 1)+Z_1;
            g.x  (l, :, 1) = Z_1;
            g.x  (l, :, 2) = Hu1;
            g.zb (l, :)    = Zb1;
            g.fix(l, 1)    = true;
            %--------------------------------------------------------------

            if plot_
                figure;
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
                close all;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if log2r(l, 1)
        if ~log2l(r, 1) || (log2l(r, 1) && Zbl(r, 1) > Z_l(R, 1))
            %--------------------------------------------------------------
            % DRY/WET
            % DRY/DRY w/ Zbl(r, 1) > Z_l(R, 1) => negative water depths
            %--------------------------------------------------------------
%             if Z_m(r, 1) > Z_m(l, 1)
%                 Z_1 = Z_m(r, 1);
%             else
                Z_1 = Z_l(R, 1);
                if Z_r(l, 1) < Z_1
                    Hu2            = Hum(l, 1);
                    Zb2            =-H_m(l, 1)+Z_1;
                    g.x  (l, :, 1) = Z_1;
                    g.x  (l, :, 2) = Hu2;
                    g.zb (l, :)    = Zb2;
                    g.fix(l, 1)    = true;
                end
%             end
            Hu1            = Hum(r, 1);
            Zb1            =-H_m(r, 1)+Z_1;
            g.x  (r, :, 1) = Z_1;
            g.x  (r, :, 2) = Hu1;
            g.zb (r, :)    = Zb1;
            g.fix(r, 1)    = true;
            %--------------------------------------------------------------

            if plot_
                figure;
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
                close all;
            end
        end
    end
end

H_dofnew = g.x(:, :, 1)-g.zbinit;

xx = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSITIVITY-PRESERVING LIMITER: NEEDS TO BE FIXED
%--------------------------------------------------------------------------
% Z_dof            = g.x(:, :, 1);
% Zbdof            = g.zb;
% H_dof            = Z_dof-Zbdof;
% log_w            = true (K, 1);
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


% Z_dof = g.x(:, :, 1);
% Hudof = g.x(:, :, 2);
% Zbdof = g.zb;
% H_dof = Z_dof-Zbdof;
% 
% if any(H_dof < 0, 'all')
%     xx = 1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end