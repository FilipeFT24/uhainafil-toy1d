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
% % Loop through "problematic" cells:
% %--------------------------------------------------------------------------
% log1 = H_l(:, 1) <= drytol & H_r(:, 1) <= drytol;
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
%     H_o            = 0;
%     Z_o            = H_o+Zbo;
%     Huo            = 0;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through "problematic" faces:
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
f               = find(log2);
n               = size(f, 1);
%----------------------------------------------------------------------

itstop = 5060;
if g.nit > itstop
    wrong = 1;
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
        aux = [1, 3, 4, 2];
    otherwise
        return
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
    if log2r(l, 1)
        if log2l(r, 1)
            %--------------------------------------------------------------
            % DRY/DRY
            %--------------------------------------------------------------
            H_L            = H_r(L, 1); %if H_L < rtol, Z_L = realmax; end
            H_R            = H_l(R, 1); %if H_R < rtol, Z_R = realmax; end
            Z_M            = min(Z_L, Z_R);
            g.x  (l, :, 1) = Z_M;
            g.x  (r, :, 1) = Z_M;
            g.x  (l, :, 2) = Hum(l, 1);
            g.x  (r, :, 2) = Hum(r, 1);
            g.zb (l, :)    =-H_m(l, 1)+Z_M;
            g.zb (r, :)    =-H_m(r, 1)+Z_M;
            g.fix(l, 1)    = true;
            g.fix(r, 1)    = true;
            %--------------------------------------------------------------

%             f1 = figure;
%             subplot(1, 2, 1);
%             hold on;
%             for j = l-1:r+1
%                 plot(g.xydc(j, aux), Z_dof(j, aux), '--ob');
%                 plot(g.xydc(j, aux), Zbdof(j, aux),  ':*k');
%             end
%             subplot(1, 2, 2);
%             hold on;
%             for j = l-1:r+1
%                 plot(g.xydc(j, aux), g.x (j, aux, 1), '--ob');
%                 plot(g.xydc(j, aux), g.zb(j, aux), ':*k');
%             end
%             close(f1);

            %--------------------------------------------------------------
        else
            %--------------------------------------------------------------
            % DRY/WET
            %--------------------------------------------------------------
            if Z_L < Z_R
                Z_X = Z_L;
            else
                Z_X = Z_R;
                if Z_r(l, 1) < Z_X
                    g.x  (l, :, 1) = Z_X;
                    g.x  (l, :, 2) = Hum(l, 1);
                    g.zb (l, :)    =-H_m(l, 1)+Z_X;
                    g.fix(l, 1)    = true;
                end
            end
            g.x  (r, :, 1) = Z_X;
            g.x  (r, :, 2) = Hum(r, 1);
            g.zb (r, :)    =-H_m(r, 1)+Z_X;
            g.fix(r, 1)    = true;
            %--------------------------------------------------------------

%             f2 = figure;
%             subplot(1, 2, 1);
%             hold on;
%             for j = l-1:r+1
%                 plot(g.xydc(j, aux), Z_dof(j, aux), '--ob');
%                 plot(g.xydc(j, aux), Zbdof(j, aux)+drytol,  ':*k');
%             end
%             subplot(1, 2, 2);
%             hold on;
%             for j = l-1:r+1
%                 plot(g.xydc(j, aux), g.x (j, aux, 1), '--ob');
%                 plot(g.xydc(j, aux), g.zb(j, aux)+drytol, ':*k');
%             end
%             close(f2);

        end
    else
        %------------------------------------------------------------------
        % WET/DRY
        %------------------------------------------------------------------
        if Z_L > Z_R
            Z_X = Z_R;
        else
            Z_X = Z_L;
            if Z_l(r, 1) < Z_X
                g.x  (r, :, 1) = Z_X;
                g.x  (r, :, 2) = Hum(r, 1);
                g.zb (r, :)    =-H_m(r, 1)+Z_X;
                g.fix(r, 1)    = true;
            end
        end
        g.x  (l, :, 1) = Z_X;
        g.x  (l, :, 2) = Hum(l, 1);
        g.zb (l, :)    =-H_m(l, 1)+Z_X;
        g.fix(l, 1)    = true;
        %------------------------------------------------------------------

%         f3 = figure;
%         subplot(1, 2, 1);
%         hold on;
%         for j = l-1:r+1
%             plot(g.xydc(j, aux), Z_dof(j, aux), '--ob');
%             plot(g.xydc(j, aux), Zbdof(j, aux)+drytol,  ':*k');
%         end
%         subplot(1, 2, 2);
%         hold on;
%         for j = l-1:r+1
%             plot(g.xydc(j, aux), g.x (j, aux, 1), '--ob');
%             plot(g.xydc(j, aux), g.zb(j, aux)+drytol, ':*k');
%         end
%         close(f3);

        %------------------------------------------------------------------
    end
    %----------------------------------------------------------------------
end


H_dofnew = g.x(:, :, 1)-g.zbinit-drytol;
a2 = find(H_dofnew(:, 1) < 0, 1, 'last');
%%p0
% if H_dofnew(a2-1, 1) > 0
%     xx = 1;
% end
%%p1
if H_dofnew(a2-1, 2) > 0
    xx = 1;
end

if any(H_dof < 0, 'all')
    xx = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end