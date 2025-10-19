function [g] = wetDryCellAverage(g)
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
%--------------------------------------------------------------------------
Z_dof  = g.x(:, :, 1);
Hudof  = g.x(:, :, 2);
Zbdof  = g.zbinit;
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

H_aux = H_dof-drytol;

if any(H_dof < drytol, 'all')
    xx = 1;
end

if g.nit == 95429
    xx = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Loop through "problematic" cells:
%--------------------------------------------------------------------------
log1 = H_l(:, 1) <= drytol & H_r(:, 1) <= drytol & H_m(:, 1) <= drytol;
c    = find(log1);
m    = size(c, 1);
%--------------------------------------------------------------------------
for i = 1:m
    %----------------------------------------------------------------------
    % VARS:
    o              = c  (i, 1);
    Zbo            = Zbm(o, 1);
    %----------------------------------------------------------------------
    % VARS MODIFIED:
    H_o            = drytol-rtol;
    Z_o            = H_o+Zbo;
    Huo            = Hum(o, 1).*(H_m(o, 1)./H_o);
    %----------------------------------------------------------------------
    % SOLUTION:
    H_dof(o, :)    = H_o;
    Hudof(o, :)    = Huo;
    Z_dof(o, :)    = Z_o;
    H_m  (o, 1)    = H_o;
    Hum  (o, 1)    = Huo;
    Z_m  (o, 1)    = Z_o;
    %
    g.x  (o, :, 1) = Z_o;
    g.x  (o, :, 2) = Huo;
    g.zb (o, :)    = Zbo;
    H_l  (o, 1)    = H_o;
    H_r  (o, 1)    = H_o;
    Z_l  (o, 1)    = Z_o;
    Z_r  (o, 1)    = Z_o;
    Zbl  (o, 1)    = Zbo;
    Zbr  (o, 1)    = Zbo;
    g.fix(o, 1)    = true;
    %----------------------------------------------------------------------
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through "problematic" faces:
%--------------------------------------------------------------------------
log2            = false(Kf, 1);
log2(1:Kf-1, 1) = log2(1:Kf-1, 1) | H_l(1:K, 1) <= drytol;
log2(2:Kf  , 1) = log2(2:Kf  , 1) | H_r(1:K, 1) <= drytol;
f               = find(log2);
n               = size(f, 1);
%--------------------------------------------------------------------------
for i = 1:n
    %----------------------------------------------------------------------
    r    = f(i, 1);
    l    = r-1;
    ldry = H_l(l, 1) <= drytol & H_r(l, 1) <= drytol;
    rdry = H_l(r, 1) <= drytol & H_r(r, 1) <= drytol;
    if ldry && rdry
        continue
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    if ldry && H_r(r, 1) > drytol
        %------------------------------------------------------------------
        % DRY/WET
        %------------------------------------------------------------------
        d = l-1;
        w = r-1;
        x = r+1;
 
%         figure;
%         hold on;
%         plot(reshape(Z_dof(1:71, :)'+drytol, [], 1), '-ob');
%         plot(reshape(Zbdof(1:71, :)'       , [], 1), '-*k');


        %precisas é de corrigir a célula a que pertence a face que tem
        %problemas!

        if Z_r(w, 1) > Zbl(w, 1)+drytol
            Z_w = Z_m(w, 1);
        else
            Z_w = Z_l(x, 1);
            if Z_w < Z_r(d, 1)
                Z_d            = Z_w;
                Hud            = Hum(d, 1);
                Zbd            =-H_m(d, 1)+Z_d;
                g.x  (d, :, 1) = Z_d;
                g.x  (d, :, 2) = Hud;
                g.zb (d, :)    = Zbd;
                g.fix(d, 1)    = true;
            end
        end
        Huw            = Hum(w, 1);
        Zbw            =-H_m(w, 1)+Z_w;
        g.x  (w, :, 1) = Z_w;
        g.x  (w, :, 2) = Huw;
        g.zb (w, :)    = Zbw;
        g.fix(w, 1)    = true;
        %------------------------------------------------------------------
    elseif rdry && H_l(l, 1) > drytol
        %------------------------------------------------------------------
        % WET/DRY
        %------------------------------------------------------------------
        d = r;
        w = l;
        x = l-1;
        if Z_l(w, 1) > Zbr(w, 1)+drytol
            Z_w = Z_m(w, 1);
        else
            Z_w = Z_r(x, 1);
            if Z_w < Z_l(d, 1)
                Z_d            = Z_w;
                Hud            = Hum(d, 1);
                Zbd            =-H_m(d, 1)+Z_d;
                g.x  (d, :, 1) = Z_d;
                g.x  (d, :, 2) = Hud;
                g.zb (d, :)    = Zbd;
                g.fix(d, 1)    = true;
            end
        end
        Huw            = Hum(w, 1);
        Zbw            =-H_m(w, 1)+Z_w;
        g.x  (w, :, 1) = Z_w;
        g.x  (w, :, 2) = Huw;
        g.zb (w, :)    = Zbw;
        g.fix(w, 1)    = true;
        %------------------------------------------------------------------
    else

        %         %------------------------------------------------------------------
%         if H_r(l, 1) < drytol && H_l(r, 1) < drytol
% 
%             Z_1            = Z_m(l, 1);
%             Hu1            = Hum(l, 1);
%             Zb1            =-H_m(l, 1)+Z_1;
%             Z_2            = Z_m(r, 1);
%             Hu2            = Hum(r, 1);
%             Zb2            =-H_m(r, 1)+Z_2;
% 
%             g.x  (l, :, 1) = Z_1;
%             g.x  (l, :, 2) = Hu1;
%             g.x  (r, :, 1) = Z_2;
%             g.x  (r, :, 2) = Hu2;
%             g.zb (l, :)    = Zb1;
%             g.zb (r, :)    = Zb2;
%             g.fix(l, 1)    = true;
%             g.fix(r, 1)    = true;
% 
%             if Z_l(l, 1) <= Zbr(l, 1)+drytol
%                 xx = 1;
%             end
%             if Z_l(r, 1) <= Zbr(r, 1)+drytol
%                 xx = 1;
%             end
%             xx = 1;
%         else
%             if H_r(l, 1) < drytol
%                 p = l;
%             else
%                 p = r;
%             end
% 
%             if Z_l(p, 1) <= Zbr(p, 1)+drytol
%                 xx = 1;
%             end
%             if n > 2
%                 yy = 1;
%             end
%  
%             Z_p            = Z_m(p, 1);
%             Hup            = Hum(p, 1);
%             Zbp            =-H_m(p, 1)+Z_p;
%             g.x  (p, :, 1) = Z_p;
%             g.x  (p, :, 2) = Hup;
%             g.zb (p, :)    = Zbp;
%             g.fix(p, 1)    = true;
%         end
%         %------------------------------------------------------------------
    end
    %----------------------------------------------------------------------
end

Z_dof = g.x(:, :, 1);
Zbdof = g.zb;
H_dof = Z_dof-Zbdof;

if any(H_dof < drytol, 'all')
    xx = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSITIVITY-PRESERVING LIMITER: NEEDS TO BE FIXED
%--------------------------------------------------------------------------
%{
log_w            =~g.fix & H_m > tol;
H_dofw           = H_dof(log_w, :);
Hudofw           = Hudof(log_w, :);
Zbdofw           = Zbdof(log_w, :);
H_mw             = H_m  (log_w, 1);
Humw             = Hum  (log_w, 1);
%--------------------------------------------------------------------------
H_quad           = H_dofw*fc';
min_H            = min(H_quad, [], 2);
theta            = min(1, (H_mw-drytol)./(H_mw-min_H)); % uhaina uses 0 instead of drytol
theta            = max(0, theta);
%--------------------------------------------------------------------------
g.x(log_w, :, 1) = H_mw+theta.*(H_dofw-H_mw)+Zbdofw;
g.x(log_w, :, 2) = Humw+theta.*(Hudofw-Humw);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



% function [g] = wetDryCellAverage(g)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Z_dof      = g.x(:, :, 1);
% Zbdof      = g.zb;
% H_dof      = Z_dof-Zbdof;
% is_dry     = all(H_dof < eps, 2);
% is_wet     = all(H_dof > eps, 2);
% is_partial =~is_dry & ~is_wet;
% K = g.numE;
% g.WD        = zeros(K, 1);
% if g.p > 0
%     g.WD(is_wet    , 1) = 0; % wet
%     g.WD(is_partial, 1) = 1; % wet and dry
%     g.WD(is_dry    , 1) = 2; % dry
% end
% 
% 
% N           = g.N;
% is_partial  = g.WD == 1;
% nb_1        = sum(is_partial, 1);
% Z_dof       = g.x(:, :, 1);
% Hudof       = g.x(:, :, 2);
% Zbdof       = g.zb;
% H_dof       = Z_dof-Zbdof;
% Z_bar       = meanval(g, Z_dof);
% Hubar       = meanval(g, Hudof);
% Zbbar       = meanval(g, Zbdof);
% H_bar       = meanval(g, H_dof);
% Z_dof_1     = Z_dof(is_partial, :);
% Zbdof_1     = Zbdof(is_partial, :);
% H_dof_1     = H_dof(is_partial, :);
% Z_bar_1     = Z_bar(is_partial, 1);
% Hubar_1     = Hubar(is_partial, 1);
% Zbbar_1     = Zbbar(is_partial, 1);
% H_bar_1     = H_bar(is_partial, 1);
% %--------------------------------------------------------------------------
% g_Z_dof     = zeros(nb_1, N);
% g_Zbdof     = zeros(nb_1, N);
% Zbdofmax    = max  (Zbdof_1, [], 2);
% %{
% Z_dofmaxwet = zeros(nb_1, 1);
% for i = 1:nb_1
%     Z_dofmaxwet(i, 1) = max(Z_dof_1(i, H_dof_1(i, :) > g.drytol), [], 2);
% end
% %}
% flag             = H_dof_1 > 0;%g.drytol;
% Z_dofmask        = Z_dof_1;
% Z_dofmask(~flag) =-realmax;
% Z_dofmaxwet      = max(Z_dofmask, [], 2);
% is_partial1      = Z_dofmaxwet <= Zbdofmax;%+g.drytol;
% is_partial2      =~is_partial1;
% 
% % is_partial1 = true(nb_1, 1);
% % is_partial2 = false(nb_1, 1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g_Z_dof(is_partial1, :) = repmat(Z_dofmaxwet(is_partial1, 1)                        , [1, N]); % #1: Z  = ZwetMAX
% g_Zbdof(is_partial1, :) = repmat(Z_dofmaxwet(is_partial1, 1)-H_bar_1(is_partial1, 1), [1, N]); % #1: Zb = ZwetMAX-Hbar
% %--------------------------------------------------------------------------
% g_Z_dof(is_partial2, :) = repmat(Z_bar_1(is_partial2, 1), [1, N]); % #2: Z  = Zbar
% g_Zbdof(is_partial2, :) = repmat(Zbbar_1(is_partial2, 1), [1, N]); % #2: Zb = Zbbar
% %--------------------------------------------------------------------------
% % Overwrite zeta, hu and zb:
% g.x (is_partial, :, 1) = g_Z_dof;
% g.x (is_partial, :, 2) = repmat(Hubar_1, [1, N]);
% g.zb(is_partial, :)    = g_Zbdof;
% 
% g.fix = g.WD == 1;
% 
% g = ZhangAndShuLimiter(g);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end