function [g] = wetDryCellAverage(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS #1:
%--------------------------------------------------------------------------
K      = g.numE;
Kf     = K+1;
fc     = g.bf;
fl     = g.bfl;
fr     = g.bfr;
rtol   = 1e-10;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through "problematic" cells:
%--------------------------------------------------------------------------
log1 = H_l(:, 1) < drytol & H_r(:, 1) < drytol & H_m(:, 1) < drytol;
c    = find(log1);
m    = size(c, 1);
%--------------------------------------------------------------------------
% for i = 1:m
%     o              = c(i, 1);
%     Zbo            = Zbdof(o, :);
%     Z_o            = Zbo;
%     H_o            = 0;
%     Huo            = 0;
%     H_dof(o, :)    = H_o;
%     H_m  (o, 1)    = H_o;
%     Hudof(o, :)    = Huo;
%     Hum  (o, 1)    = Huo;
%     g.x  (o, :, 1) = Z_o;
%     g.x  (o, :, 2) = Huo;
%     g.zb (o, :)    = Zbo;
%     H_l  (o, 1)    = H_o;
%     H_r  (o, 1)    = H_o;
%     Z_l  (o, 1)    = Z_o*fl';
%     Z_r  (o, 1)    = Z_o*fr';
%     Zbl  (o, 1)    = Zbo*fl';
%     Zbr  (o, 1)    = Zbo*fr';
%     g.fix(o, 1)    = true;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through "problematic" faces:
%--------------------------------------------------------------------------
log2            = false(Kf, 1);
log2(1:Kf-1, 1) = log2(1:Kf-1, 1) | H_l(1:K, 1) < drytol;
log2(2:Kf  , 1) = log2(2:Kf  , 1) | H_r(1:K, 1) < drytol;
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
    if ldry && H_r(r, 1) > drytol
        %------------------------------------------------------------------
        % DRY/WET
        %------------------------------------------------------------------
        d = l;
        w = r;
        if Z_r(w, 1) > Zbl(w, 1)+drytol %max(Zbl(w, 1), Zbr(d, 1))+drytol
            Z_w = Z_m(w, 1);
        else
            Z_w = Z_l(w+1, 1);
            if Zbr(d, 1) < Z_w
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
        if Z_l(w, 1) > Zbr(w, 1)+drytol %max(Zbr(w, 1), Zbl(d, 1))+drytol
            Z_w = Z_m(w, 1);
        else
            Z_w = Z_r(w-1, 1);
            if Zbl(d, 1) < Z_w
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
        %------------------------------------------------------------------
        if H_l(r, 1) < drytol
            p = r;
            w = l;
        else
            p = l;
            w = r;
        end



        xx = 1;
        %------------------------------------------------------------------
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