function [g] = wetDryCellAverage(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS:
%--------------------------------------------------------------------------
K      = g.numE;
N      = g.N;
fc     = g.bf;
fl     = g.bfl;
fr     = g.bfr;
tol1   = eps;
tol2   = g.drytol;
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
%--------------------------------------------------------------------------
log_dr = H_r(1:K-1, 1) < tol1;
log_dl = H_l(2:K  , 1) < tol1;
%--------------------------------------------------------------------------
idl    = find(log_dl);
idr    = find(log_dr)+1;
ndl    = size(idl, 1);
ndr    = size(idr, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ?/Dry face:
%--------------------------------------------------------------------------
for i = 1:ndl
    %----------------------------------------------------------------------
    w = idl(i, 1);
    d = w+1;
    if H_m(w, 1) < tol1
        continue
    end
    %----------------------------------------------------------------------
    if Z_r(w-1, 1) > Zbl(d, 1)%+tol2
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
            Z_l  (d, 1)    = Z_d;
            Z_r  (d, 1)    = Z_d;
            Zbl  (d, 1)    = Zbd;
            Zbr  (d, 1)    = Zbd;
        end
    end
    Huw            = Hum(w, 1);
    Zbw            =-H_m(w, 1)+Z_w;
    g.x  (w, :, 1) = Z_w;
    g.x  (w, :, 2) = Huw;
    g.zb (w, :)    = Zbw;
    g.fix(w, 1)    = true;
    Z_l  (w, 1)    = Z_w;
    Z_r  (w, 1)    = Z_w;
    Zbl  (w, 1)    = Zbw;
    Zbr  (w, 1)    = Zbw;
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dry/? face:
%--------------------------------------------------------------------------
for i = 1:ndr
    %----------------------------------------------------------------------
    w = idr(i, 1);
    d = w-1;
    if H_m(w, 1) < tol1
        continue
    end
    %----------------------------------------------------------------------
    if Z_l(w+1, 1) > Zbr(d, 1)
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
            %{
            Z_l  (d, 1)    = Z_d;
            Z_r  (d, 1)    = Z_d;
            Zbl  (d, 1)    = Zbd;
            Zbr  (d, 1)    = Zbd;
            %}
        end
    end
    Huw            = Hum(w, 1);
    Zbw            =-H_m(w, 1)+Z_w;
    g.x  (w, :, 1) = Z_w;
    g.x  (w, :, 2) = Huw;
    g.zb (w, :)    = Zbw;
    g.fix(w, 1)    = true;
    %{
    Z_l  (w, 1)    = Z_w;
    Z_r  (w, 1)    = Z_w;
    Zbl  (w, 1)    = Zbw;
    Zbr  (w, 1)    = Zbw;
    %}
    %----------------------------------------------------------------------
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
% N           = g.N;
% is_partial  = g.WD == 1;
% nb_1        = sum(is_partial, 1);
% Z_dof       = g.x(:, :, 1);
% Hudof       = g.x(:, :, 2);
% Zbdof       = g.zbinit;
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