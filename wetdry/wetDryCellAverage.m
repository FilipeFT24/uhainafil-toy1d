function [g] = wetDryCellAverage(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS:
%--------------------------------------------------------------------------
K                = g.numE;
N                = g.N;
fc               = g.bf;
fl               = g.bfl;
fr               = g.bfr;
tol              = g.drytol; % round-off tolerance
%--------------------------------------------------------------------------
Z_dof            = g.x(:, :, 1);
Hudof            = g.x(:, :, 2);
Zbdof            = g.zbinit;
H_dof            = Z_dof-Zbdof;
Hum              = meanval(g, Hudof);
H_m              = meanval(g, H_dof);
Z_l              = Z_dof*fl';
Z_r              = Z_dof*fr';
Zbl              = Zbdof*fl';
Zbr              = Zbdof*fr';
H_l              = H_dof*fl';
H_r              = H_dof*fr';
%--------------------------------------------------------------------------
log_wd           = H_m(1:K-1, 1) > tol & H_m(2:K, 1) < tol; % wet/dry
log_dw           = H_m(1:K-1, 1) < tol & H_m(2:K, 1) > tol; % dry/wet
wet1             = find(log_wd);
wet2             = find(log_dw)+1;
g.wt_wd(wet1, 1) = true;
g.wt_dw(wet2, 1) = true;
log_w            = H_m > tol & ~(g.wt_wd | g.wt_dw);
n1               = size(wet1, 1);
n2               = size(wet2, 1);
Z_1              = zeros(n1, 1);
Z_2              = zeros(n2, 1);
%--------------------------------------------------------------------------
% DEBUG:

% aa = (H_dof(:, 1) < eps & H_dof(:, 2) > eps) | (H_dof(:, 1) > eps & H_dof(:, 2) < eps);
% if sum(aa, 1) > 1
%     xx = 1;
% end
% bb = H_m(1:K-1, 1)-H_m(2:K, 1);
% if any(bb(dry2:dry2+30, 1) > 0, 1)
%     wrong = 1;
% end

if n1 > 1
    dw = 1;
end
if n2 > 1
    wd = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WET/DRY TRANSITION:
%--------------------------------------------------------------------------
if n1 > 0
    %----------------------------------------------------------------------
    % VARS:
    l_1               = Z_l(wet1, 1) <= Z_r(wet1-1, 1) | (Z_l(wet1, 1) > Z_r(wet1-1, 1) & Zbl(wet1, 1) < Zbr(wet1-1, 1));
    Z_1( l_1, 1)      = Z_r(wet1( l_1, 1)-1, 1);
    Z_1(~l_1, 1)      = Z_l(wet1(~l_1, 1)  , 1);
    Hu1               = Hum(wet1, 1);
    Zb1               =-H_m(wet1, 1)+Z_1;
    %----------------------------------------------------------------------
    % P0 fix:
    g.x  (wet1, :, 1) = repmat(Z_1, 1, N);
    g.x  (wet1, :, 2) = repmat(Hu1, 1, N);
    g.zb (wet1, :)    = repmat(Zb1, 1, N);
    g.fix(wet1, 1)    = true;
    %----------------------------------------------------------------------
    % BATHY fix:
    lb_1              = Z_l(wet1+1, 1) < Z_1;
    dry1              = wet1(lb_1, 1)+1;
    g.x  (dry1, :, 1) = repmat(Z_1(lb_1, 1), 1, N);
    g.x  (dry1, :, 2) = 0;
    g.zb (dry1, :)    = repmat(Z_1(lb_1, 1), 1, N);
    g.fix(dry1, 1)    = true;
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRY/WET TRANSITION:
%--------------------------------------------------------------------------
if n2 > 0
    %----------------------------------------------------------------------
    % VARS:
    l_2               = Z_r(wet2, 1) <= Z_l(wet2+1, 1) | (Z_r(wet2, 1) > Z_l(wet2+1, 1) & Zbr(wet2, 1) < Zbl(wet2+1, 1));
    Z_2( l_2, 1)      = Z_l(wet2( l_2, 1)+1, 1);
    Z_2(~l_2, 1)      = Z_r(wet2(~l_2, 1)  , 1);
    Zb2               =-H_m(wet2, 1)+Z_2;
    Hu2               = Hum(wet2, 1);
    %----------------------------------------------------------------------
    % P0 fix:
    g.x  (wet2, :, 1) = repmat(Z_2, 1, N);
    g.x  (wet2, :, 2) = repmat(Hu2, 1, N);
    g.zb (wet2, :)    = repmat(Zb2, 1, N);
    g.fix(wet2, 1)    = true;
    %----------------------------------------------------------------------
    % BATHY fix:
    lb_2              = Z_r(wet2-1, 1) < Z_2;
    dry2              = wet2(lb_2, 1)-1;
    g.x  (dry2, :, 1) = repmat(Z_2(lb_2, 1), 1, N);
    g.x  (dry2, :, 2) = 0;
    g.zb (dry2, :)    = repmat(Z_2(lb_2, 1), 1, N);
    g.fix(dry2, 1)    = true;
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSITIVITY-PRESERVING LIMITER:
%--------------------------------------------------------------------------
Hudofw           = Hudof(log_w, :);
Zbdofw           = Zbdof(log_w, :);
H_dofw           = H_dof(log_w, :);
Humw             = Hum  (log_w, 1);
H_mw             = H_m  (log_w, 1);
%--------------------------------------------------------------------------
H_quad           = H_dofw*fc';
min_H            = min(H_quad, [], 2);
theta            = min(1, (H_mw-tol)./(H_mw-min_H)); % uhaina uses 0 instead of drytol
theta            = max(0, theta);
%--------------------------------------------------------------------------
g.x(log_w, :, 1) = H_mw+theta.*(H_dofw-H_mw)+Zbdofw;
g.x(log_w, :, 2) = Humw+theta.*(Hudofw-Humw);
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
% flag             = H_dof_1 > g.drytol;
% Z_dofmask        = Z_dof_1;
% Z_dofmask(~flag) =-realmax;
% Z_dofmaxwet      = max(Z_dofmask, [], 2);
% is_partial1      = Z_dofmaxwet <= Zbdofmax+g.drytol;
% is_partial2      =~is_partial1;
% 
% is_partial1 = true(nb_1, 1);
% is_partial2 = false(nb_1, 1);
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




% %
% function [g] = wetDryCellAverage(g)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N           = g.N;
% drytol      = g.drytol;
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
% H_bar(H_bar < drytol) = 0;
% %--------------------------------------------------------------------------
% Z_dof_1     = Z_dof(is_partial, :);
% Zbdof_1     = Zbdof(is_partial, :);
% H_dof_1     = H_dof(is_partial, :);
% Z_bar_1     = Z_bar(is_partial, 1);
% Hubar_1     = Hubar(is_partial, 1);
% Zbbar_1     = Zbbar(is_partial, 1);
% H_bar_1     = H_bar(is_partial, 1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g_Z_dof     = zeros(nb_1, N);
% g_Zbdof     = zeros(nb_1, N);
% Zbdofmax    = max  (Zbdof_1, [], 2);
% %{
% Z_dofmaxwet = zeros(nb_1, 1);
% for i = 1:nb_1
%     Z_dofmaxwet(i, 1) = max(Z_dof_1(i, H_dof_1(i, :) > drytol), [], 2);
% end
% %}
% flag             = H_dof_1 > drytol;
% Z_dofmask        = Z_dof_1;
% Z_dofmask(~flag) =-realmax;
% Z_dofmaxwet      = max(Z_dofmask, [], 2);
% is_partial1      = Z_dofmaxwet <= Zbdofmax;%+drytol;
% is_partial2      =~is_partial1;
% 
% % is_partial1 = true(nb_1, 1);
% % is_partial2 = false(nb_1, 1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g_Z_dof(is_partial1, :) = repmat(Z_dofmaxwet(is_partial1, 1)                        , [1, N]); % #1: Z  = ZwetMAX
% g_Zbdof(is_partial1, :) = repmat(Z_dofmaxwet(is_partial1, 1)-H_bar_1(is_partial1, 1), [1, N]); % #1: Zb = ZwetMAX-Hbar
% %--------------------------------------------------------------------------
% % g_Z_dof(is_partial2, :) = repmat(Zbdofmax(is_partial2, 1)+H_bar_1(is_partial2, 1), [1, N]); % #1: Z  = ZwetMAX
% % g_Zbdof(is_partial2, :) = repmat(Zbdofmax(is_partial2, 1), [1, N]); % #1: Zb = ZwetMAX-Hbar
% 
% %
% g_Z_dof(is_partial2, :) = repmat(Z_bar_1(is_partial2, 1), [1, N]); % #2: Z  = Zbar
% g_Zbdof(is_partial2, :) = repmat(Zbbar_1(is_partial2, 1), [1, N]); % #2: Zb = Zbbar
% %
% 
% if any(is_partial2, 1)
%     diff1 = Z_bar_1(is_partial2, 1)-(Zbdofmax(is_partial2, 1)+H_bar_1(is_partial2, 1))
% 
%     xx = 1;
% end
% 
% %--------------------------------------------------------------------------
% % Overwrite zeta, hu and zb:
% g.x(is_partial, :, 1) = g_Z_dof;
% g.x(is_partial, :, 2) = ones(1, N).*Hubar_1;
% g.zb(is_partial, :)    = g_Zbdof;
% %--------------------------------------------------------------------------
% g.fix = g.WD == 1;
% 
% g = ZhangAndShuLimiter(g);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
% %
