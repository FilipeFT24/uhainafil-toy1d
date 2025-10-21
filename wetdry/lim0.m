function [g] = lim0(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS #1:
%--------------------------------------------------------------------------
K     = g.numE;
Kf    = K+1;
fl    = g.bfl;
fr    = g.bfr;
rtol  = eps;
g.fix = false(K, 1);
%--------------------------------------------------------------------------
Z_dof = g.x(:, :, 1);
Hudof = g.x(:, :, 2);
Zbdof = g.zbinit;
H_dof = Z_dof-Zbdof;
H_m   = meanval(g, H_dof);
Hum   = meanval(g, Hudof);
Z_l   = Z_dof*fl';
Z_r   = Z_dof*fr';
Zbl   = Zbdof*fl';
Zbr   = Zbdof*fr';
Zbmax = max(Zbl(2:K, 1), Zbr(1:K-1, 1));
H_l   = Z_l-Zbl;
H_r   = Z_r-Zbr;
%
H_l   = Z_l(2:K  , 1)-Zbmax;
H_r   = Z_r(1:K-1, 1)-Zbmax;
%

% figure;
% hold on;
% %plot(H_l-H_r, '-ob');
% plot(Zbl(2:K  , 1)-Zbr(1:K-1, 1), '-*r');
% 
% HH = min(Z_l(2:K  , 1), Z_r(1:K-1, 1))-Zbmax;
% 
% 
% xx = H_l-H_r;
% 
% yy = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through "problematic" faces:
%--------------------------------------------------------------------------
n = 0;
if g.p > 0
    log1 = abs(H_l-H_r) > rtol;
    f    = find(log1);
    n    = size(f, 1);
end
%--------------------------------------------------------------------------
for i = 1:n
    %----------------------------------------------------------------------
    F = f(i, 1)+1;
    r = F;
    l = r-1;
    %----------------------------------------------------------------------
    if Z_r(l, 1) > Z_l(r, 1)
        p   = l;
        Z_p = Z_l(r, 1);
    else
        p   = r;
        Z_p = Z_r(l, 1);
    end
    Hup            = Hum(p, 1);
    Zbp            =-H_m(p, 1)+Z_p;
    g.x  (p, :, 1) = Z_p;
    g.x  (p, :, 2) = Hup;
    g.zb (p, :)    = Zbp;
    g.fix(p, 1)    = true;
    Z_l  (p, 1)    = Z_p;
    Z_r  (p, 1)    = Z_p;
    %----------------------------------------------------------------------
end


Z_dof = g.x(:, :, 1);
Zbdof = g.zb;
H_dof = Z_dof-Zbdof;

xx = 1;


% % isto é capaz de resolver o problema disto:         if log2l(r, 1) % modify wet cell instead
%             p = p+1;
%             d = d+1;
%             x = x+1;
%         end
% 
% numa dry/wet é sempre a dry. numa dry/dry nao sei...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end