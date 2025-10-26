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
log1 = (H_l(:, 1) < drytol & H_r(:, 1) < drytol); 
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
for i = 2:Kf-1
    if log2(i, 1) && (log2l(i-1, 1) && log2r(i, 1)) % face surrounded by (externally) "dry" cells
        log2(i, 1) = false;
    end
end
log2( 1, 1)     = false;
log2(Kf, 1)     = false;
f               = find(log2);
n               = size(f, 1);
%--------------------------------------------------------------------------

wrong = false;
itstop = 150360000;
if g.nit > itstop
    wrong = true;
    
end

H_x = [H_l, H_r];
H_aux = [H_l, H_r]-drytol;


H_dofold = g.x(:, :, 1)-g.zbinit-drytol;
a2_old = find(H_dofold(:, 1) < 0, 1, 'last');
if n > 2
    stop = 1;
end
flag = 0;

%--------------------------------------------------------------------------
for i = 1:n
    %----------------------------------------------------------------------
    r   = f(i, 1);
    l   = r-1;
    L   = l-1;
    R   = r+1;
    Z_L = g.x (L, :, 1)*fr';
    Z_R = g.x (R, :, 1)*fl';
    ZbL = g.zb(L, :, 1)*fr';
    ZbR = g.zb(R, :, 1)*fl';
    H_L = Z_L-ZbL;
    H_R = Z_R-ZbR;
    Z_X = min(Z_L, Z_R);
    %----------------------------------------------------------------------
    if log2r(l, 1)
        if log2l(r, 1)
            %if Z_R < Z_L %H_L < rtol && H_R > rtol

               % if Zbl(r, 1) < Z_R
                    Z_X            = min(Z_L, Z_l(r, 1));
                    if Z_X > Z_L
                        zz = 1;
                    end

                    g.x  (l, :, 1) = Z_X;
                    g.x  (l, :, 2) = Hum(l, 1);
                    g.zb (l, :)    =-H_m(l, 1)+Z_X;
                    g.fix(l, 1)    = true;


%                     g.x  (r, :, 1) = Z_X;
%                     g.x  (r, :, 2) = Hum(r, 1);
%                     g.zb (r, :)    =-H_m(r, 1)+Z_X;
%                     g.fix(r, 1)    = true;
                %end



            %else
%                 Z_X            = Z_L;
%                 g.x  (r, :, 1) = Z_X;
%                 g.x  (r, :, 2) = Hum(r, 1);
%                 g.zb (r, :)    =-H_m(r, 1)+Z_X;
%                 g.fix(r, 1)    = true;
%                 if Zbr(l, 1) < Z_X
%                     g.x  (l, :, 1) = Z_X;
%                     g.x  (l, :, 2) = Hum(l, 1);
%                     g.zb (l, :)    =-H_m(l, 1)+Z_X;
%                     g.fix(l, 1)    = true;
% 
%                     flag = 0;
%                 end
            %end



%             if Zbl(r, 1) < Z_R
%                 flag = 0;
%                 xx = 1;
%             end
% 
%             Z_X = min(Z_L, Z_l(r, 1));
%             g.x  (l, :, 1) = Z_X;
%             g.x  (l, :, 2) = Hum(l, 1);
%             g.zb (l, :)    =-H_m(l, 1)+Z_X;
%             g.fix(l, 1)    = true;
        else
            %--------------------------------------------------------------
            % DRY/WET
            %--------------------------------------------------------------
            g.x  (r, :, 1) = Z_X;
            g.x  (r, :, 2) = Hum(r, 1);
            g.zb (r, :)    =-H_m(r, 1)+Z_X;
            g.fix(r, 1)    = true;
            if Z_r(l, 1) < Z_X
                g.x  (l, :, 1) = Z_X;
                g.x  (l, :, 2) = Hum(l, 1);
                g.zb (l, :)    =-H_m(l, 1)+Z_X;
                g.fix(l, 1)    = true;
            end
            %--------------------------------------------------------------
        end
    else
        if log2l(r, 1)
            %--------------------------------------------------------------
            % WET/DRY
            %--------------------------------------------------------------
            g.x  (l, :, 1) = Z_X;
            g.x  (l, :, 2) = Hum(l, 1);
            g.zb (l, :)    =-H_m(l, 1)+Z_X;
            g.fix(l, 1)    = true;
            if Z_l(r, 1) < Z_X
                g.x  (r, :, 1) = Z_X;
                g.x  (r, :, 2) = Hum(r, 1);
                g.zb (r, :)    =-H_m(r, 1)+Z_X;
                g.fix(r, 1)    = true;
            end
            %--------------------------------------------------------------
        end
    end
    PLOT(flag, g, l, r, Z_dof, Zbdof, drytol);
    %----------------------------------------------------------------------
end

Z_m = meanval(g, g.x(:, :, 1));

H_dofnew = g.x(:, :, 1)-g.zbinit-drytol;
%a2 = find(H_dofnew(:, 1) < 0, 1, 'last');

% if ~isempty(a2_old)
%     if H_dofnew(a2_old, 1) > 0 && H_dofold(a2_old, 1) < 0 && any(H_dofnew < 0, 'all')
%         xx = 1;
%     end
% end

%
%%p0
% if H_dofnew(a2-1, 1) > 0
%     xx = 1;
% end
%%p1

% if ~isempty(a1)
%     if ~all(H_dofnew(1:a1-1, 1) > 0, 1)
%         xx = 1;
%     end
% end
% if ~isempty(a2)
%     if H_dofnew(a2-2, 2) < 0
%         xx = 1;
%     end
% end


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
        plot(g.xydc(j, aux), Zbdof(j, aux),  ':*k');
    end
    subplot(1, 2, 2);
    hold on;
    for j = l-1:r+1
        plot(g.xydc(j, aux), g.x (j, aux, 1), '--ob');
        plot(g.xydc(j, aux), g.zbinit(j, aux), ':*k');
    end
    close(f1);
end
end