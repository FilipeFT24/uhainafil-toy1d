function [g] = dryLimiter(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drytol              = g.drytol;
veltol              = g.velcutoff;
Z_dof               = g.x(:, :, 1);
Hudof               = g.x(:, :, 2);
Zbdofold            = g.zb;
Zbdofnew            = g.zbinit;
H_dofold            = Z_dof-Zbdofold;
H_dofnew            = Z_dof-Zbdofnew;
H_mold              = meanval(g, H_dofold);
%--------------------------------------------------------------------------
isdry_drytol        = H_dofnew <= drytol;
isdry_veltol        = H_dofnew <= veltol;
isdry               = isdry_drytol | isdry_veltol;
Z_dof(isdry_drytol) = Zbdofnew(isdry_drytol);
Hudof(isdry       ) = 0;
%--------------------------------------------------------------------------
H_aux               = Z_dof-Zbdofnew;
H_maux              = meanval(g, H_aux);
% wet                 = find(g.wt_dw | g.wt_wd);
% nw                  = size(wet, 1);
% dH_m                = H_mold(wet, 1)-H_maux(wet, 1);
% for i = 1:nw
%     if dH_m(i, 1) < 0 && min(H_aux(wet(i, 1), :)+dH_m(i, 1)-drytol, [], 2) < 0
%         dH_m(i, 1) = 0;
%     end
%     Z_dof(wet(i, 1), :) = Z_dof(wet(i, 1), :)+dH_m(i, 1);
% end




% ADD TO REST (IF NEEDED)


%--------------------------------------------------------------------------
g.x(:, :, 1)        = Z_dof;
g.x(:, :, 2)        = Hudof;
g.zb                = g.zbinit;
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end