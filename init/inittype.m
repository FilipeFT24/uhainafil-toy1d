function [F_dof] = inittype(type, f, xydc, xyqc, Wfi)
switch type
    case 0
        % INTERPOLATION (NODAL!):
        F_dof = f(xydc);
    case 1
        % PROJECTION:
        F_dof = f(xyqc)*Wfi;
    otherwise
        return
end
end