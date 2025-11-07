function [obj] = PlotA2(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color = linspecer(9, 'Qualitative');
grey  = repmat(0.80, 1, 3);
LW    = [1.0, 2.5];
%--------------------------------------------------------------------------
K     = g.numE;
N     = g.N;
test  = g.test;
xv    = g.data.xv;
xc    = zeros(K, 1);
for i = 1:K
    xc(i, 1) = 1./2.*(xv(i, 1)+xv(i+1, 1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch test
    case 1
        path        = "Plot/S0/T0 - Convergence (Green-Nagdhi)/DATA/C_";
        xlim        = [0, 200];
        xtick       = 0:50:200;
        xticklabel  = ["0", "50", "100", "150", "200"];
        ylim        = [-0.01, 0.25];
        ytick       = [-0.01, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25];
        yticklabel  = ["-0.01", "0.00", "0.05", "0.10", "0.15", "0.20", "0.25"];
        xtitle      = "$x$";
        ytitle      = "$\zeta^{\prime}$";
        legendlabel = ["$\zeta_{a}^{\prime}$", "$\zeta_{h}^{\prime}$"];
        numcols     = 1;
    otherwise
        return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ha = g.data.H(xc, 0);
Za = g.data.Z(xc, 0);
H_ = sum(pagemtimes(g.x(:, :, 1)-g.zbinit, g.Wbf'), 2)./sum(g.W, 2);
Z_ = sum(pagemtimes(g.x(:, :, 1)         , g.Wbf'), 2)./sum(g.W, 2);
h0 = g.data.h0;
Za = Za./h0;
Z_ = Z_./h0;
%--------------------------------------------------------------------------
ph = cell(1, 2);
Ph = cell(1, 2);
figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
hold on;
PlotX1(...
    xlim, ylim, ...
    xtick, ytick, ...
    xticklabel, yticklabel, ...
    xtitle, ytitle);
ph{1, 1} = plot(nan, nan, 'Color', Color(1, :), 'LineWidth', LW(1, 1));
ph{1, 2} = plot(nan, nan, 'Color', Color(2, :), 'LineWidth', LW(1, 2));
Ph{1, 1} = plot(xc , Za , 'Color', Color(1, :), 'LineWidth', LW(1, 1));
Ph{1, 2} = plot(xc , Z_ , 'Color', Color(2, :), 'LineWidth', LW(1, 2));
PlotX2(...
    ph(1, :), legendlabel, 'northeast', numcols);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CFL  = g.CFL;
str2 = split(sprintf('%.15g', CFL), '.');
%--------------------------------------------------------------------------
switch test
    case 1
    otherwise
        return
end
num1 = '';
if numel(str2) < 2
    num2 = '0000';
else
    num2 = [str2{2}, '0000'];
    num2 = num2(1, 1:2);
end
%--------------------------------------------------------------------------
fid = sprintf('%sP%d%s_%s_GN.mat', path, g.p, num1, num2);
if isfile(fid)
    delete(fid);
end
nf  = g.data.nf+1;
ndt = g.data.tend./g.data.nf;
obj = Output(Ph, fid, ndt, nf, K, N, 0, xc, nan, nan, nan);
%--------------------------------------------------------------------------
obj.Write1(0, g.x, 0);
obj.Write2(0, g.x, nan, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end