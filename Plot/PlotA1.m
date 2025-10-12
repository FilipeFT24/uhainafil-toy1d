function [obj] = PlotA1(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color = linspecer(9, 'Qualitative');
LW    = 2.5;
%--------------------------------------------------------------------------
h0    = g.data.h0;
xg    = g.data.xg;
xv    = g.data.xv;
CFL   = g.CFL;
K     = g.numE;
N     = g.N;
BF    = g.BF;
test  = g.test;
n     = size(xg, 1);
xc    = zeros(K, 1);
for i = 1:K
    xc(i, 1) = 1./2.*(xv(i, 1)+xv(i+1, 1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch test
    case 1
        %------------------------------------------------------------------
        xlim       = [xv(1, 1), xv(K+1, 1)]./h0;
        xtick      =-50:10:30;
        xticklabel = ["-50", "-40", "-30", "-20", "-10", "0", "10", "20", "30"];
        xtitle     = "$x^{\prime}$";
        xc         = xc./h0;
        %------------------------------------------------------------------
        switch g.data.opt
            case 1
                path       = sprintf("Plot/T01 - Grilli/DATA/0.10/C%d/Grilli_", g.data.bathy);
                ylim       = [-0.01, 0.20];
                ytick      = [-0.01, 0.00:0.05:0.20];
                yticklabel = ["-0.01", "0.00", "0.05", "0.10", "0.15", "0.20"];
                ytitle     = "$\zeta^{\prime}$";
            case 2
                path       = sprintf("Plot/T01 - Grilli/DATA/0.15/C%d/Grilli_", g.data.bathy);
                ylim       = [-0.01, 0.30];
                ytick      = [-0.01, 0.00:0.05:0.30];
                yticklabel = ["-0.01", "0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30"];
                ytitle     = "$\zeta^{\prime}$";
            case 3
                path       = sprintf("Plot/T01 - Grilli/DATA/0.20/C%d/Grilli_", g.data.bathy);
                ylim       = [-0.02, 0.40];
                ytick      = [-0.02, 0.00:0.10:0.40];
                yticklabel = ["-0.02", "0.00", "0.10", "0.20", "0.30", "0.40"];
                ytitle     = "$\zeta^{\prime}$";
            case 4
                path       = sprintf("Plot/T01 - Grilli/DATA/0.25/C%d/Grilli_", g.data.bathy);
                ylim       = [-0.02, 0.50];
                ytick      = [-0.02, 0.00:0.10:0.50];
                yticklabel = ["-0.02", "0.00", "0.10", "0.20", "0.30", "0.40", "0.50"];
                ytitle     = "$\zeta^{\prime}$";
            case 5
                path       = sprintf("Plot/T01 - Grilli/DATA/0.30/C%d/Grilli_", g.data.bathy);
                ylim       = [-0.05, 0.60];
                ytick      = [-0.05, 0.00:0.10:0.60];
                yticklabel = ["-0.05", "0.00", "0.10", "0.20", "0.30", "0.40", "0.50", "0.60"];
                ytitle     = "$\zeta^{\prime}$";
            case 6
                path       = sprintf("Plot/T01 - Grilli/DATA/0.40/C%d/Grilli_", g.data.bathy);
                ylim       = [-0.05, 0.80];
                ytick      = [-0.05, 0.00:0.20:0.80];
                yticklabel = ["-0.05", "0.00", "0.20", "0.40", "0.60", "0.80"];
                ytitle     = "$\zeta^{\prime}$";
            otherwise
                return
        end
        %------------------------------------------------------------------
    case 2
        %------------------------------------------------------------------
        path       = "Plot/T02 - LEGI/DATA/LEGI_";
        %------------------------------------------------------------------
    case 3
        %------------------------------------------------------------------
        path       = "Plot/T03 - Composite beach (NOAA)/DATA/Compbeach_";
        %------------------------------------------------------------------
        xlim       = [xv(1, 1), xv(K+1, 1)];
        xtick      =-15:5:20;
        xticklabel = ["-15", "-10", "-5", "0", "5", "10", "15", "20"];
        xtitle     = "$x$";
        %------------------------------------------------------------------
        switch g.data.opt
            case 1
                ylim       = [-0.01, 0.12];
                ytick      = [-0.01, 0.00:0.02:0.12];
                yticklabel = ["-0.01", "0.00", "0.02", "0.04", "0.06", "0.08", "0.10", "0.12"];
            case 2
                ylim       = [-0.10, 0.90];
                ytick      = [-0.10, 0.00:0.10:0.90];
                yticklabel = ["-0.10", "0.00", "0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70", "0.80", "0.90"];
            otherwise
                return
        end
        ytitle     = "$\zeta^{\prime}$";
        %------------------------------------------------------------------
    case 4
        %------------------------------------------------------------------
        xlim       = [-55, 20];
        xtick      = [-55, -50:10:20];
        xticklabel = ["-55", "50", "-40", "-30", "-20", "-10", "0", "10", "20"];
        xtitle     = "$x^{\prime}$";
        %------------------------------------------------------------------
        switch g.data.opt
            case 1
                path       = "Plot/T04 - Reflection of shoaling waves/DATA/0.1000/Ref_";
                ylim       = [-0.02, 0.30];
                ytick      = [-0.02, 0.00:0.05:0.30];
                yticklabel = ["-0.02", "0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30"];
            case 2
                path       = "Plot/T04 - Reflection of shoaling waves/DATA/0.1834/Ref_";
                ylim       = [-0.05, 0.70];
                ytick      = [-0.05, 0.00:0.10:0.70];
                yticklabel = ["-0.05", "0.00", "0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70"];
            otherwise
                return
        end
        ytitle     = "$\zeta^{\prime}$";
        %------------------------------------------------------------------
    case 5
        %%%
    otherwise
        return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_ = sum(pagemtimes(g.x(:, :, 1), g.Wbf'), 2)./sum(g.W, 2);
Z_ = Z_./h0;
%--------------------------------------------------------------------------
Ph = cell(1, 2);
%--------------------------------------------------------------------------
figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
hold on;
PlotX1(...
    xlim, ylim, ...
    xtick, ytick, ...
    xticklabel, yticklabel, ...
    xtitle, ytitle);
Ph{1, 2} = plot(xc, Z_, 'Color', Color(2, :), 'LineWidth', LW);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str2 = split(sprintf('%.15g', CFL), '.');
num3 = '';
%--------------------------------------------------------------------------
switch test
    case 1
        switch g.data.opt
            case 1
                num1 = "_10";
            case 2
                num1 = "_15";
            case 3
                num1 = "_20";
            case 4
                num1 = "_25";
            case 5
                num1 = "_30";
            case 6
                num1 = "_40";
            otherwise
                return
        end
    case 2
        switch g.data.opt
            case 1
                num1 = "_0960";
            case 2
                num1 = "_2975";
            case 3
                num1 = "_4560";
            case 4
                num1 = "_5343";
            otherwise
                return
        end
    case 3
        switch g.data.opt
            case 1
                num1 = "_039";
            case 2
                num1 = "_264";
            otherwise
                return
        end
    case 4
        switch g.data.opt
            case 1
                num1 = "_1000";
            case 2
                num1 = "_1834";
            otherwise
                return
        end
    case 5
        %%%
    otherwise
        return
end
if numel(str2) < 2
    num2 = '0000';
else
    num2 = [str2{2}, '0000'];
    num2 = num2(1, 1:2);
end
if g.data.alpha > 1+eps
    num3 = "_159";
end
%--------------------------------------------------------------------------
fid = sprintf("%sP%d%s_%s_GN%s.mat", path, g.p, num1, num2, num3);
if isfile(fid)
    delete(fid);
end
nf  = g.data.nf+1;
ndt = g.data.tend./g.data.nf;
obj = Output(Ph, fid, ndt, nf, K, N, n, xc, xv, xg, BF);
%--------------------------------------------------------------------------
obj.Write1(0, g.x, g.WD);
if test == 1 && g.data.opt == 3
    phi = zeros(K, N);
else
    phi = nan;
end
obj.Write2(g.t, g.x, phi, g.WD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end