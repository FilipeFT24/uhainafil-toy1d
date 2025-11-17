function [obj] = PlotA2(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color  = linspecer(9, 'Qualitative');
grey   = repmat(0.80, 1, 3);
LW     = [1.0, 2.5];
%--------------------------------------------------------------------------
K      = g.numE;
N      = g.N;
flag   = 0; % Plot auxiliary vars.
test   = g.test;
wetdry = g.data.wetdry;
xv     = g.data.xv;
xc     = zeros(K, 1);
for i = 1:K
    xc(i, 1) = 1./2.*(xv(i, 1)+xv(i+1, 1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch test
    case 1
        %------------------------------------------------------------------
        path        = "Plot/S0/T0 - Convergence (SW)/DATA/S00";
        xlim        = [0.0, 1.0];
        xtick       = 0.0:0.2:1.0;
        xticklabel  = ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"];
        ylim        = [0, 5];
        ytick       = 0:1:5;
        yticklabel  = ["0", "1", "2", "3", "4", "5"];
        xtitle      = "$x$";
        ytitle      = "$\zeta$";
        legendlabel = ["$\zeta_{a}^{\prime}$", "$\zeta_{h}^{\prime}$"];
        numcols     = 1;
        %------------------------------------------------------------------
    case 2
        %------------------------------------------------------------------
        path        = "Plot/S0/T0 - Convergence (GN)/DATA/S00";
        xlim        = [-50, 50];
        xtick       =-50:25:50;
        xticklabel  = ["-50", "-25", "0", "25", "50"];
        ylim        = [-0.01, 0.25];
        ytick       = [-0.01, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25];
        yticklabel  = ["-0.01", "0.00", "0.05", "0.10", "0.15", "0.20", "0.25"];
        xtitle      = "$x$";
        ytitle      = "$\zeta^{\prime}$";
        legendlabel = ["$\zeta_{a}^{\prime}$", "$\zeta_{h}^{\prime}$"];
        numcols     = 1;
        %------------------------------------------------------------------
    case {7, 8}
        %------------------------------------------------------------------
        switch test
            case 7
                path       = "Plot/S2/T1 - Head-on-colission/DATA/Hoc_";
                xlim       = [-100, 100];
                xtick      =-100:50:100;
                xticklabel = ["-100", "-50", "0", "50", "100"];
            case 8
                path       = "Plot/S2/T2 - Wall/DATA/Wall_";
                xlim       = [-100, 0];
                xtick      =-100:25:0;
                xticklabel = ["-100", "-75", "-50", "-25", "0"];
            otherwise
                return
        end
        switch g.data.opt
            case 1
                ylim       = [-0.45, 2.65];
                ytick      = [-0.45, 0.00:0.50:2.50, 2.65];
                yticklabel = ["-0.45", "0.00", "0.50", "1.00", "1.50", "2.00", "2.50", "2.65"];
            case 2
                ylim       = [-0.05, 0.45];
                ytick      = [-0.05, 0.00:0.10:0.40, 0.45];
                yticklabel = ["-0.05", "0.00", "0.10", "0.20", "0.30", "0.40", "0.45"];
            otherwise
                return
        end
        xtitle      = "$x$";
        ytitle      = "";
        legendlabel = ["$\zeta_{a}^{\prime}$", "$\zeta_{h}^{\prime}$"];
        numcols     = 1;
        %------------------------------------------------------------------
    case 9
        %------------------------------------------------------------------
        path        = "Plot/S2/T3 - Overtaking colission/DATA/Oc_";
        xlim        = [-200, 200];
        xtick       =-200:50:200;
        xticklabel  = ["-200", "-150", "-100", "-50", "0", "50", "100", "150", "200"];
        ylim        = [-0.10, 1.15];
        ytick       = [-0.10, 0.00:0.25:1.00, 1.15];
        yticklabel  = ["-0.10", "0.00", "0.25", "0.50", "0.75", "1.00", "1.15"];
        xtitle      = "$x$";
        ytitle      = "$\zeta^{\prime}$";
        legendlabel = ["$\zeta_{a}^{\prime}$", "$\zeta_{h}^{\prime}$"];
        numcols     = 1;
        %------------------------------------------------------------------
    case 10
        %------------------------------------------------------------------
        path        = "Plot/S2/T4 - Breakup of a Gaussian hump/DATA/Bgh_";
        xlim        = [-100, 100];
        xtick       =-100:50:100;
        xticklabel  = ["-100", "-50", "0", "50", "100"];
        %------------------------------------------------------------------
        switch g.data.opt
            case 1
                ylim       = [-0.25, 1.65];
                ytick      = [-0.25, 0.00:0.50:1.50, 1.65];
                yticklabel = ["-0.25", "0.00", "0.50", "1.00", "1.50", "1.65"];
            case 2
                ylim       = [-0.45, 2.10];
                ytick      = [-0.45, 0.00:0.50:2.10];
                yticklabel = ["-0.45", "0.00", "0.50", "1.00", "1.50", "2.00", "2.10"];
            otherwise
                return
        end
        xtitle      = "$x$";
        ytitle      = "$\zeta^{\prime}$";
        %------------------------------------------------------------------
    case 11
        %------------------------------------------------------------------
        path        = "Plot/S2/T5 - Dispersive dam-break/DATA/Ddb_";
        %------------------------------------------------------------------
        xlim        = [-300, 300];
        xtick       =-300:100:300;
        xticklabel  = ["-300", "-200", "-100", "0", "100", "200", "300"];
        ylim        = [0.95, 2.00];
        ytick       = [0.95, 1.00:0.25:2.00];
        yticklabel  = ["0.95", "1.00", "1.25", "1.50", "1.75", "2.00"];
        xtitle      = "$x$";
        ytitle      = "$h$";
        %------------------------------------------------------------------
    otherwise
        return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test == 11
    Xa = g.data.H(xc, 0);
    Xh = sum(pagemtimes(g.x(:, :, 1)-g.zbinit, g.Wbf'), 2)./sum(g.W, 2);
else
    h0 = g.data.h0;
    Xa = g.data.N(xc, 0)-h0;
    Xh = sum(pagemtimes(g.x(:, :, 1)-h0, g.Wbf'), 2)./sum(g.W, 2);
    Xa = Xa./h0;
    Xh = Xh./h0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test == 2 && flag
    %----------------------------------------------------------------------
    ph = cell(4, 2);
    Ph = cell(4, 2);
    for i = 1:2
        figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        for j = 1:2
            %--------------------------------------------------------------
            k = 2*(i-1)+j;
            %--------------------------------------------------------------
            subplot(1, 2, j);
            PlotX1(...
                xlim, ylim(k, :), ...
                xtick, ytick{k, 1}, ...
                xticklabel, yticklabel{k, 1}, ...
                xtitle, ytitle);
            if i ~= 1 || j ~= 1
                Xa = zeros(K, 1);
                Xh = zeros(K, 1);
            end
            ph{k, 1} = plot(nan, nan, 'Color', Color(1, :), 'LineWidth', LW(1, 1));
            ph{k, 2} = plot(nan, nan, 'Color', Color(2, :), 'LineWidth', LW(1, 2));
            Ph{k, 1} = plot(xc , Xa , 'Color', Color(1, :), 'LineWidth', LW(1, 1));
            Ph{k, 2} = plot(xc , Xh , 'Color', Color(2, :), 'LineWidth', LW(1, 2));
            PlotX2(...
                ph(k, :), legendlabel(k, :), 'northeast', numcols);
            %--------------------------------------------------------------
            grid on;
            grid minor;
            %--------------------------------------------------------------
        end
    end
    %----------------------------------------------------------------------
else
    %----------------------------------------------------------------------
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
    Ph{1, 1} = plot(xc , Xa , 'Color', Color(1, :), 'LineWidth', LW(1, 1));
    Ph{1, 2} = plot(xc , Xh , 'Color', Color(2, :), 'LineWidth', LW(1, 2));
    if ismembc(test, [1, 2, 7, 8, 9])
        PlotX2(...
            ph(1, :), legendlabel, 'northeast', numcols);
    end
    %----------------------------------------------------------------------
end
if test == 1 || wetdry
    xbathy = sum(g.coordV0T, 2)./2;
    ybathy = meanval(g, g.zbinit);
    patch(...
        xbathy, ybathy, grey, 'EdgeColor', grey);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CFL  = g.CFL;
str2 = split(sprintf('%.15g', CFL), '.');
%--------------------------------------------------------------------------
switch test
    case 1
    case 2
    case 7
    case 8
    case 9
    case 10
    case 11
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
fid = sprintf('%s_0512_P%d.mat', path, g.p);
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