function [obj] = PlotA2(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color  = linspecer(9, 'Qualitative');
grey   = repmat(0.80, 1, 3);
LW     = [1.0, 2.5];
%--------------------------------------------------------------------------
K      = g.numE;
p      = g.p;
h0     = g.data.h0;
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
        opt         = g.data.opt;
        switch opt
            case 1
                str = "0.21";
            case 2
                str = "0.96";
            otherwise
                return
        end
        %------------------------------------------------------------------
        switch test
            case 7
                path       = sprintf("Plot/S2/T1 - Head-on-colission/DATA/%s/Hoc", str);
                xlim       = [-100, 100];
                xtick      =-100:50:100;
                xticklabel = ["-100", "-50", "0", "50", "100"];
            case 8
                path       = sprintf("Plot/S2/T2 - Wall/DATA/%s/Wall", str);
                xlim       = [-100, 0];
                xtick      =-100:25:0;
                xticklabel = ["-100", "-75", "-50", "-25", "0"];
            otherwise
                return
        end
        switch opt
            case 1
                ylim       = [-0.05, 0.45];
                ytick      = [-0.05, 0.00:0.10:0.40, 0.45];
                yticklabel = ["-0.05", "0.00", "0.10", "0.20", "0.30", "0.40", "0.45"];
            case 2
                ylim       = [-0.45, 2.65];
                ytick      = [-0.45, 0.00:0.50:2.50, 2.65];
                yticklabel = ["-0.45", "0.00", "0.50", "1.00", "1.50", "2.00", "2.50", "2.65"];
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
        path        = "Plot/S2/T3 - Overtaking colission/DATA/Oc";
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
        opt         = g.data.opt;
        switch opt
            case 1
                str = "01";
            case 2
                str = "20";
            otherwise
                return
        end
        path       = sprintf("Plot/S2/T4 - Breakup of a Gaussian hump/DATA/%s/Bgh", str);
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
        path        = "Plot/S2/T5 - Dispersive dam-break/DATA/Ddb";
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
    case 12
        %------------------------------------------------------------------
        bathy       = g.data.bathy;
        xlim        = [xv(1, 1), xv(K+1, 1)]./h0;
        xtick       =-50:10:30;
        xticklabel  = ["-50", "-40", "-30", "-20", "-10", "0", "10", "20", "30"];
        xtitle      = "$x^{\prime}$";
        %------------------------------------------------------------------
        switch g.data.opt
            case 1
                path       = sprintf("Plot/S3/T1 - Grilli/DATA/0.10/C%d/Grilli", bathy);
                ylim       = [-0.01, 0.20];
                ytick      = [-0.01, 0.00:0.05:0.20];
                yticklabel = ["-0.01", "0.00", "0.05", "0.10", "0.15", "0.20"];
                ytitle     = "$\zeta^{\prime}$";
            case 2
                path       = sprintf("Plot/S3/T1 - Grilli/DATA/0.15/C%d/Grilli", bathy);
                ylim       = [-0.01, 0.30];
                ytick      = [-0.01, 0.00:0.05:0.30];
                yticklabel = ["-0.01", "0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30"];
                ytitle     = "$\zeta^{\prime}$";
            case 3
                path       = sprintf("Plot/S3/T1 - Grilli/DATA/0.20/C%d/Grilli", bathy);
                ylim       = [-0.02, 0.40];
                ytick      = [-0.02, 0.00:0.10:0.40];
                yticklabel = ["-0.02", "0.00", "0.10", "0.20", "0.30", "0.40"];
                ytitle     = "$\zeta^{\prime}$";
            case 4
                path       = sprintf("Plot/S3/T1 - Grilli/DATA/0.25/C%d/Grilli", bathy);
                ylim       = [-0.02, 0.50];
                ytick      = [-0.02, 0.00:0.10:0.50];
                yticklabel = ["-0.02", "0.00", "0.10", "0.20", "0.30", "0.40", "0.50"];
                ytitle     = "$\zeta^{\prime}$";
            otherwise
                return
        end
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
if test == 12
    xc = xc./h0;
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
            Ph{k, 1} = plot(xc , Xa , 'Color', Color(1, :), 'LineWidth', LW(1, 1));
            ph{k, 2} = plot(nan, nan, 'Color', Color(2, :), 'LineWidth', LW(1, 2));
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
    if ~ismembc(test, [12, 13, 14, 15, 16, 17])
        ph{1, 1} = plot(nan, nan, 'Color', Color(1, :), 'LineWidth', LW(1, 1));
        Ph{1, 1} = plot(xc , Xa , 'Color', Color(1, :), 'LineWidth', LW(1, 1));
    end
    ph{1, 2} = plot(nan, nan, 'Color', Color(2, :), 'LineWidth', LW(1, 2));
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
switch test
    case 1
    case 2
    case {7, 8, 9, 10, 11}
        fid = sprintf('%s_P%d_%d.mat', path, p, K);
    case 12
    otherwise
        return
end
if isfile(fid)
    delete(fid);
end
obj = Output(g, fid, Ph);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end