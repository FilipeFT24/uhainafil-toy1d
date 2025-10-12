function [] = PlotT03C(g, export)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color    = linspecer(9, 'Qualitative');
grey     = repmat(0.80, 1, 3);
MS       = 10;
M        = ["o", "-", "--", ":"];
LW       = 2.5;
flag     = [0, 1];
%--------------------------------------------------------------------------
G        = g.data.G;
h0       = g.data.h0;
sqrth0_G = sqrt(h0./G);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1)
    %----------------------------------------------------------------------
    fid1      = [...
        "compbeach_g04_Pk_039_05_GN.pdf", ...
        "compbeach_g05_Pk_039_05_GN.pdf", ...
        "compbeach_g06_Pk_039_05_GN.pdf", ...
        "compbeach_g07_Pk_039_05_GN.pdf", ...
        "compbeach_g08_Pk_039_05_GN.pdf", ...
        "compbeach_g09_Pk_039_05_GN.pdf", ...
        "compbeach_g10_Pk_039_05_GN.pdf";
        "compbeach_g04_Pk_264_05_GN.pdf", ...
        "compbeach_g05_Pk_264_05_GN.pdf", ...
        "compbeach_g06_Pk_264_05_GN.pdf", ...
        "compbeach_g07_Pk_264_05_GN.pdf", ...
        "compbeach_g08_Pk_264_05_GN.pdf", ...
        "compbeach_g09_Pk_264_05_GN.pdf", ...
        "compbeach_g10_Pk_264_05_GN.pdf"];
    path      = "Plot/T03 - Composite beach (NOAA)/DATA/";
    n         = 3;
    ng        = 7;
    Fe        = cell(2, 1);
    Fh        = cell(2, n);
    p1        = cell(2, ng, 3+n);
    Fe{1, 1}  = readmatrix('ts3a.txt');
    Fe{2, 1}  = readmatrix('ts3b.txt');
    Fh_data11 = load('Plot/T03 - Composite beach (NOAA)/DATA/compbeach_P1_039_05_GN2_159.mat'); Fh{1, 1} = Fh_data11.data;
    Fh_data12 = load('Plot/T03 - Composite beach (NOAA)/DATA/compbeach_P2_039_05_GN2_159.mat'); Fh{1, 2} = Fh_data12.data;
    Fh_data13 = load('Plot/T03 - Composite beach (NOAA)/DATA/compbeach_P3_039_05_GN2_159.mat'); Fh{1, 3} = Fh_data13.data;
    Fh_data21 = load('Plot/T03 - Composite beach (NOAA)/DATA/compbeach_P1_264_05_GN2_159.mat'); Fh{2, 1} = Fh_data21.data;
    Fh_data22 = load('Plot/T03 - Composite beach (NOAA)/DATA/compbeach_P2_264_05_GN2_159.mat'); Fh{2, 2} = Fh_data22.data;
    Fh_data23 = load('Plot/T03 - Composite beach (NOAA)/DATA/compbeach_P3_264_05_GN2_159.mat'); Fh{2, 3} = Fh_data23.data;
    for i = 1:2
        Fe{i, 1}(:, 1)      = (Fe{i, 1}(:, 1)-265)./sqrth0_G;
        Fe{i, 1}(:, 2:ng+1) = (Fe{i, 1}(:, 2:ng+1))./h0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xlim11       = [0, 200];
    xtick11      = 0:50:200;
    xticklabel11 = ["0", "50", "100", "150", "200"];
    ylim11       = [-0.01-eps, 0.08];
    ytick11      = [-0.01, 0:0.01:0.08];
    yticklabel11 = ["-0.01", "0.00", "0.01", "0.02", "0.03", "0.04", "0.05", "0.06", "0.07", "0.08"];
    %----------------------------------------------------------------------
    xlim12       = [0, 200];
    xtick12      = 0:50:200;
    xticklabel12 = ["0", "50", "100", "150", "200"];
    ylim12       = [-0.1-eps, 0.6];
    ytick12      = [-0.1, 0:0.1:0.6];
    yticklabel12 = ["-0.10", "0.00", "0.10", "0.20", "0.30", "0.40", "0.50", "0.60"];
    %----------------------------------------------------------------------
    xlim1        = {xlim11, xlim12};
    xtick1       = {xtick11, xtick12};
    xticklabel1  = {xticklabel11, xticklabel12};
    ylim1        = {ylim11, ylim12};
    ytick1       = {ytick11, ytick12};
    yticklabel1  = {yticklabel11, yticklabel12};
    %----------------------------------------------------------------------
    xtitle1      = "$t^{\prime}$";
    ytitle1      = "$\zeta^{\prime}$";
    numcols1     = 2;
    legendlabelg = ["$g_{4}$", "$g_{5}$", "$g_{6}$", "$g_{7}$", "$g_{8}$", "$g_{9}$", "$g_{10}$"];
    legendlabel1 = ["", "", "$k=1$", "$k=2$", "$k=3$"];
    %----------------------------------------------------------------------
    for i = 2
        for j = 1:ng
            %--------------------------------------------------------------
            fig1 = figure('Color', 'w', 'Windowstate', 'Maximized');
            hold on;
            for k = 2:3
                p1{i, j, k} = plot(nan, nan, M(1, 1), ...
                    'Color', 'w', ...
                    'LineWidth', LW(1, 1), ...
                    'MarkerFaceColor', 'w', ...
                    'MarkerSize', MS(1, 1));
            end
            for k = 4:6
                p1{i, j, k} = plot(nan, nan, M(1, k-2), ...
                    'Color', Color(k-2, :), ...
                    'LineWidth', LW(1, 2));
            end
            p1{i, j, 1} = plot(Fe{i, 1}(:, 1), Fe{i, 1}(:, 1+j), M(1, 1), 'Color', Color(1, :), 'LineWidth', LW(1, 1), 'MarkerFaceColor', 'w', 'MarkerSize', MS(1, 1));
            for k = 1:3
                %plot(Xh, Yh, Mh, 'Color', Color, 'LineWidth', LW);
                funch(Fh{i, k}(:, 1), Fh{i, k}(:, 1+j), Color(1+k, :), LW(1, 2), M(1, 1+k), xlim1{1, i}, ylim1{1, i}, xtick1{1, i}, ytick1{1, i}, xticklabel1{1, i}, yticklabel1{1, i}, xtitle1, ytitle1);
            end
            legend(cat(2, p1{i, j, :}), [legendlabelg(1, j), legendlabel1], ...
                'AutoUpdate', 'Off', ...
                'FontSize', 20.5, ...
                'Location', 'NorthEast', ...
                'NumColumns', numcols1);
            %--------------------------------------------------------------
            if export(1, 1)
                exportgraphics(fig1, fid1(i, j), 'ContentType', 'Vector', 'Resolution', 600);
                movefile      (fid1, path);
            end
            %--------------------------------------------------------------
        end
    end
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xm          =-15.04;
    xM          = 23.23;
    K           = 3827;
    xv          = linspace(xm, xM, K+1)';
    xc          = zeros(K, 1);
    for i = 1:K
        xc(i, 1) = 1./2.*(xv(i, 1)+xv(i+1, 1));
    end
    %----------------------------------------------------------------------
    path        = "Plot/T03 - Composite beach (NOAA)/";
    src         = [...
        join([path, "DATA/VIDEO_039.mat"], '');
        join([path, "DATA/VIDEO_264.mat"], '')];
    png         = [...
        join([path, "FRAMES/0.039/"], '');
        join([path, "FRAMES/0.264/"], '')];
    mp4         = [...
        join([path, "Compbeach_039"], '');
        join([path, "Compbeach_264"], '')];
    yc          = cell(2, 1);
    for i = 1:2
        YC       = load(src(i, 1));
        yc{i, 1} = YC.data1(:, 2:K+1);
    end
    nf          = size(yc{1, 1}, 1);
    ns          = 15;
    mod0        = round(nf./ns, 0);
    f1          = 0.10;
    f2          = 0.50;
    alpha       = linspace(f1, f2, ns);
    %----------------------------------------------------------------------
    aux.color   = zeros(ns, 4);
    for i = 1:ns
        aux.color(i, :) = [Color(2, :), alpha(1, i)];
    end
    aux.nf      = nf;
    aux.ns      = ns;
    aux.mod0    = mod0;
    aux.fps     = g.data.fps;
    aux.quality = 100;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    slope1      = 1./53;
    slope2      = 1./150;
    slope3      = 1./13;
    y0          =-1.00;
    y1          = y0+slope1./h0.*(19.40-15.04);
    y2          = y1+slope2./h0.*(22.33-19.40);
    y3          = y2+slope3./h0.*(23.23-22.33);
    I1          = [15.04, 19.40, 19.40, 15.04; -1.00, -1.00, y1, y0];
    I2          = [19.40, 22.33, 22.33, 19.40; -1.00, -1.00, y2, y1];
    I3          = [22.33, 23.23, 23.23, 22.33; -1.00, -1.00, y3, y2];
    %----------------------------------------------------------------------
    xlim        = [xv(1, 1), xv(K+1, 1)];
    xtick       =-15:5:20;
    xticklabel  = ["-15", "-10", "-5", "0", "5", "10", "15", "20"];
    %----------------------------------------------------------------------
    ylim1       = [-1.00, 0.15];
    ytick1      = [-1.00:0.20:0.00, 0.05:0.05:0.15];
    yticklabel1 = ["-1.00", "-0.80", "-0.60", "-0.40", "-0.20", "0.00", "0.05", "0.10", "0.15"];
    ylim2       = [-1.00, 1.00];
    ytick2      =-1.00:0.20:1.00;
    yticklabel2 = ["-1.00", "-0.80", "-0.60", "-0.40", "-0.20", "0.00", "0.20", "0.40", "0.60", "0.80", "1.00"];
    %----------------------------------------------------------------------
    ylim        = {ylim1; ylim2};
    ytick       = {ytick1; ytick2};
    yticklabel  = {yticklabel1; yticklabel2};
    %----------------------------------------------------------------------
    xtitle      = "$x$";
    ytitle      = "$z^{\prime}$";
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1
        %------------------------------------------------------------------
        if ~export(1, 2)
            figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        else
            figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
        end
        hold on;
        patch(I1(1, :), I1(2, :), grey, 'EdgeColor', grey);
        patch(I2(1, :), I2(2, :), grey, 'EdgeColor', grey);
        patch(I3(1, :), I3(2, :), grey, 'EdgeColor', grey);
        set(gca, ...
            'Box', 'on', ...
            'Clipping', 'on', ...
            'FontSize', 19.50, ...
            'Layer', 'top', ...
            'TickLabelInterpreter', 'latex', ...
            'XLim', xlim, ...
            'YLim', ylim{i, 1}, ...
            'XTick', xtick, ...
            'YTick', ytick{i, 1}, ...
            'XTickLabel', xticklabel, ...
            'YTickLabel', yticklabel{i, 1}, ...
            'XMinorTick', 'on', ...
            'YMinorTick', 'on');
        axis square;
        xlabel(xtitle, 'Fontsize', 25.00); xtickangle(0);
        ylabel(ytitle, 'Fontsize', 25.00); ytickangle(0);
        %------------------------------------------------------------------
        p1     = plot(xc, yc{i, 1}(1, :), '-', 'Color', Color(2, :), 'LineWidth', LW);
        p2     = copyobj(p1, gca);
        aux.yc = yc{i, 1};
        aux.p1 = p1;
        %------------------------------------------------------------------
        set(p2, 'Color', aux.color(1, :));
        %------------------------------------------------------------------
        PlotMP4(aux, png(i, 1), mp4(i, 1), export(1, 2));
        %------------------------------------------------------------------
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end