function [] = PlotS0T0(g, export)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





Color  = linspecer(9, 'Qualitative');
grey   = repmat(0.80, 1, 3);






M      = ["-o", "-s"];
lw     = 1.5;
LW     = 2.5;
flag   = [1, 1, 1];
path1  = "Plot/S0/T0 - Convergence (SW)/DATA";
path2  = "Plot/S0/T0 - Convergence (SW)";
path3  = "Plot/S0/T0 - Convergence (SW)/";


xm     = 0;
xM     = 1;
K      = 5000;
xv     = linspace(xm, xM, K+1)';
xc     = zeros(K, 1);
for i = 1:K
    xc(i, 1) = 1./2.*(xv(i, 1)+xv(i+1, 1));
end



n   = 3;
v   = 2;
Fh1 = zeros(6, 1+v, n);
Fh2 = zeros(K, v);
P1  = cell (v, n);
P2  = cell (v, 1);
%--------------------------------------------------------------------------
for i = 1:n
    Fh1(:, :, i) = readmatrix(sprintf("%s/P%d_CONVERGENCE.txt", path1, i));
    for j = 2:6
        Fh1(j, 1, i) = Fh1(1, 1, i)./Fh1(j, 1, i);
        for k = 2:1+v
            Fh1(j, k, i) = Fh1(j, k, i)./Fh1(1, k, i);
        end
    end
    Fh1(1, :, i) = 1;
end
%--------------------------------------------------------------------------
aux = geth0(1);
fh  = load(sprintf("%s/S00_5000_P1.mat", path1));
for i = 1:v
    Fh2(:, i) = sum(pagemtimes(fh.uh(:, :, i), aux.Wbf'), 2)./sum(aux.W, 2);
end
Fh3 = abs(fh.ec);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim1        = [0.50, 1.00];
xtick1       = 0.50:0.10:1.00;
xticklabel1  = ["0.50", "0.60", "0.70", "0.80", "0.90", "1.00"];
ylim1        = [0.05, 1.00];
ytick1       = [0.05, 0.10, 0.25, 0.50, 1.00];
yticklabel1  = ["0.05", "0.10", "0.25", "0.50", "1.00"];
xtitle1      = "$h/h_{0}$";
ytitle1      = "$\|e/e_{0}\|_{\ell_{2}}$";
legendlabel1 = ["$k=1$", "$k=2$", "$k=3$"];
numcols1     = 1;
fid1         = ["SW (z, 0.01s).pdf", "SW (q, 0.01s).pdf"];
%--------------------------------------------------------------------------
xlim2        = [0.00, 1.00];
xtick2       = 0.00:0.25:1.00;
xticklabel2  = ["0.00", "0.25", "0.50", "0.75", "1.00"];
ylim21       = [0, 5];
ytick21      = 0:1:5;
yticklabel21 = ["0", "1", "2", "3", "4", "5"];
ylim22       = [0, 60];
ytick22      = 0:10:600;
yticklabel22 = ["0", "10", "20", "30", "40", "50", "60"];
ytitle21     = "$\zeta$";
ytitle22     = "$q$";
ylim2        = {ylim21, ylim22};
ytick2       = {ytick21, ytick22};
yticklabel2  = {yticklabel21, yticklabel22};
xtitle2      = "$x$";
ytitle2      = {ytitle21, ytitle22};
fid2         = ["SW_P1 (z, 0.01s).pdf", "SW_P1 (q, 0.01s).pdf"];
%--------------------------------------------------------------------------
xlim3        = [0.00, 1.00];
xtick3       = 0.00:0.25:1.00;
xticklabel3  = ["0.00", "0.25", "0.50", "0.75", "1.00"];
ylim3        = 10.^[-9, -3];
ytick3       = 10.^(-9:-3);
yticklabel3  = ["-9", "-8", "-7", "-6", "-5", "-4", "-3"];
xtitle3      = "$x$";
ytitle3      = "$log(|e_{c}|)$";
legendlabel3 = ["$\zeta$", "$q$"];
numcols3     = 1;
fid3         = "SW_P1 (ec, 0.01s).pdf";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1)
    for i = 1:v
        %------------------------------------------------------------------
        fig1 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        hold on;
        %------------------------------------------------------------------
        hi = 1./4;
        p  = 2:1:4;
        A  = exp(1).^(linspace(0, log(0.5), 100));
        B  = [A(60), A(41)];
        Z  = {@(x) hi.*x.^2, @(x) hi.*x.^3, @(x) hi.*x.^4};
        for j = 1:n
            fplot(@(x) Z{1, j}(x), '-', 'Color', Color(j, :), 'LineWidth', 0.50);
            plot(...
                [B(1, 1), B(1, 1); B(1, 1), B(1, 1); B(1, 2), B(1, 1)], [Z{1, j}(B(1, 1)), Z{1, j}(B(1, 2)); Z{1, j}(B(1, 2)), Z{1, j}(B(1, 2)); Z{1, j}(B(1, 2)), Z{1, j}(B(1, 1))], ...
                'Color', Color(j, :), 'LineWidth', 0.50);
            text(min(B, [], 2)-5.0e-03, sqrt(Z{1, j}(B(1, 1)).*Z{1, j}(B(1, 2))), num2str(p(1, j)), 'Color', 'k', 'FontSize', 19.5);
        end
        %------------------------------------------------------------------
        for j = 1:n
            P1{i, j} = plot(Fh1(:, 1, j), Fh1(:, 1+i, j), M(1, i), 'Color', Color(j, :), 'MarkerFaceColor', Color(j, :), 'LineWidth', LW);
        end
        PlotX1(...
            xlim1, ylim1, ...
            xtick1, ytick1, ...
            xticklabel1, yticklabel1, ...
            xtitle1, ytitle1);
        set(gca, ...
            'XScale', 'log', ...
            'YScale', 'log', ...
            'XDir','reverse');
        grid on;
        set(gca, ...
            'GridLineStyle', '-', ...
            'GridAlpha', 0.15, ...
            'MinorGridLineStyle', ':', ...
            'MinorGridAlpha', 0.10, ...
            'XMinorGrid', 'on', ...
            'YMinorGrid', 'on');
        PlotX2(...
            P1(i, :), legendlabel1, 'northeast', numcols1);
        %------------------------------------------------------------------
        if export
            exportgraphics(fig1, fid1(1, i), 'ContentType', 'Vector', 'Resolution', 600);
            movefile(fid1(1, i), path3);
        end
        %------------------------------------------------------------------
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 2)
    for i = 1:v
        %------------------------------------------------------------------
        fig2 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        hold on;
        if i == 1
            xb = sum(g.coordV0T, 2)./2;
            yb = meanval(g, g.zbinit);
            patch(...
                xb, yb, grey, 'EdgeColor', grey);
        end
        %------------------------------------------------------------------
        P2{i, 1} = plot(xc, Fh2(:, i), '-', 'Color', Color(2, :), 'LineWidth', LW);
        %------------------------------------------------------------------
        PlotX1(...
            xlim2, ylim2{1, i}, ...
            xtick2, ytick2{1, i}, ...
            xticklabel2, yticklabel2{1, i}, ...
            xtitle2, ytitle2{1, i});
        %------------------------------------------------------------------
        if export
            exportgraphics(fig2, fid2(1, i), 'ContentType', 'Vector', 'Resolution', 600);
            movefile(fid2(1, i), path3);
        end
        %------------------------------------------------------------------
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 3)
    %----------------------------------------------------------------------
    fig3 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
    hold on;
    %----------------------------------------------------------------------
    for i = 1:v
        P2{i, 2} = plot(xc, Fh3(:, i), '-', 'Color', Color(i, :), 'LineWidth', LW);
    end
    PlotX1(...
        xlim3, ylim3, ...
        xtick3, ytick3, ...
        xticklabel3, yticklabel3, ...
        xtitle3, ytitle3);
    set(gca, ...
        'YScale', 'log');
    grid on;
    set(gca, ...
        'GridLineStyle', '-', ...
        'GridAlpha', 0.15, ...
        'MinorGridLineStyle', ':', ...
        'MinorGridAlpha', 0.10, ...
        'XMinorGrid', 'on', ...
        'YMinorGrid', 'on');
    PlotX2(...
        P2(:, 2)', legendlabel3, 'northeast', numcols3);
    %----------------------------------------------------------------------
    if export
        exportgraphics(fig3, fid3, 'ContentType', 'Vector', 'Resolution', 600);
        movefile(fid3, path3);
    end
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end