function [] = PlotT15C(g, export)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color = linspecer(9, 'Qualitative');
grey  = repmat(0.80, 1, 3);
nc    = 50;
ns    = [10, 8];
nf    = g.data.nf+1;
R1    = linspace(1, Color(1, 1), nc);
G1    = linspace(1, Color(1, 2), nc);
B1    = linspace(1, Color(1, 3), nc);
R2    = linspace(1, Color(2, 1), nc);
G2    = linspace(1, Color(2, 2), nc);
B2    = linspace(1, Color(2, 3), nc);
lred  = [R1', G1', B1'];
lblue = [R2', G2', B2'];
M     = ["o", ":", "-", "-.", "--"];
lw    = 1.5;
LW    = 2.5;
flag  = [1, 0];
path1 = "Plot/T15 - Synolakis/DATA";
path2 = "Plot/T15 - Synolakis";
path3 = "Plot/T15 - Synolakis/";
%--------------------------------------------------------------------------
G     = g.data.G;
h0    = g.data.h0;
slope = 1./19.85;
xz    =-20;
xM    = 80;
K     = 1000;
xv    = linspace(xz, xM, K+1)';
xc    = zeros(K, 1);
for i = 1:K
    xc(i, 1) = 1./2.*(xv(i, 1)+xv(i+1, 1));
end
%--------------------------------------------------------------------------
m     = 2;
n     = 3;
Fh1   = zeros(nf, K, 1+n, m, 2);
th1   = zeros(nf, 1+n, m, 2);
Fe    = cell (m, 1);
str1  = ["019", "040"];
str2  = {...
    ["synt25.dat", "synt40.dat", "synt55.dat", "synt65.dat"];
    ["synt20.dat", "synt26.dat", "synt32.dat", "synt38.dat", "synt44.dat", "synt50.dat", "synt56.dat", "synt62.dat"]};
for i = 1:m
    Fe{i, 1} = cell(1, size(str2{i, 1}, 2));
    for j = 1:size(str2{i, 1}, 2)
        Fe{i, 1}{1, j}       = readmatrix(sprintf("%s/0.%s/%s", path2, str1(1, i), str2{i, 1}(1, j)));
        Fe{i, 1}{1, j}(:, 2) = Fe{i, 1}{1, j}(:, 2)./h0;
    end
    for j = 1:1+n
        Fh1(:, :, j, i, 1) = getdata1(sprintf("%s/0.%s/S_P%d_%s_SW.mat", path1, str1(1, i), j-1, str1(1, i)), j-1, 1);
        th1(:, j, i, 1)    = getdata3(sprintf("%s/0.%s/S_P%d_%s_SW.mat", path1, str1(1, i), j-1, str1(1, i)));
        %{
        if j ~= 1 % EXCLUDE P0
            Fh1(:, :, j, i, 2) = getdata1(sprintf("%s/0.%s/S_P%d_%s_GN.mat", path1, str1(1, i), j-1, str1(1, i)), j-1, 1);
            th1(:, j, i, 2)    = getdata3(sprintf("%s/0.%s/S_P%d_%s_GN.mat", path1, str1(1, i), j-1, str1(1, i)));
        end
        %}
    end
end
Fh1 = Fh1./h0;
th1 = th1./sqrt(h0./G);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim11        = [-2, 20];
xtick11       =-2:2:20;
xticklabel11  = ["-2", "0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20"];
xlim12        = [-5, 20];
xtick12       =-5:5:20;
xticklabel12  = ["-5", "0", "5", "10", "15", "20"];
xlim2         = [-20, 80];
xtick2        =-20:10:80;
xticklabel2   = ["-20", "-10", "0", "10", "20", "30", "40", "50", "60", "70", "80"];
%--------------------------------------------------------------------------
ylim11        = [-0.04, 0.10];
ytick11       =-0.04:0.02:0.10;
yticklabel11  = ["-0.04", "-0.02", "0.00", "0.02", "0.04", "0.06", "0.08", "0.10"];
ylim12        = [-0.05, 0.20];
ytick12       =-0.05:0.05:0.20;
yticklabel12  = ["-0.05", "0.00", "0.05", "0.10", "0.15", "0.20"];
%--------------------------------------------------------------------------
legendlabel1  = ["$\textrm{SW}$", "$\textrm{GN}$", "", "", "$k=0$", "$k=1$", "$k=2$", "$k=3$"];
legendlabel2  = ["$\epsilon = 0.019$", "$\epsilon = 0.040$"];
legendlabel3  = ["$t^{\prime} = t^{\prime}_{0}$", "$t^{\prime} = t^{\prime}_{p}$", "$t^{\prime} = t^{\prime}_{f}$"];
numcols1      = 1;
xtitle1       = "$x$";
ytitle1       = "$z^{\prime}$";
%--------------------------------------------------------------------------
xlim1         = {xlim11; xlim12};
xtick1        = {xtick11; xtick12};
xticklabel1   = {xticklabel11; xticklabel12};
ylim1         = {ylim11; ylim12};
ytick1        = {ytick11; ytick12};
yticklabel1   = {yticklabel11; yticklabel12};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlimz11       = [-3, 2];
xtickz11      =-3:1:2;
xticklabelz11 = ["-3", "-2", "-1", "0", "1", "2"];
xlimz12       = [-3, 2];
xtickz12      =-3:1:2;
xticklabelz12 = ["-3", "-2", "-1", "0", "1", "2"];
%--------------------------------------------------------------------------
ylimz11       = [-0.02, 0.18];
ytickz11      = [-0.02, 0.00:0.05:0.15, 0.18];
yticklabelz11 = ["-0.02", "0.00", "0.05", "0.10", "0.15", "0.18"];
ylimz12       = [-0.04, 0.15];
ytickz12      = [-0.04, 0.00:0.05:0.15];
yticklabelz12 = ["-0.04", "0.00", "0.05", "0.10", "0.15"];
%--------------------------------------------------------------------------
xm1           = (xlim11(1, 1)+xlim11(1, 2))./2;
xm2           = (xlim12(1, 1)+xlim12(1, 2))./2;
xf            = 1.025./2;
xz11          = (xlimz11-xm1).*xf+xm1;
xz12          = (xlimz12-xm2).*xf+xm2;
xr11          = ([4, 11.5]-xm1).*xf+xm1;
xr12          = ([4, 11.5]-xm2).*xf+xm2;
%--------------------------------------------------------------------------
yz11          = ylimz11;
yz12          = ylimz12;
yr11          = [0.050, 0.225];
yr12          = [0.075, 0.240];
%--------------------------------------------------------------------------
xlimz1        = {xlimz11; xlimz12};
xtickz1       = {xtickz11; xtickz12};
xticklabelz1  = {xticklabelz11; xticklabelz12};
ylimz1        = {ylimz11; ylimz12};
ytickz1       = {ytickz11; ytickz12};
yticklabelz1  = {yticklabelz11; yticklabelz12};
xz1           = {xz11; xz12};
yz1           = {yz11; yz12};
xr1           = {xr11; xr12};
yr1           = {yr11; yr12};
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P1            = cell(1, 2*(1+n));
P2            = cell(m, 3);
tg            = {...
    [25, 40, 55, 70];
    [20, 26, 32, 38, 44, 50, 56, 62]};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = m
        for j = 1:size(str2{i, 1}, 2)
            %--------------------------------------------------------------
            fig1 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
            cp1 = get(gca, 'Position');
            hold on;
            %--------------------------------------------------------------
            patch(...
                [xlim1{i, 1}(1, 1), xlim1{i, 1}(1, 1), -1./slope.*ylim1{i, 1}(1, 1)], [ylim1{i, 1}(1, 1), -slope.*xlim1{i, 1}(1, 1), ylim1{i, 1}(1, 1)], grey, 'EdgeColor', grey);
            %--------------------------------------------------------------
            plot(Fe{i, 1}{1, j}(:, 1), Fe{i, 1}{1, j}(:, 2), M(1, 1), 'Color', Color(1, :), 'LineWidth', lw, 'MarkerFaceColor', 'w', 'MarkerSize', 10);
            for k = 1:1+n
                [~, idx] = ...
                    min(abs(th1(:, k, i, 1)-tg{i, 1}(1, j)));
                plot(xc, Fh1(idx, :, k, i, 1)', M(1, 1+k), 'Color', Color(k, :), 'LineWidth', LW);
            end
            %--------------------------------------------------------------
            P1{1, 1} = plot(nan, nan, ':', 'Color', 'k', 'LineWidth', LW);
            P1{1, 2} = plot(nan, nan, '-', 'Color', 'k', 'LineWidth', LW);
            for k = 3:4
                P1{1, k} = plot(nan, nan, '-', 'Color', 'w', 'LineWidth', LW);
            end


            %--------------------------------------------------------------
            PlotX1(...
                xlim1{i, 1}, ylim1{i, 1}, ...
                xtick1{i, 1}, ytick1{i, 1}, ...
                xticklabel1{i, 1}, yticklabel1{i, 1}, ...
                xtitle1, ytitle1);
%             PlotX2(...
%                 P1(1, :), legendlabel1{i, 1}, 'northeast', numcols1);
            %--------------------------------------------------------------
%             MagInset(gcf, -1, ...
%                 [xz1{i, 1}, yz1{i, 1}], ...
%                 [xr1{i, 1}, yr1{i, 1}], {'SE', 'SE'; 'SW', 'SW'}, cp1, xlimz1{i, 1}, ylimz1{i, 1}, xtickz1{i, 1}, ytickz1{i, 1}, xticklabelz1{i, 1}, yticklabelz1{i, 1}, 10.5);
            %--------------------------------------------------------------
            if export(1, 1)
                exportgraphics(fig1, fid1(j, k), 'ContentType', 'Vector', 'Resolution', 600);
                movefile      (fid1(j, k), path3);
            end
            %--------------------------------------------------------------
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if flag(1, 2)
    %----------------------------------------------------------------------
    for j = 1:m
        %------------------------------------------------------------------
        if ~export(1, 1)
            fig21 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        else
            fig21 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
        end
        hold on;
        patch(...
            [xlim2(1, 1), xlim2(1, 1), -19.85.*ylim2(1, 1)], [ylim2(1, 1), -xlim2(1, 1)./19.85, ylim2(1, 1)], grey, 'EdgeColor', grey);
        %------------------------------------------------------------------
        plot(xc, Fh1(   1, :, 1, j)', ':', 'Color', Color(j, :), 'LineWidth', LW);
        plot(xc, Fh1(1000, :, 1, j)', '-', 'Color', Color(j, :), 'LineWidth', LW);
        plot(xc, Fh1( end, :, 1, j)', '-', 'Color', Color(j, :), 'LineWidth', LW);
        for k = 1:3
            P2{j, k} = plot(nan, nan, '-', 'Color', Color(j, :), 'LineWidth', LW);
        end
        %------------------------------------------------------------------
        PlotX1(...
            xlim2, ylim2, ...
            xtick2, ytick2, ...
            xticklabel2, yticklabel2, ...
            xtitle1, ytitle1);
        PlotX2(...
            P2(1, :), legendlabel3, 'northeast', numcols1);
        %------------------------------------------------------------------
        if export(1, 1)
            exportgraphics(fig21, "aa", 'ContentType', 'Vector', 'Resolution', 600);
            movefile      ("aa", path3);
        end
        %------------------------------------------------------------------
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end