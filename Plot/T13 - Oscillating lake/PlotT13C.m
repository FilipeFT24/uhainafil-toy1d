function [] = PlotT13C(g, export)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color = linspecer(9, 'Qualitative');
grey  = repmat(0.80, 1, 3);
nc    = 50;
ns    = 5;
nf    = g.data.nf+1;
nk    = linspace(1, nf, ns);
R     = linspace(1, Color(2, 1), nc);
G     = linspace(1, Color(2, 2), nc);
B     = linspace(1, Color(2, 3), nc);
lblue = [R', G', B'];
M     = [":", ":", "-", "-.", "--"];
LW    = 2.5;
flag  = [1, 0];
path1 = "Plot/T13 - Oscillating lake/DATA (-12)";
path2 = "Plot/T13 - Oscillating lake";
path3 = "Plot/T13 - Oscillating lake/";
%--------------------------------------------------------------------------
xm    =-1.5;
xM    = 1.5;
K     = 3000;
xv    = linspace(xm, xM, K+1)';
xc    = zeros(K, 1);
for i = 1:K
    xc(i, 1) = 1./2.*(xv(i, 1)+xv(i+1, 1));
end
%--------------------------------------------------------------------------
m     = 2;
n     = 4;
p     = 3;
Fh1   = zeros(nf, K, n, m);
str1  = ["zeta", "hu"];
tend  = 2.*pi./sqrt(2.*0.1005);
tk    = [1./4, 1./2, 3./4, 1].*tend;
for i = 1:n
    for j = 1:m
        Fh1(:, :, i, j) = getdata1(sprintf("%s/Ol_P%d.mat", path1, i-1), i-1, j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim            = [-1.5, 1.5];
xtick           =-1.5:0.5:1.5;
xticklabel      = ["-1.5", "-1.0", "-0.5", "0.0", "0.5", "1.0", "1.5"];
%--------------------------------------------------------------------------
ylim1           = [0.00, 0.25];
ytick1          = 0.00:0.05:0.25;
yticklabel1     = ["\phantom{-}0.00", "0.05", "0.10", "0.15", "0.20", "0.25"];
ylim2           = [-10, 10].*10.^-3;
ytick2          = (-10:5:10).*10.^-3;
yticklabel2     = ["-10", "-5", "0", "5", "10"];
ylim            = {ylim1; ylim2};
yticklabel      = {yticklabel1; yticklabel2};
ytick           = {ytick1; ytick2};
%--------------------------------------------------------------------------
legendlabel     = ["$\mathrm{Analytical}$", "$k = 0$", "$k = 1$", "$k = 2$", "$k = 3$"];
numcols         = 1;
xtitle          = "$x$";
ytitle1         = "$\zeta$";
ytitle2         = "$hu$";
ytitle          = {ytitle1; ytitle2};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlimz11_1       = [-1.0025, -0.9775];
xtickz11_1      =-1.00:0.01:-0.98;
xticklabelz11_1 = ["-1.00", "-0.99", "-0.98"];
xlimz12_1       = [0.9775, 1.0025];
xtickz12_1      = 0.98:0.01:1.00;
xticklabelz12_1 = ["0.98", "0.99", "1.00"];
xlimz11_2       = [-1.105, -1.080];
xtickz11_2      =-1.10:0.01:-1.08;
xticklabelz11_2 = ["-1.10", "-1.09", "-1.08"];
xlimz12_2       = [0.885, 0.910];
xtickz12_2      = 0.89:0.01:0.91;
xticklabelz12_2 = ["0.89", "0.90", "0.91"];
xlimz11_3       = xlimz11_1;
xtickz11_3      = xtickz11_1;
xticklabelz11_3 = xticklabelz11_1;
xlimz12_3       = xlimz12_1;
xtickz12_3      = xtickz12_1;
xticklabelz12_3 = xticklabelz12_1;
xlimz11_4       = [-0.910, -0.885];
xtickz11_4      =-0.91:0.01:-0.89;
xticklabelz11_4 = ["-0.91", "-0.90", "-0.89"];
xlimz12_4       = [1.080, 1.105];
xtickz12_4      = 1.08:0.01:1.10;
xticklabelz12_4 = ["1.08", "1.09", "1.10"];
%--------------------------------------------------------------------------
ylimz11_1       = [0.0980, 0.1020];
ytickz11_1      = 0.0980:0.0010:0.1020;
yticklabelz11_1 = ["0.0980", "0.0990", "0.1000", "0.1010", "0.1020"];
ylimz12_1       = [0.1000, 0.1015];
ytickz12_1      = 0.1000:0.0005:0.1015;
yticklabelz12_1 = ["0.1000", "0.1005", "0.1010", "0.1015"];
ylimz11_2       = [0.1195, 0.1220];
ytickz11_2      = 0.1195:0.0005:0.1220;
yticklabelz11_2 = ["0.1195", "0.1200", "0.1205", "0.1210", "0.1215", "0.1220"];
ylimz12_2       = [0.0810, 0.0825];
ytickz12_2      = 0.0810:0.0005:0.0825;
yticklabelz12_2 = ["0.0810", "0.0815", "0.0820", "0.0825"];
ylimz11_3       = ylimz12_1;
ytickz11_3      = ytickz12_1;
yticklabelz11_3 = yticklabelz12_1;
ylimz12_3       = ylimz11_1;
ytickz12_3      = ytickz11_1;
yticklabelz12_3 = yticklabelz11_1;
ylimz11_4       = ylimz12_2;
ytickz11_4      = ytickz12_2;
yticklabelz11_4 = yticklabelz12_2;
ylimz12_4       = ylimz11_2;
ytickz12_4      = ytickz11_2;
yticklabelz12_4 = yticklabelz11_2;
%--------------------------------------------------------------------------
zm              = 0;
zf              = 1.025./2;
%--------------------------------------------------------------------------
xz11_1          = (xlimz11_1-zm).*zf+zm;
xz12_1          = (xlimz12_1-zm).*zf+zm;
xz11_2          = (xlimz11_2-zm).*zf+zm;
xz12_2          = (xlimz12_2-zm).*zf+zm;
xz11_3          = (xlimz11_3-zm).*zf+zm;
xz12_3          = (xlimz12_3-zm).*zf+zm;
xz11_4          = (xlimz11_4-zm).*zf+zm;
xz12_4          = (xlimz12_4-zm).*zf+zm;
yz11_1          = ylimz11_1;
yz12_1          = ylimz12_1;
yz11_2          = ylimz11_2;
yz12_2          = ylimz12_2;
yz11_3          = ylimz11_3;
yz12_3          = ylimz12_3;
yz11_4          = ylimz11_4;
yz12_4          = ylimz12_4;
%--------------------------------------------------------------------------
xr              = [-0.5000, 0.5000].*zf+zm;
dx              = 1.0./3.0.*0.25;
yr1             = [ 0.1165   , 0.1165+dx];
yr2             = [ 0.0925-dx, 0.0925];
%--------------------------------------------------------------------------
xlimz1          = {xlimz11_1, xlimz11_2, xlimz11_3, xlimz11_4};
xlimz2          = {xlimz12_1, xlimz12_2, xlimz12_3, xlimz12_4};
ylimz1          = {ylimz11_1, ylimz11_2, ylimz11_3, ylimz11_4};
ylimz2          = {ylimz12_1, ylimz12_2, ylimz12_3, ylimz12_4};
xtickz1         = {xtickz11_1, xtickz11_2, xtickz11_3, xtickz11_4};
xtickz2         = {xtickz12_1, xtickz12_2, xtickz12_3, xtickz12_4};
ytickz1         = {ytickz11_1, ytickz11_2, ytickz11_3, ytickz11_4};
ytickz2         = {ytickz12_1, ytickz12_2, ytickz12_3, ytickz12_4};
xticklabelz1    = {xticklabelz11_1, xticklabelz11_2, xticklabelz11_3, xticklabelz11_4};
xticklabelz2    = {xticklabelz12_1, xticklabelz12_2, xticklabelz12_3, xticklabelz12_4};
yticklabelz1    = {yticklabelz11_1, yticklabelz11_2, yticklabelz11_3, yticklabelz11_4};
yticklabelz2    = {yticklabelz12_1, yticklabelz12_2, yticklabelz12_3, yticklabelz12_4};
%--------------------------------------------------------------------------
xz1             = {xz11_1, xz11_2, xz11_3, xz11_4};
xz2             = {xz12_1, xz12_2, xz12_3, xz12_4};
xr1             = {xr, xr, xr, xr};
xr2             = {xr, xr, xr, xr};
yz1             = {yz11_1, yz11_2, yz11_3, yz11_4};
yz2             = {yz12_1, yz12_2, yz12_3, yz12_4};
yr1             = {yr1, yr1, yr1, yr1};
yr2             = {yr2, yr2, yr2, yr2};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid1            = ["Ol (025T).pdf", "Ol (050T).pdf", "Ol (075T).pdf", "Ol (100T).pdf"];
%--------------------------------------------------------------------------
P1              = cell(m, 1+n);
P2              = cell(2, m);
P3              = cell(2, m, ns);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:n
        %------------------------------------------------------------------
        fig1 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        cp = get(gca, 'Position');
        hold on;
        %------------------------------------------------------------------
        P1{1, 1} = fplot(@(x) g.data.Z (x, tk(1, j)), M(1, 1), 'Color', Color(1, :), 'LineWidth', LW);
        for k = 1:n
            P1{1, 1+k} = plot(xc, Fh1(nk(1, 1+j), :, k, 1)', M(1, 1+k), 'Color', Color(1+k, :), 'LineWidth', LW);
        end
        xbathy = linspace(xv(1, 1), xv(K+1, 1), 10e3);
        ybathy = g.data.Zb(xbathy);
        patch(...
            [xv(1, 1), xbathy, xv(K+1, 1)], [0, ybathy, 0], grey, 'EdgeColor', grey);
        %------------------------------------------------------------------
        PlotX1(...
            xlim, ylim{1, 1}, ...
            xtick, ytick{1, 1}, ...
            xticklabel, yticklabel{1, 1}, ...
            xtitle, ytitle{1, 1});
        PlotX2(...
            P1(1, :), legendlabel, 'northeast', numcols);
        %------------------------------------------------------------------
        MagInset(gcf, -1, ...
            [xz1{1, j}, yz1{1, j}], ...
            [xr1{1, j}, yr1{1, j}], {'NW', 'NW'; 'SW', 'SW'}, cp, xlimz1{1, j}, ylimz1{1, j}, xtickz1{1, j}, ytickz1{1, j}, xticklabelz1{1, j}, yticklabelz1{1, j}, 10.5);
        MagInset(gcf, -1, ...
            [xz2{1, j}, yz2{1, j}], ...
            [xr2{1, j}, yr2{1, j}], {'NE', 'NE'; 'SE', 'SE'}, cp, xlimz2{1, j}, ylimz2{1, j}, xtickz2{1, j}, ytickz2{1, j}, xticklabelz2{1, j}, yticklabelz2{1, j}, 10.5);
        %------------------------------------------------------------------
        if export(1, 1)
            exportgraphics(fig1, fid1(1, j), 'ContentType', 'Vector', 'Resolution', 600);
            movefile      (fid1(1, j), path3);
        end
        %------------------------------------------------------------------
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:2
        for j = 1:2
            %--------------------------------------------------------------
            if ~export(1, 2)
                fig2 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
            else
                fig2 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
            end
            hold on;
            %--------------------------------------------------------------
            for k = 1:ns
                P3{i, j, k} = plot(xc, Fh1(nk(1, k), :, p, j)', '-', 'Color', lblue(11, :), 'LineWidth', LW, 'Visible', 'off');
            end
            P2{i, j} = plot(xc, Fh1(1, :, p, j)', '-', 'Color', Color(2, :), 'LineWidth', LW);
            %--------------------------------------------------------------
            if j == 1
                xbathy = linspace(xv(1, 1), xv(K+1, 1), 10e3);
                ybathy = g.data.Zb(xbathy);
                patch(...
                    [xv(1, 1), xbathy, xv(K+1, 1)], [0, ybathy, 0], grey, 'EdgeColor', grey);
            end
            %--------------------------------------------------------------
            if i == 2
                set(P3{i, j, 1}, 'Visible', 'on');
            end
            %--------------------------------------------------------------
            PlotX1(...
                xlim, ylim{j, 1}, ...
                xtick, ytick{j, 1}, ...
                xticklabel, yticklabel{j, 1}, ...
                xtitle, ytitle{j, 1});
            if j == 2
                ax = gca;
                xp = ax.XLim(1);
                yp = ax.YLim(2);
                text(xp, yp, "$\times 10^{-3}$", ...
                    'FontSize', 19.5, 'Interpreter', 'latex', ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
            end
            %--------------------------------------------------------------
            video.path_mp4 = sprintf("%s/Ol%d_P%d (%s)", path2, i, p, str1(1, j));
            video.path_png = sprintf("%s/FRAMES/", path2);
            video.fig      = fig2;
            video.P1       = P2(i, j);
            video.P2       = P3(i, j, :);
            video.Fh1      = Fh1(:, :, p, j);
            video.nf       = nf;
            video.nk       = nk;
            %--------------------------------------------------------------
            switch i
                case 1
                    WriteVideo1(video, export(1, 2));
                case 2
                    WriteVideo2(video, export(1, 2));
                otherwise
                    return
            end
            %--------------------------------------------------------------
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end