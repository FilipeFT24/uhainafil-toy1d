function [] = PlotT04C(g, export)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color = linspecer(9, 'Qualitative');
grey  = repmat(0.80, 1, 3);
nc    = 50;
ns    = 3;
nf    = g.data.nf+numel(g.data.tk);
R1    = linspace(1, Color(1, 1), nc);
G1    = linspace(1, Color(1, 2), nc);
B1    = linspace(1, Color(1, 3), nc);
R2    = linspace(1, Color(2, 1), nc);
G2    = linspace(1, Color(2, 2), nc);
B2    = linspace(1, Color(2, 3), nc);
lred  = [R1', G1', B1'];
lblue = [R2', G2', B2'];
M     = [":", "-", "-.", "--"];
LW    = 2.5;
flag  = [0, 1];
path1 = "Plot/T04 - Reflection of shoaling waves/DATA";
path2 = "Plot/T04 - Reflection of shoaling waves";
path3 = "Plot/T04 - Reflection of shoaling waves/";
type  = 0;
switch type
    case 0
        path4 = "GN";
    case 1
        path4 = "GN_159";
    otherwise
        return
end
%--------------------------------------------------------------------------
h0    = g.data.h0;
xm    =-55;
xM    = 20;
K     = 750;
xv    = linspace(xm, xM, K+1)';
xc    = zeros(K, 1);
for i = 1:K
    xc(i, 1) = 1./2.*(xv(i, 1)+xv(i+1, 1));
end
A     = [0.1000; 0.1834].*h0;
c     = sqrt(g.data.G.*(h0+A));
%--------------------------------------------------------------------------
m     = 2;
n     = 3;
ng    = 3;
p     = 3;
Fh1   = zeros(nf, K   , n, m);
Fh2   = zeros(ns, 1+K , n, m);
Fh3   = zeros(nf, 1+ng, n, m);
Fe    = cell (ng, m);
th1   = zeros(nf, n, m);
str1  = ["1000", "1834"];
str2  = ["walkleya1.txt", "walkleya2.txt", "walkleya3.txt"; "walkleyb1.txt", "walkleyb2.txt", "walkleyb3.txt"];
for i = 1:n
    for j = 1:m
        Fh1(:, :, i, j) = getdata1(sprintf("%s/0.%s/Ref_P%d_%s_05_%s.mat", path1, str1(1, j), i, str1(1, j), path4), i, 1);
        Fh2(:, :, i, j) = getdata2(sprintf("%s/0.%s/Ref_P%d_%s_05_%s.mat", path1, str1(1, j), i, str1(1, j), path4), i, 1, ns);
        Fh3(:, :, i, j) = getdatag(sprintf("%s/0.%s/Ref_P%d_%s_05_%s.mat", path1, str1(1, j), i, str1(1, j), path4));
        th1(:, i, j)    = getdata3(sprintf("%s/0.%s/Ref_P%d_%s_05_%s.mat", path1, str1(1, j), i, str1(1, j), path4));
    end
end
for i = 1:m
    for j = 1:ng
        Fe{i, j}       = readmatrix(str2(i, j));
        Fe{i, j}(:, 2) = Fe{i, j}(:, 2)./h0;
    end
end
%--------------------------------------------------------------------------
Fh1                  = Fh1./h0;
Fh2(:, 2:K +1, :, :) = Fh2(:, 2:K +1, :, :)./h0;
Fh3(:, 2:ng+1, :, :) = Fh3(:, 2:ng+1, :, :)./h0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim1           = [0, 30];
xtick1          = 0:5:30;
xticklabel1     = ["0", "5", "10", "15", "20", "25", "30"];
xlim2           = [-55, 21.875];
xtick2          = [-55, -50:10:20];
xticklabel2     = ["-55", "-50", "-40", "-30", "-20", "-10", "0", "10", "20"];
%--------------------------------------------------------------------------
ylim11          = [-0.02, 0.20];
ytick11         = [-0.02, 0.00:0.05:0.20];
yticklabel11    = ["-0.02", "0.00", "0.05", "0.10", "0.15", "0.20"];
ylim12          = [-0.05, 0.40];
ytick12         = [-0.05, 0.00:0.10:0.40];
yticklabel12    = ["-0.05", "0.00", "0.10", "0.20", "0.30", "0.40"];
ylim21          = [-1.00, 0.70];
ytick21         = [-1.00:0.20:0.60, 0.70];
yticklabel21    = ["-1.00", "-0.80", "-0.60", "-0.40", "-0.20", "0.00", "0.20", "0.40", "0.60", "0.70"];
ylim22          = [-0.05, 0.70];
ytick22         = [-0.05, 0.00:0.10:0.70];
yticklabel22    = ["-0.05", "0.00", "0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70"];
%--------------------------------------------------------------------------
ylim1           = {ylim11; ylim12};
ytick1          = {ytick11; ytick12};
yticklabel1     = {yticklabel11; yticklabel12};
%--------------------------------------------------------------------------
legendlabel11   = ["$g_{0}$", "", "", "$k=1$", "$k=2$", "$k=3$"];
legendlabel12   = ["$g_{1}$", "", "", "$k=1$", "$k=2$", "$k=3$"];
legendlabel13   = ["$g_{2}$", "", "", "$k=1$", "$k=2$", "$k=3$"];
%--------------------------------------------------------------------------
legendlabel1    = {legendlabel11; legendlabel12; legendlabel13};
%--------------------------------------------------------------------------
legendlabel2    = ["$\epsilon = 0.07$", "$\epsilon = 0.12$"];
numcols1        = 2;
numcols2        = 1;
xtitle1         = "$t$";
xtitle2         = "$x$";
ytitle1         = "$\zeta^{\prime}$";
ytitle2         = "$z^{\prime}$";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlimz11_1       = [10.0, 12.0];
xtickz11_1      = 10.0:0.5:12.0;
xticklabelz11_1 = ["10.0", "10.5", "11.0", "11.5", "12.0"];
xlimz12_1       = [27.0, 29.0];
xtickz12_1      = 27.0:0.5:29.0;
xticklabelz12_1 = ["27.0", "27.5", "28.0", "28.5", "29.0"];
xlimz11_2       = [16, 22];
xtickz11_2      = 16:1:22;
xticklabelz11_2 = ["16", "17", "18", "19", "20", "21", "22"];
xlimz11_3       = [16, 22];
xtickz11_3      = 16:1:22;
xticklabelz11_3 = ["16", "17", "18", "19", "20", "21", "22"];
xlimz11_4       = [9.5, 11.5];
xtickz11_4      = 9.5:0.5:11.5;
xticklabelz11_4 = ["9.5", "10.0", "10.5", "11.0", "11.5"];
xlimz12_4       = [26.0, 28.0];
xtickz12_4      = 26.0:0.5:28.0;
xticklabelz12_4 = ["26.0", "26.5", "27.0", "27.5", "28.0"];
xlimz11_5       = [16, 21];
xtickz11_5      = 16:1:21;
xticklabelz11_5 = ["16", "17", "18", "19", "20", "21"];
xlimz11_6       = [16, 21];
xtickz11_6      = 16:1:21;
xticklabelz11_6 = ["16", "17", "18", "19", "20", "21"];
%--------------------------------------------------------------------------
ylimz11_1       = [0.085, 0.110];
ytickz11_1      = 0.085:0.005:0.110;
yticklabelz11_1 = ["0.085", "0.090", "0.095", "0.100", "0.105", "0.110"];
ylimz12_1       = [0.000, 0.125];
ytickz12_1      = 0.000:0.025:0.125;
yticklabelz12_1 = ["0.000", "0.025", "0.050", "0.075", "0.100", "0.125"];
ylimz11_2       = [0.08, 0.15];
ytickz11_2      = 0.08:0.01:0.15;
yticklabelz11_2 = ["0.08", "0.09", "0.10", "0.11", "0.12", "0.13", "0.14", "0.15"];
ylimz11_3       = [0.09, 0.16];
ytickz11_3      = 0.09:0.01:0.16;
yticklabelz11_3 = ["0.09", "0.10", "0.11", "0.12", "0.13", "0.14", "0.15", "0.16"];
ylimz11_4       = [0.15, 0.20];
ytickz11_4      = 0.15:0.01:0.20;
yticklabelz11_4 = ["0.15", "0.16", "0.17", "0.18", "0.19", "0.20"];
ylimz12_4       = [-0.05, 0.25];
ytickz12_4      =-0.05:0.05:0.25;
yticklabelz12_4 = ["-0.05", "0.00", "0.05", "0.10", "0.15", "0.20", "0.25"];
ylimz11_5       = [0.150, 0.300];
ytickz11_5      = 0.150:0.025:0.300;
yticklabelz11_5 = ["0.150", "0.175", "0.200", "0.225", "0.250", "0.275", "0.300"];
ylimz11_6       = [0.175, 0.325];
ytickz11_6      = 0.175:0.025:0325;
yticklabelz11_6 = ["0.175", "0.200", "0.225", "0.250", "0.275", "0.300", "0.325"];
%--------------------------------------------------------------------------
zm              = 15;
zf              = 1.025./2;
%--------------------------------------------------------------------------
xz11_1          = (xlimz11_1-zm).*zf+zm;
xz12_1          = (xlimz12_1-zm).*zf+zm;
xz11_2          = (xlimz11_2-zm).*zf+zm;
xz11_3          = (xlimz11_3-zm).*zf+zm;
xz11_4          = (xlimz11_4-zm).*zf+zm;
xz12_4          = (xlimz12_4-zm).*zf+zm;
xz11_5          = (xlimz11_5-zm).*zf+zm;
xz11_6          = (xlimz11_6-zm).*zf+zm;
yz11_1          = ylimz11_1;
yz12_1          = ylimz12_1;
yz11_2          = ylimz11_2;
yz11_3          = ylimz11_3;
yz11_4          = ylimz11_4;
yz12_4          = ylimz12_4;
yz11_5          = ylimz11_5;
yz11_6          = ylimz11_6;
%--------------------------------------------------------------------------
xr11_1          = [ -9.00,  1.00].*zf+zm;
xr12_1          = [  0.25, 10.25].*zf+zm;
xr11_2          = [-12.50,  0.00].*zf+zm;
xr11_3          = [-12.50,  0.00].*zf+zm;
xr11_4          = [ -9.50,  0.50].*zf+zm;
xr12_4          = [ -0.75,  9.25].*zf+zm;
xr11_5          = [-12.50,  0.00].*zf+zm;
xr11_6          = [-12.50,  0.00].*zf+zm;
yr11_1          = [0.120, 0.120+(10.0./30.0.*0.22)];
yr12_1          = [0.015, 0.015+(10.0./30.0.*0.22)];
yr11_2          = [-1, 1].*(6.25./30.0.*0.22)+0.115;
yr11_3          = [-1, 1].*(6.25./30.0.*0.22)+0.125;
yr11_4          = [0.225, 0.225+(10.0./30.0.*0.45)];
yr12_4          = [0.035, 0.035+(10.0./30.0.*0.45)];
yr11_5          = [-1, 1].*(6.25./30.0.*0.45)+0.225;
yr11_6          = [-1, 1].*(6.25./30.0.*0.45)+0.250;
%--------------------------------------------------------------------------
xlimz1          = {xlimz11_1, xlimz11_2, xlimz11_3; xlimz11_4, xlimz11_5, xlimz11_6};
xlimz2          = {xlimz12_1; xlimz12_4};
ylimz1          = {ylimz11_1, ylimz11_2, ylimz11_3; ylimz11_4, ylimz11_5, ylimz11_6};
ylimz2          = {ylimz12_1; ylimz12_4};
xtickz1         = {xtickz11_1, xtickz11_2, xtickz11_3; xtickz11_4, xtickz11_5, xtickz11_6};
xtickz2         = {xtickz12_1; xtickz12_4};
ytickz1         = {ytickz11_1, ytickz11_2, ytickz11_3; ytickz11_4, ytickz11_5, ytickz11_6};
ytickz2         = {ytickz12_1; ytickz12_4};
xticklabelz1    = {xticklabelz11_1, xticklabelz11_2, xticklabelz11_3; xticklabelz11_4, xticklabelz11_5, xticklabelz11_6};
xticklabelz2    = {xticklabelz12_1; xticklabelz12_4};
yticklabelz1    = {yticklabelz11_1, yticklabelz11_2, yticklabelz11_3; yticklabelz11_4, yticklabelz11_5, yticklabelz11_6};
yticklabelz2    = {yticklabelz12_1; yticklabelz12_4};
%--------------------------------------------------------------------------
xz1             = {xz11_1, xz11_2, xz11_3; xz11_4, xz11_5, xz11_6};
xz2             = {xz12_1; xz12_4};
xr1             = {xr11_1, xr11_2, xr11_3; xr11_4, xr11_5, xr11_6};
xr2             = {xr12_1; xr12_4};
yz1             = {yz11_1, yz11_2, yz11_3; yz11_4, yz11_5, yz11_6};
yz2             = {yz12_1; yz12_4};
yr1             = {yr11_1, yr11_2, yr11_3; yr11_4, yr11_5, yr11_6};
yr2             = {yr12_1; yr12_4};
%--------------------------------------------------------------------------
fid01           = ["Ref2 (t0).pdf", "Ref2 (tf).pdf"];
fid02           = ["Ref4 (t0).pdf", "Ref4 (tf).pdf"];
switch type
    case 0
        fid1            = [...
            "Ref0 (gn, 0.07).pdf", "Ref1 (gn, 0.07).pdf", "Ref2 (gn, 0.07).pdf";
            "Ref0 (gn, 0.12).pdf", "Ref1 (gn, 0.12).pdf", "Ref2 (gn, 0.12).pdf"];
    case 1
        fid1            = [...
            "Ref0 (e-gn, 0.07).pdf", "Ref1 (e-gn, 0.07).pdf", "Ref2 (e-gn, 0.07).pdf";
            "Ref0 (e-gn, 0.12).pdf", "Ref1 (e-gn, 0.12).pdf", "Ref2 (e-gn, 0.12).pdf"];
    otherwise
        return
end
%--------------------------------------------------------------------------
P1              = cell(2, m, 1+2+n);
P21_1           = cell(1, m);
P21_2           = cell(1, m);
P22_1           = cell(m, ns);
P22_2           = cell(m, ns);
P23_1           = cell(m, ns);
P23_2           = cell(m, ns);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:m
        for j = 1:ng
            %--------------------------------------------------------------
            fig1 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
            cp = get(gca, 'Position');
            hold on;
            %--------------------------------------------------------------
            xline(50./c(i, 1), '--', ...
                'Color', 'k');
            for k = 1:n-1
                xline(Fh2(3, 1, k, i), '-', ...
                    'Color', 'k', ...
                    'LineWidth', 0.5);
            end
            xline(Fh2(3, 1, n, i), '-', "$t_{r}$", ...
                'Color', 'k', ...
                'LineWidth', 0.5, ...
                'FontSize', 22.5, ...
                'LabelVerticalAlignment', 'bottom', ...
                'LabelOrientation', 'horizontal');
            %--------------------------------------------------------------
            P1{i, j, 1} = plot(Fe{i, j}(:, 1), Fe{i, j}(:, 2), M(1, 1), 'Color', Color(1, :), 'LineWidth', LW);
            for k = 2:3
                P1{i, j, k} = plot(nan, nan, M(1, 1), 'Color', 'w', 'LineWidth', LW);
            end
            for k = 1:n
                P1{i, j, 3+k} = plot(Fh3(:, 1, k, i), Fh3(:, 1+j, k, i), M(1, 1+k), 'Color', Color(1+k, :), 'LineWidth', LW);
            end
            %--------------------------------------------------------------
            PlotX1(...
                xlim1, ylim1{i, 1}, ...
                xtick1, ytick1{i, 1}, ...
                xticklabel1, yticklabel1{i, 1}, ...
                xtitle1, ytitle1);
            PlotX2(...
                P1(i, j, :), legendlabel1{j, 1}, 'northeast', numcols1);
            %--------------------------------------------------------------
            if j == 1
                MagInset(gcf, -1, ...
                    [xz1{i, j}, yz1{i, j}], ...
                    [xr1{i, j}, yr1{i, j}], {'SE', 'SE'; 'SW', 'SW'}, cp, xlimz1{i, j}, ylimz1{i, j}, xtickz1{i, j}, ytickz1{i, j}, xticklabelz1{i, j}, yticklabelz1{i, j}, 10.5);
                MagInset(gcf, -1, ...
                    [xz2{i, 1}, yz2{i, 1}], ...
                    [xr2{i, 1}, yr2{i, 1}], {'NE', 'NE'; 'NW', 'NW'}, cp, xlimz2{i, j}, ylimz2{i, j}, xtickz2{i, j}, ytickz2{i, j}, xticklabelz2{i, j}, yticklabelz2{i, j}, 10.5);
            else
                MagInset(gcf, -1, ...
                    [xz1{i, j}, yz1{i, j}], ...
                    [xr1{i, j}, yr1{i, j}], {'NE', 'NE'; 'SE', 'SE'}, cp, xlimz1{i, j}, ylimz1{i, j}, xtickz1{i, j}, ytickz1{i, j}, xticklabelz1{i, j}, yticklabelz1{i, j}, 10.5);
            end
            %--------------------------------------------------------------
            if export(1, 1)
                exportgraphics(fig1, fid1(i, j), 'ContentType', 'Vector', 'Resolution', 600);
                movefile      (fid1(i, j), path3);
            end
            %--------------------------------------------------------------
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PB = cell(1, ns);
    TB = cell(1, ns);
    for i = 1:2
        if ~export(1, 2)
            fig21 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        else
            fig21 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
        end
        hold on;
        %------------------------------------------------------------------
        for j = 1:ns
            P22_1{1, j} = plot(xc, Fh2(j, 2:K+1, p, 1)', '-', 'Color', lred (11, :), 'LineWidth', LW, 'Visible', 'off');
            P22_1{2, j} = plot(xc, Fh2(j, 2:K+1, p, 2)', '-', 'Color', lblue(11, :), 'LineWidth', LW, 'Visible', 'off');
        end
        for j = 1:m
            P21_1{1, j} = plot(xc, Fh1(1, :, p, j)', '-', 'Color', Color(j, :), 'LineWidth', LW);
            P23_1{1, j} = plot(nan, nan, '-', 'Color', Color(j, :), 'LineWidth', LW);
        end
        patch(...
            [-55, -55, 0, 20, 20, 21.875, 21.875], [-0.70, -0.70, -0.70, -0.30, 0.70, 0.70, -0.70]./h0, grey, 'EdgeColor', grey);
        if i == 2
            for j = 1:m
                set(P22_1{j, 1}, 'Visible', 'on');
            end
        end
        %------------------------------------------------------------------
        PlotX1(...
            xlim2, ylim21, ...
            xtick2, ytick21, ...
            xticklabel2, yticklabel21, ...
            xtitle2, ytitle2);
        PlotX2(...
            P23_1, legendlabel2, 'northwest', numcols2);
        %------------------------------------------------------------------
        if i == 2 && export(1, 1)
            xg = [0.00, 16.25, 17.75];
            yg = [0.00,  0.00,  0.00];
            PA = plot(xg, yg, 'o', 'Color', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 5);
            for j = 1:ns
                PB{1, j} = plot(...
                    [xg(1, j), xg(1, j)], [0.00, -1.00], '--', 'Color', 'k', 'LineWidth', 0.01);
                TB{1, j} = text(xg(1, j), yg(1, j)+0.05, sprintf("$g_{%d}$", j-1), 'FontSize', 12.5, 'HorizontalAlignment', 'center');
            end
            for j = 1:m
                for k = 2:ns
                    set(P22_1{j, k}, 'Visible', 'on');
                end
            end
            exportgraphics(fig21, fid01(1, 1), 'ContentType', 'Vector', 'Resolution', 600);
            movefile      (fid01(1, 1), path3);
            for j = 1:m
                for k = 2:ns
                    set(P22_1{j, k}, 'Visible', 'off');
                end
            end
            set(PA, 'Visible', 'off');
            for j = 1:ns
                set(PB{1, j}, 'Visible', 'off');
                set(TB{1, j}, 'Visible', 'off');
            end
        end
        %------------------------------------------------------------------
        video.path_mp4 = sprintf("%s/Ref%d_P%d", path2, i, p);
        video.path_png = sprintf("%s/FRAMES/", path2);
        video.fig      = fig21;
        video.P1       = P21_1;
        video.P2       = P22_1;
        video.Fh1      = permute(Fh1(:, :, p, :), [1, 2, 4, 3]);
        video.nf       = nf;
        video.nk       = [sum(Fh2(:, 1, p, 1)' > th1(:, p, 1), 1); sum(Fh2(:, 1, p, 2)' > th1(:, p, 2), 1)];
        video.nk(:, 2) = nf+1;
        %------------------------------------------------------------------
        switch i
            case 1
                WriteVideo1(video, export(1, 2));
            case 2
                WriteVideo2(video, export(1, 2));
            otherwise
                return
        end
        %------------------------------------------------------------------
        if i == 2 && export(1, 1)
            set(PA, 'Visible', 'off');
            for j = 1:ns
                set(PB{1, j}, 'Visible', 'off');
                set(TB{1, j}, 'Visible', 'off');
            end
            exportgraphics(fig21, fid01(1, 2), 'ContentType', 'Vector', 'Resolution', 600);
            movefile      (fid01(1, 2), path3);
        end
        %------------------------------------------------------------------
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:2
        if ~export(1, 2)
            fig22 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        else
            fig22 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
        end
        hold on;
        %------------------------------------------------------------------
        for j = 1:ns
            P22_2{1, j} = plot(xc, Fh2(j, 2:K+1, p, 1)', '-', 'Color', lred (11, :), 'LineWidth', LW, 'Visible', 'off');
            P22_2{2, j} = plot(xc, Fh2(j, 2:K+1, p, 2)', '-', 'Color', lblue(11, :), 'LineWidth', LW, 'Visible', 'off');
        end
        for j = 1:m
            P21_2{1, j} = plot(xc, Fh1(1, :, p, j)', '-', 'Color', Color(j, :), 'LineWidth', LW);
            P23_2{1, j} = plot(nan, nan, '-', 'Color', Color(j, :), 'LineWidth', LW);
        end
        patch(...
            [20, 20, 21.875, 21.875], [-0.05, 0.70, 0.70, -0.05], grey, 'EdgeColor', grey);
        %------------------------------------------------------------------
        if i == 2
            for j = 1:m
                set(P22_2{j, 1}, 'Visible', 'on');
            end
        end
        %------------------------------------------------------------------
        PlotX1(...
            xlim2, ylim22, ...
            xtick2, ytick22, ...
            xticklabel2, yticklabel22, ...
            xtitle2, ytitle1);
        PlotX2(...
            P23_2, legendlabel2, 'northwest', numcols2);
        %------------------------------------------------------------------
        if i == 2 && export(1, 1)
            xg = [0.00, 16.25, 17.75];
            yg = [0.00,  0.00,  0.00];
            PA = plot(xg, yg, 'o', 'Color', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 5);
            for j = 1:ns
                PB{1, j} = plot(...
                    [xg(1, j), xg(1, j)], [0.00, -0.05], '--', 'Color', 'k', 'LineWidth', 0.01);
                TB{1, j} = text(xg(1, j), yg(1, j)+0.02, sprintf("$g_{%d}$", j-1), 'FontSize', 12.5, 'HorizontalAlignment', 'center');
            end
            for j = 1:m
                for k = 2:ns
                    set(P22_2{j, k}, 'Visible', 'on');
                end
            end
            exportgraphics(fig22, fid02(1, 1), 'ContentType', 'Vector', 'Resolution', 600);
            movefile      (fid02(1, 1), path3);
            for j = 1:m
                for k = 2:ns
                    set(P22_2{j, k}, 'Visible', 'off');
                end
            end
            set(PA, 'Visible', 'off');
            for j = 1:ns
                set(PB{1, j}, 'Visible', 'off');
                set(TB{1, j}, 'Visible', 'off');
            end
        end
        %------------------------------------------------------------------
        video.path_mp4 = sprintf("%s/Ref%d_P%d", path2, 2+i, p);
        video.path_png = sprintf("%s/FRAMES/", path2);
        video.fig      = fig22;
        video.P1       = P21_2;
        video.P2       = P22_2;
        video.Fh1      = permute(Fh1(:, :, p, :), [1, 2, 4, 3]);
        video.nf       = nf;
        video.nk       = [sum(Fh2(:, 1, p, 1)' > th1(:, p, 1), 1); sum(Fh2(:, 1, p, 2)' > th1(:, p, 2), 1)];
        video.nk(:, 2) = nf+1;
        %------------------------------------------------------------------
        switch i
            case 1
                WriteVideo1(video, export(1, 2));
            case 2
                WriteVideo2(video, export(1, 2));
            otherwise
                return
        end
        %------------------------------------------------------------------
        if i == 2 && export(1, 1)
            exportgraphics(fig22, fid02(1, 2), 'ContentType', 'Vector', 'Resolution', 600);
            movefile      (fid02(1, 2), path3);
        end
        %------------------------------------------------------------------
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end