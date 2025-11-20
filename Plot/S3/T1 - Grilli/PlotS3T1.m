function [] = PlotS3T1(g, export)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color  = linspecer(9, 'Qualitative');
grey   = repmat(0.80, 1, 3);
nc     = 50;
ns     = 5;
nf     = g.data.nf+numel(g.data.tk);
R1     = linspace(1, Color(1, 1), nc);
G1     = linspace(1, Color(1, 2), nc);
B1     = linspace(1, Color(1, 3), nc);
R2     = linspace(1, Color(2, 1), nc);
G2     = linspace(1, Color(2, 2), nc);
B2     = linspace(1, Color(2, 3), nc);
R3     = linspace(1, Color(3, 1), nc);
G3     = linspace(1, Color(3, 2), nc);
B3     = linspace(1, Color(3, 3), nc);
lred   = [R1', G1', B1'];
lblue  = [R2', G2', B2'];
lgreen = [R3', G3', B3'];
M      = ["o", "-", "-.", "--"];
lw     = 1.5;
LW     = 2.5;
flag   = [0, 0, 1, 0];
path1  = "Plot/S3/T1 - Grilli/DATA";
path2  = "Plot/S3/T1 - Grilli";
path3  = "Plot/S3/T1 - Grilli/";
bathy  = 0;
type   = 1;
switch type
    case 0
        path4 = "GN";
    case 1
        path4 = "GN_159";
    otherwise
        return
end
path5  = sprintf("C%d", bathy);
%--------------------------------------------------------------------------
h0     = 0.44;
slope  = 1./34.70;
xm     =-52.32.*h0;
xM     = h0./slope;
K      = 383;
xv     = linspace(xm, xM, K+1)';
xc     = zeros(K, 1);
for i = 1:K
    xc(i, 1) = 1./2.*(xv(i, 1)+xv(i+1, 1));
end
%--------------------------------------------------------------------------
ds     = 0.5./h0;
f1     = bathyC1(slope, ds, 1);
f2     = bathyC2(slope, ds, 1);
f3     = bathyC3(slope, ds, 1);
ZbC0   = @(x) -1+heaviside(x).*slope.*x;
ZbC1   = @(x) -1+...
    (heaviside(x+ds)-heaviside(x-ds)).*f1(x)+...
    (heaviside(x-ds)).*slope.*x;
ZbC2   = @(x) -1+...
    (heaviside(x+ds)-heaviside(x-ds)).*f2(x)+...
    (heaviside(x-ds)).*slope.*x;
ZbC3   = @(x) -1+...
    (heaviside(x+ds)-heaviside(x-ds)).*f3(x)+...
    (heaviside(x-ds)).*slope.*x;
%--------------------------------------------------------------------------
n     = 3;
ne    = 3;
ng    = 6;
Fh1   = zeros(nf, K, n, ne);
Fh2   = cell (ne, n);
Fh3   = zeros(nf, 1+ng, n, ne);
Fe1   = cell (ng, 1);
Fe2   = cell (ne, 1);
Xc1   = cell (ne, 1);
th1   = zeros(nf, n);
str1  = ["10", "15", "20", "25"];
str2  = ["Exp. data (g0).txt", "Exp. data (g1).txt", "Exp. data (g3).txt", "Exp. data (g5).txt", "Exp. data (g7).txt", "Exp. data (g9).txt"];
str3  = ["Exp. data (a).txt", "Exp. data (b).txt", "Exp. data (c).txt", "Exp. data (d).txt"];
lim2  = [...
    30.0, ... % (a)
    28.0, ... % (b)
    26.0, ... % (c)
    24.5];    % (d)
for i = 1:n
    for j = 3:ne
        %if j ~= 3
            Fh1(:, :, i, j) = getdata1(sprintf("%s/0.%s/C0/Grilli_P%d_383.mat", path1, str1(1, j), i), i, 1);
            Fh3(:, :, i, j) = getdatag(sprintf("%s/0.%s/C0/Grilli_P%d_383.mat", path1, str1(1, j), i));
            th1(:, i, j)    = getdata3(sprintf("%s/0.%s/C0/Grilli_P%d_383.mat", path1, str1(1, j), i));
        %else
        %    Fh1(:, :, i, j) = getdata1(sprintf("%s/0.%s/%s/Grilli_P%d_383.mat", path1, str1(1, j), path5, i), i, 1);
        %    Fh3(:, :, i, j) = getdatag(sprintf("%s/0.%s/%s/Grilli_P%d_383.mat", path1, str1(1, j), path5, i));
        %    th1(:, i, j)    = getdata3(sprintf("%s/0.%s/%s/Grilli_P%d_383.mat", path1, str1(1, j), path5, i));
        %end
    end
end
Fh1                  = Fh1-h0;
Fh3(:, 2:ng+1, :, :) = Fh3(:, 2:1+ng, :, :)-h0;
Fh1                  = Fh1./h0;
Fh3(:, 2:ng+1, :, :) = Fh3(:, 2:1+ng, :, :)./h0;
Fh3(:, 1     , :, :) = Fh3(:, 1     , :, :)./(sqrt(h0./g.data.G));
for i = 1:ng
    Fe1{i, 1} = readmatrix(str2(1, i));
end
for i = 1:ne
    Fe2{i, 1} = readmatrix(str3(1, i));
end
for j = 1:ne
    log       = xc./h0 > 0 & xc./h0 < lim2(1, j);
    Xc1{j, 1} = xc(log)./h0;
    switch bathy
        case 0
            zb = ZbC0(xc(log, 1)./h0);
        case 3
            zb = ZbC3(xc(log, 1)./h0);
        otherwise
            return
    end
    for i = 1:n
        aux       = max(Fh1(:, :, i, j), [], 1)';
        Fh2{j, i} =-aux(log, :)./zb;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim10          = [xm, xM]./h0;
xtick10         =-50:10:30;
xticklabel10    = ["-50", "-40", "-30", "-20", "-10", "0", "10", "20", "30"];
xlim11          = [5, 45];
xtick11         = 5:5:45;
xticklabel11    = ["5", "10", "15", "20", "25", "30", "35", "40", "45"];
xlim12          = [37, 47];
xtick12         = 37:47;
xticklabel12    = ["37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47"];
xlim2           = [19, 30];
xtick2          = 19:1:30;
xticklabel2     = ["19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30"];
%--------------------------------------------------------------------------
ylim10          = [-1.05, 0.00];
ytick10         = [-1.05, -1.00:0.20:0.00];
yticklabel10    = ["-1.05", "-1.00", "-0.80", "-0.60", "-0.40", "-0.20", "0.00"];
ylim11          = [-0.02, 0.40];
ytick11         = [-0.02, 0.00:0.05:0.40];
yticklabel11    = ["-0.02", "0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40"];
ylim2           = [0.2, 2.0];
ytick2          = 0.2:0.2:2.0;
yticklabel2     = ["0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0"];
ylim3           = [-1.00, 0.40];
ytick3          =-1.00:0.20:0.40;
yticklabel3     = ["-1.00", "-0.80", "-0.60", "-0.40", "-0.20", "0.00", "0.20", "0.40"];
ylim4           = [-0.15, 0.15];
ytick4          =-0.15:0.05:0.15;
yticklabel4     = ["-0.15", "-0.10", "-0.05", "0.00", "0.05", "0.10", "0.15"];
%--------------------------------------------------------------------------
legendlabel10   = ["$b$", "$\tilde{b}_{C^{1}}$", "$\tilde{b}_{C^{2}}$", "$\tilde{b}_{C^{3}}$"];
legendlabel11   = ["$g_{0}$", "", "", "$k=1$", "$k=2$", "$k=3$"];
legendlabel12   = ["$g_{1}$", "$g_{3}$", "$g_{5}$", "$g_{7}$", "$g_{9}$", "$k=1$", "$k=2$", "$k=3$"];
legendlabel2    = ["$(a)$", "$(b)$", "$(c)$", "$(d)$", "$k=1$", "$k=2$", "$k=3$"];
legendlabel3    = ["$k = 1$", "$k = 2$", "$k = 3$"];
numcols1        = 1;
numcols2        = 2;
xtitle1         = "$x^{\prime}$";
ytitle10        = "$z^{\prime}$";
ytitle11        = "$\zeta^{\prime}$";
ytitle2         = "$\zeta^{\prime}/\left|b\right|$";
ytitle4         = "$\phi$";
%--------------------------------------------------------------------------
xlimz10         = [-2.0, 2.0];
xtickz10        =-2.0:0.5:2.0;
xticklabelz10   = ["-2.0", "-1.5", "-1.0", "-0.5", "0.0", "0.5", "1.0", "1.5", "2.0"];
xlimz11_1       = [12, 16];
xtickz11_1      = 12:1:16;
xticklabelz11_1 = ["12", "13", "14", "15", "16"];
xlimz11_2       = [24, 36];
xtickz11_2      = 24:2:36;
xticklabelz11_2 = ["24", "26", "28", "30", "32", "34", "36"];
xlimz4          = [-5, 5];
xtickz4         =-5.0:2.5:5.0;
xticklabelz4    = ["-5.0", "-2.5", "0.0", "2.5", "5.0"];
%--------------------------------------------------------------------------
ylimz10         = [-1.005, -0.950];
ytickz10        = [-1.005, -1.000:0.010:-0.950];
yticklabelz10   = ["-1.005", "-1.000", "-0.990", "-0.980", "-0.970", "-0.960", "-0.950"];
ylimz11_1       = [0.180, 0.210];
ytickz11_1      = 0.180:0.005:0.210;
yticklabelz11_1 = ["0.180", "0.185", "0.190", "0.195", "0.200", "0.205", "0.210"];
ylimz11_2       = [-10, 10].*1e-3;
ytickz11_2      = (-10:5:10).*1e-3;
yticklabelz11_2 = ["-10", "-5", "0", "5", "10"];
ylimz3          = [-1.00, -0.95];
ytickz3         =-1.00:0.01:-0.95;
yticklabelz3    = ["-1.00", "-0.99", "-0.98", "-0.97", "-0.96", "-0.95"];
ylimz4          = [-5.0, 5.0].*1e-02;
ytickz4         = (-5.0:2.5:5.0).*1e-02;
yticklabelz4    = ["-5.0", "-2.5", "0.0", "2.5", "5.0"];
%--------------------------------------------------------------------------
L               = xM-xm;
zm1             = 1./2.*(xm+xM)./h0;
zm2             = 25;
zf              = 1.025./2;
%--------------------------------------------------------------------------
xz10            = (xlimz10  -zm1).*zf+zm1;
xz11_1          = (xlimz11_1-zm2).*zf+zm2;
xz11_2          = (xlimz11_2-zm2).*zf+zm2;
xz4             = (xlimz4   -zm1).*zf+zm1;
yz10            = ylimz10;
yz11_1          = ylimz11_1;
yz11_2          = ylimz11_2;
yz3             = ylimz3;
yz4             = ylimz4;
%--------------------------------------------------------------------------
xr10            = ([-40.0,  0.0]-zm1).*zf+zm1;
xr11_1          = ([ -2.5, 12.5]-0  ).*zf+zm2;
xr11_2          = ([ -5.0, 15.0]-0  ).*zf+zm2;
xr4             = ([-45.0, -5.0]-zm1).*zf+zm1;
dx_0            = (40.*h0./L.*1.05)./2;
dx_1            = (15./40.0 .*0.42)./2;
dx_2            = (20./40.0 .*0.42)./2;
dx_3            = (40.*h0./L.*1.40)./2;
yr10            = [-1, 1].*dx_0   -0.525;
yr11_1          = [-1, 1].*dx_1   +0.195;
yr11_2          = [-1, 1].*dx_2./5+0.060;
yr3             = [-1, 1].*dx_3   -0.500;
yr4             = [-0.125, -0.025];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P10             = cell(1, 4);
P11             = cell(1, 1+2+n);
P12             = cell(1, ng-1+n);
P2              = cell(1, n);
P31             = cell(1, n);
P32             = cell(1, n);
P33             = cell(n, ns);
P4              = cell(1, 4);
%--------------------------------------------------------------------------
fid10           = "bathymetry (0.5).pdf";
switch type
    case 0
        fid11 = sprintf("Grilli0 (%s, gn).pdf", path5);
        fid12 = sprintf("Grillik (%s, gn).pdf", path5);
        fid13 = sprintf("Grilli_env (%s, gn).pdf", path5);
    case 1
        fid11 = sprintf("Grilli0 (%s, e-gn).pdf", path5);
        fid12 = sprintf("Grillik (%s, e-gn).pdf", path5);
        fid13 = sprintf("Grilli_env (%s, e-gn).pdf", path5);
    otherwise
        return
end
fid14 = ["Grilli_phi (t1).pdf", "Grilli_phi (t2).pdf"];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fig10 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
%     cp = get(gca, 'Position');
%     hold on;
%     %----------------------------------------------------------------------
%     patch(...
%         [xm./h0, xm./h0, 0.00, xM./h0, xM./h0], [-1.05, -1.00, -1.00, 0.0, -1.05], grey, 'EdgeColor', grey);
%     plot(xv./h0, -ones(1, K+1), ...
%         '-|', 'Color', 'k', 'LineWidth', 0.01);
%     P10{1, 2} = fplot(@(x) ZbC1(x), 'Color', Color(2, :), 'LineWidth', LW);
%     P10{1, 3} = fplot(@(x) ZbC2(x), 'Color', Color(3, :), 'LineWidth', LW);
%     P10{1, 4} = fplot(@(x) ZbC3(x), 'Color', Color(4, :), 'LineWidth', LW);
%     P10{1, 1} = fplot(@(x) ZbC0(x), 'Color', Color(1, :), 'LineWidth', lw);
%     %----------------------------------------------------------------------
%     PlotX1(...
%         xlim10, ylim10, ...
%         xtick10, ytick10, ...
%         xticklabel10, yticklabel10, ...
%         xtitle1, ytitle10);
%     PlotX2(...
%         P10, legendlabel10, 'northeast', numcols1);
%     %----------------------------------------------------------------------
%     MagInset(gcf, -1, ...
%         [xz10, yz10], ...
%         [xr10, yr10], {'SE', 'SE'; 'SW', 'SW'}, cp, xlimz10, ylimz10, xtickz10, ytickz10, xticklabelz10, yticklabelz10, 10.5);
%     %------------------------------------------------------------------
%     if export(1, 1)
%         exportgraphics(fig10, fid10, 'ContentType', 'Vector', 'Resolution', 600);
%         movefile      (fid10, path3);
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig11 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
    cp = get(gca, 'Position');
    hold on;
    %----------------------------------------------------------------------
    P11{1, 1} = plot(Fe1{1, 1}(:, 1), Fe1{1, 1}(:, 2), M(1, 1), 'Color', Color(1, :), 'LineWidth', lw, 'MarkerFaceColor', 'w', 'MarkerSize', 10);
    for i = 2:3
        P11{1, i} = plot(nan, nan, M(1, 1), 'Color', 'w', 'LineWidth', LW);
    end
    for i = 1:n
        P11{1, 3+i} = plot(Fh3(:, 1, i, 3), Fh3(:, 2, i, 3), M(1, 1+i), 'Color', Color(1, :), 'LineWidth', LW);
    end
    %----------------------------------------------------------------------
    PlotX1(...
        xlim11, ylim11, ...
        xtick11, ytick11, ...
        xticklabel11, yticklabel11, ...
        xtitle1, ytitle11);
    PlotX2(...
        P11, legendlabel11, 'northeast', numcols2);
    %----------------------------------------------------------------------
    MagInset(gcf, -1, ...
        [xz11_1, yz11_1], ...
        [xr11_1, yr11_1], {'NW', 'NW'; 'SW', 'SW'}, cp, xlimz11_1, ylimz11_1, xtickz11_1, ytickz11_1, xticklabelz11_1, yticklabelz11_1, 10.5);
    MagInset(gcf, -1, ...
        [xz11_2, yz11_2], ...
        [xr11_2, yr11_2], {'SE', 'SE'; 'SW', 'SW'}, cp, xlimz11_2, ylimz11_2, xtickz11_2, ytickz11_2, xticklabelz11_2, yticklabelz11_2, 10.5);
        ax = gca;
    xp = ax.XLim(1);
    yp = ax.YLim(2);
    text(xp, yp, "$\times 10^{-3}$", ...
        'FontSize', 10.5, 'Interpreter', 'latex', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    %----------------------------------------------------------------------
    if export(1, 1)
        exportgraphics(fig11, fid11, 'ContentType', 'Vector', 'Resolution', 600);
        movefile      (fid11, path3);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~export(1, 1)
        fig12 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
    else
        fig12 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
    end
    hold on;
    %----------------------------------------------------------------------
    xline(44.52, '--', "$t_{b}^{\prime}$", ...
        'Color', 'k', ...
        'LineWidth', 0.5, ...
        'FontSize', 22.5, ...
        'LabelVerticalAlignment', 'bottom', ...
        'LabelOrientation', 'horizontal');
    %----------------------------------------------------------------------
    for i = 2:ng
        P12{1, i-1} = plot(Fe1{i, 1}(:, 1), Fe1{i, 1}(:, 2), M(1, 1), 'Color', Color(i, :), 'LineWidth', lw, 'MarkerFaceColor', 'w', 'MarkerSize', 10);
    end
    for i = ng:ng+2
        P12{1, i} = plot(nan, nan, M(1, i-4), ...
            'Color', 'k', ...
            'LineWidth', LW);
    end
    for i = 3:1+ng
        for j = 1:n
            plot(Fh3(:, 1, j, 3), Fh3(:, i, j, 3), M(1, 1+j), 'Color', Color(i-1, :), 'LineWidth', LW);
        end
    end
    %----------------------------------------------------------------------
    PlotX1(...
        xlim12, ylim11, ...
        xtick12, ytick11, ...
        xticklabel12, yticklabel11, ...
        xtitle1, ytitle11);
    PlotX2(...
        P12, legendlabel12, 'northeast', numcols1);
    %----------------------------------------------------------------------
    if export(1, 1)
        exportgraphics(fig12, fid12, 'ContentType', 'Vector', 'Resolution', 600);
        movefile      (fid12, path3);
    end
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 2)
    if ~export(1, 1)
        fig2 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
    else
        fig2 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
    end
    hold on;
    for i = 1:ne
        P2{1, i} = plot(Fe2{i, 1}(:, 1), Fe2{i, 1}(:, 2), M(1, 1), 'Color', Color(i, :), 'LineWidth', lw, 'MarkerFaceColor', 'w', 'MarkerSize', 10);
    end
    for i = 1:n
        P2{1, ne+i} = plot(nan, nan, M(1, 1+i), 'Color', 'k', 'LineWidth', LW);
    end
    for i = 1:ne
        for j = 1:n
            plot(Xc1{i, 1}(:, 1), Fh2{i, j}(:, 1), M(1, 1+j), 'Color', Color(i, :), 'LineWidth', LW);
        end
    end
    %----------------------------------------------------------------------
    PlotX1(...
        xlim2, ylim2, ...
        xtick2, ytick2, ...
        xticklabel2, yticklabel2, ...
        xtitle1, ytitle2);
    PlotX2(...
        P2(1, :), legendlabel2, 'northeast', numcols2);
    %----------------------------------------------------------------------
    if export(1, 1)
        exportgraphics(fig2, fid13, 'ContentType', 'Vector', 'Resolution', 600);
        movefile      (fid13, path3);
    end
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 3)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fig3 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
%     cp = get(gca, 'Position');
%     hold on;
%     %----------------------------------------------------------------------
%     xbathy = linspace(-0.5./h0, 0.5./h0, 1000);
%     switch bathy
%         case 0
%             ybathy = ZbC0(xbathy);
%         case 3
%             ybathy = ZbC3(xbathy);
%         otherwise
%             return
%     end
%     patch(...
%         [xm./h0, -0.5./h0, xbathy, 0.5./h0, xM./h0, xM./h0], [-1.0, -1.0, ybathy, -1.0+0.5./h0.*slope, 0.0, -1.0], grey, 'EdgeColor', grey);
%     %----------------------------------------------------------------------
%     NK       = linspace(0, 150, ns)';
%     nk       = round(NK(1:ns-1, 1), 0);
%     nk(1, 1) = 1;
%     for i = 1:ns-1
%         %P33{1, i} = plot(xc./h0, Fh1(nk(i, 1), :, 1, 3)', '-', 'Color', lred  (11, :), 'LineWidth', LW, 'Visible', 'off');
%         P33{1, i} = plot(xc./h0, Fh1(nk(i, 1), :, 3, 3)', '-', 'Color', lblue (11, :), 'LineWidth', LW, 'Visible', 'off');
%         %P33{3, i} = plot(xc./h0, Fh1(nk(i, 1), :, 3, 3)', '-', 'Color', lgreen(11, :), 'LineWidth', LW, 'Visible', 'off');
%     end
%     %for i = 1:n
%         P31{1, 1} = plot(xc./h0, Fh1(1, :, 3, 3)', '-', 'Color', Color(2, :), 'LineWidth', LW);
%         P32{1, 1} = plot(nan   , nan          , '-', 'Color', Color(2, :), 'LineWidth', LW);
%     %end
%     %----------------------------------------------------------------------
%     for i = 1:n
%         set(P33{1, 1}, 'Visible', 'on');
%     end
%     %----------------------------------------------------------------------
%     PlotX1(...
%         xlim10, ylim3, ...
%         xtick10, ytick3, ...
%         xticklabel10, yticklabel3, ...
%         xtitle1, ytitle10);
%     %PlotX2(...
%     %    P32(1, :), legendlabel3, 'northwest', numcols1);
%     %----------------------------------------------------------------------
%     MagInset(gcf, -1, ...
%         [xz10, yz3], ...
%         [xr10, yr3], {'SE', 'SE'; 'SW', 'SW'}, cp, xlimz10, ylimz3, xtickz10, ytickz3, xticklabelz10, yticklabelz3, 10.5);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     video.path_mp4 = sprintf("%s/Grilli1 (C%d)", path2, bathy);
%     video.path_png = sprintf("%s/FRAMES/", path2);
%     video.fig      = fig3;
%     video.P1       = P31;
%     video.P2       = P33;
%     video.Fh1      = Fh1(:, :, :, 3);
%     video.nf       = nf;
%     video.nk       = repmat(nk', n, 1);
%     %----------------------------------------------------------------------
%     WriteVideo2(video, export(1, 2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~export(1, 2)
        fig3 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
    else
        fig3 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
    end
    hold on;
    %----------------------------------------------------------------------
    NK       = linspace(0, 150, ns)';
    nk       = round(NK(1:ns-1, 1), 0);
    nk(1, 1) = 1;
    for i = 1:ns-1
    %    P33{1, i} = plot(xc./h0, Fh1(nk(i, 1), :, 1, 3)', '-', 'Color', lred  (11, :), 'LineWidth', LW, 'Visible', 'off');
        P33{1, i} = plot(xc./h0, Fh1(nk(i, 1), :, 3, 3)', '-', 'Color', lblue (11, :), 'LineWidth', LW, 'Visible', 'off');
    %    P33{3, i} = plot(xc./h0, Fh1(nk(i, 1), :, 3, 3)', '-', 'Color', lgreen(11, :), 'LineWidth', LW, 'Visible', 'off');
    end
    %for i = 1:n
        P31{1, 1} = plot(xc./h0, Fh1(1, :, 3, 3)', '-', 'Color', Color(2, :), 'LineWidth', LW);
        %P32{1, 1} = plot(nan   , nan          , '-', 'Color', Color(1, :), 'LineWidth', LW);
    %end
    %----------------------------------------------------------------------
    %for i = 1:n
        set(P33{1, 1}, 'Visible', 'on');
    %end
    %----------------------------------------------------------------------
    PlotX1(...
        xlim10, ylim11, ...
        xtick10, ytick11, ...
        xticklabel10, yticklabel11, ...
        xtitle1, ytitle11);
    %PlotX2(...
    %    P32(1, :), legendlabel3, 'northwest', numcols1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    video.path_mp4 = sprintf("%s/Grilli3 (C%d)", path2, bathy);
    video.path_png = sprintf("%s/FRAMES/", path2);
    video.fig      = fig3;
    video.P1       = P31;
    video.P2       = P33;
    video.Fh1      = Fh1(:, :, :, 3);
    video.nf       = nf;
    video.nk       = repmat(nk', n, 1);
    %----------------------------------------------------------------------
    WriteVideo2(video, export(1, 2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 4)
    %----------------------------------------------------------------------
    G  = 9.81;
    A  = 0.20;
    c  = sqrt(G.*(h0+A));
    xs =-5.2.*h0-(c.*13.56.*(sqrt(h0./G))); 
    ts =-xs./c;
    is = sum(th1(:, 3, 3) <= ts, 1); % note: Fh4 starts in t = t1.
    %----------------------------------------------------------------------
    fig41 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
    cp = get(gca, 'Position');
    hold on;
    %----------------------------------------------------------------------
    for i = 1:4
        P4{1, i} = plot(xc./h0, Fh4(1250, :, i)', '-', 'Color', Color(i, :), 'LineWidth', LW);
    end
    %----------------------------------------------------------------------
    PlotX1(...
        xlim10, ylim4, ...
        xtick10, ytick4, ...
        xticklabel10, yticklabel4, ...
        xtitle1, ytitle4);
    PlotX2(...
        P4(1, :), legendlabel10, 'northwest', numcols1);
    %----------------------------------------------------------------------
    MagInset(gcf, -1, ...
        [xz4, yz4], ...
        [xr4, yr4], {'SE', 'SE'; 'SW', 'SW'}, cp, xlimz4, ylimz4, xtickz4, ytickz4, xticklabelz4, yticklabelz4, 10.5);
    %----------------------------------------------------------------------
    if export(1, 1)
        exportgraphics(fig41, fid14, 'ContentType', 'Vector', 'Resolution', 600);
        movefile      (fid14, path3);
    end
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end