function [] = PlotT11C(g, export)
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
M     = ["-", "-.", "--"];
LW    = 2.5;
flag  = [0, 1];
bathy = 3;
path1 = "Plot/T11 - Reversibility check/DATA";
path2 = "Plot/T11 - Reversibility check";
path3 = "Plot/T11 - Reversibility check/";
path4 = sprintf("C%d", bathy);
%--------------------------------------------------------------------------
xm    = 0;
xM    = 240;
K     = 2400;
xv    = linspace(xm, xM, K+1)';
xc    = zeros(K, 1);
for i = 1:K
    xc(i, 1) = 1./2.*(xv(i, 1)+xv(i+1, 1));
end
%--------------------------------------------------------------------------
n     = 3;
p     = 3;
Fh1   = zeros(nf, K, n);
for i = 1:n
    Fh1(:, :, i) = getdata1(sprintf("%s/%s/Shelf_P%d_01_GN.mat", path1, path4, i), i, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim         = [0, 240];
xtick        = 0:30:240;
xticklabel   = ["0", "30", "60", "90", "120", "150", "180", "210", "240"];
%--------------------------------------------------------------------------
ylim1        = [-0.02, 0.40];
ytick1       = [-0.02, 0.00:0.05:0.40];
yticklabel1  = ["-0.02", "0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40"];
ylim2        = [-1.00, 0.40];
ytick2       =-1.00:0.20:0.40;
yticklabel2  = ["-1.00", "-0.80", "0.60", "-0.40", "-0.20", "0.00", "0.20", "0.40"];
%--------------------------------------------------------------------------
legendlabel  = ["$k=1$", "$k=2$", "$k=3$"];
numcols      = 1;
xtitle       = "$x$";
ytitle1      = "$\zeta^{\prime}$";
ytitle2      = "$z^{\prime}$";
%--------------------------------------------------------------------------
fid01        = [sprintf("Shelf2_P%d (t0).pdf", p), sprintf("Shelf2_P%d (tf).pdf", p)];
fid02        = [sprintf("Shelf4_P%d (t0).pdf", p), sprintf("Shelf4_P%d (tf).pdf", p)];
fid1         = "Shelf.pdf";
%--------------------------------------------------------------------------
P1           = cell(1, n);
P2           = cell(2, 1);
P3           = cell(2, ns);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlimz1       = [120, 150];
xtickz1      = 120:10:150;
xticklabelz1 = ["120", "130", "140", "150"];
xlimz2       = [228, 236];
xtickz2      = 228:2:236;
xticklabelz2 = ["228", "230", "232", "234", "235", "236"];
%--------------------------------------------------------------------------
ylimz1       = [-5.0, 5.0].*10.^-3;
ytickz1      = (-5.0:2.5:5.0).*10.^-3;
yticklabelz1 = ["-5.0", "-2.5", "0.0", "2.5", "5.0"];
ylimz2       = [0.00, 0.35];
ytickz2      = 0.00:0.05:0.35;
yticklabelz2 = ["0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30", "0.35"];
%--------------------------------------------------------------------------
zm           = 120;
zf           = 1.025./2;
%--------------------------------------------------------------------------
xz1          = [  0,  30].*zf+zm;
xz2          = [108, 116].*zf+zm;
xr1          = [-30,  60].*zf+zm;
xr2          = [-50,  80].*zf+zm;
%--------------------------------------------------------------------------
yz1          = ylimz1;
yz2          = ylimz2;
yr1          = [0.025, 0.025+( 90./240.*0.42)./2];
yr2          = [0.150, 0.150+(130./240.*0.42)   ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig1 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
    cp = get(gca, 'Position');
    hold on;
    %----------------------------------------------------------------------
    for j = 1:n
        P1{1, j} = plot(xc, Fh1(nf, :, j)', M(1, j), 'Color', Color(1+j, :), 'LineWidth', LW);
    end
    %----------------------------------------------------------------------
    PlotX1(...
        xlim, ylim1, ...
        xtick, ytick1, ...
        xticklabel, yticklabel1, ...
        xtitle, ytitle1);
    PlotX2(...
        P1(1, :), legendlabel, 'northwest', numcols);
    %----------------------------------------------------------------------
    MagInset(gcf, -1, ...
        [xz1, yz1], ...
        [xr1, yr1], {'SE', 'SE'; 'SW', 'SW'}, cp, xlimz1, ylimz1, xtickz1, ytickz1, xticklabelz1, yticklabelz1, 10.5);
    ax = gca;
    xp = ax.XLim(1);
    yp = ax.YLim(2);
    text(xp, yp, "$\times 10^{-3}$", ...
        'FontSize', 10.5, 'Interpreter', 'latex', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    MagInset(gcf, -1, ...
        [xz2, yz2], ...
        [xr2, yr2], {'NE', 'NE'; 'NW', 'NW'}, cp, xlimz2, ylimz2, xtickz2, ytickz2, xticklabelz2, yticklabelz2, 10.5);
    %----------------------------------------------------------------------
    if export(1, 1)
        exportgraphics(fig1, fid1, 'ContentType', 'Vector', 'Resolution', 600);
        movefile      (fid1, path3);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:2
        if ~export(1, 2)
            fig21 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        else
            fig21 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
        end
        hold on;
        %------------------------------------------------------------------
        for j = 1:ns
            P3{i, j} = plot(xc, Fh1(nk(1, j), :, p)', '-', 'Color', lblue(11, :), 'LineWidth', LW, 'Visible', 'off');
        end
        P2{i, 1} = plot(xc, Fh1(1, :, p)', '-', 'Color', Color(2, :), 'LineWidth', LW);
        patch(...
            [ 130,  140,  240,  240], ...
            [-1.0, -0.5, -0.5, -1.0], grey, 'EdgeColor', grey);
        %------------------------------------------------------------------
        if i == 2
            set(P3{i, 1}, 'Visible', 'on');
        end
        %------------------------------------------------------------------
        PlotX1(...
            xlim, ylim2, ...
            xtick, ytick2, ...
            xticklabel, yticklabel2, ...
            xtitle, ytitle2);
        %------------------------------------------------------------------
        if i == 2 && export(1, 1)
            for j = 2:ns
                set(P3{i, j}, 'Visible', 'on');
            end
            exportgraphics(fig21, fid01(1, 1), 'ContentType', 'Vector', 'Resolution', 600);
            movefile      (fid01(1, 1), path3);
            for j = 2:ns
                set(P3{i, j}, 'Visible', 'off');
            end
        end
        %------------------------------------------------------------------
        video.path_mp4 = sprintf("%s/Shelf%d_P%d", path2, i, p);
        video.path_png = sprintf("%s/FRAMES/", path2);
        video.fig      = fig21;
        video.P1       = P2(i, 1);
        video.P2       = P3(i, :);
        video.Fh1      = Fh1(:, :, p);
        video.nf       = nf;
        video.nk       = nk;
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
            exportgraphics(fig21, fid01(1, 2), 'ContentType', 'Vector', 'Resolution', 600);
            movefile      (fid01(1, 2), path3);
        end
        %------------------------------------------------------------------
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:2
        %------------------------------------------------------------------
        if ~export(1, 2)
            fig22 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        else
            fig22 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
        end
        hold on;
        %------------------------------------------------------------------
        for j = 1:ns
            P3{i, j} = plot(xc, Fh1(nk(1, j), :, p)', '-', 'Color', lblue(11, :), 'LineWidth', LW, 'Visible', 'off');
        end
        P2{i, 1} = plot(xc, Fh1(1, :, p)', '-', 'Color', Color(2, :), 'LineWidth', LW);
        %------------------------------------------------------------------
        if i == 2
            set(P3{i, 1}, 'Visible', 'on');
        end
        %------------------------------------------------------------------
        PlotX1(...
            xlim, ylim1, ...
            xtick, ytick1, ...
            xticklabel, yticklabel1, ...
            xtitle, ytitle1);
        %------------------------------------------------------------------
        if i == 2 && export(1, 1)
            for j = 2:ns
                set(P3{i, j}, 'Visible', 'on');
            end
            exportgraphics(fig22, fid02(1, 1), 'ContentType', 'Vector', 'Resolution', 600);
            movefile      (fid02(1, 1), path3);
            for j = 2:ns
                set(P3{i, j}, 'Visible', 'off');
            end
        end
        %------------------------------------------------------------------
        video.path_mp4 = sprintf("%s/Shelf%d_P%d", path2, 2+i, p);
        video.path_png = sprintf("%s/FRAMES/", path2);
        video.fig      = fig22;
        video.P1       = P2(i, 1);
        video.P2       = P3(i, :);
        video.Fh1      = Fh1(:, :, p);
        video.nf       = nf;
        video.nk       = nk;
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