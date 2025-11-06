function [] = PlotT08C(g, export)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color = linspecer(9, 'Qualitative');
nc    = 50;
nf    = g.data.nf+1;
nk    = linspace(1, nf, 10);
R     = linspace(1, Color(2, 1), nc);
G     = linspace(1, Color(2, 2), nc);
B     = linspace(1, Color(2, 3), nc);
lblue = [R', G', B'];
M     = [":", "-"];
LW    = 2.5;
flag  = [1, 1];
path1 = "Plot/T08 - Overtaking colission/DATA";
path2 = "Plot/T08 - Overtaking colission";
path3 = "Plot/T08 - Overtaking colission/";
%--------------------------------------------------------------------------
xm    =-200;
xM    = 200;
K     = 4000;
xv    = linspace(xm, xM, K+1)';
xc    = zeros(K, 1);
for i = 1:K
    xc(i, 1) = 1./2.*(xv(i, 1)+xv(i+1, 1));
end
%--------------------------------------------------------------------------
p     = 3;
Fh1   = getdata1(sprintf("%s/Oc_P%d_01_GN.mat", path1, p), p, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim        = [-200, 200];
xtick       =-200:50:200;
xticklabel  = ["-200", "-150", "-100", "-50", "0", "50", "100", "150", "200"];
%--------------------------------------------------------------------------
ylim        = [-0.10, 1.15];
ytick       = [-0.10, 0.00:0.25:1.00, 1.15];
yticklabel  = ["-0.10", "0.00", "0.25", "0.50", "0.75", "1.00", "1.15"];
%--------------------------------------------------------------------------
legendlabel = ["$t^{\prime} = t_{0}^{\prime}$", "$t^{\prime} = t_{f}^{\prime}$"];
numcols     = 1;
xtitle1     = "$x-ct^{\prime}$";
xtitle2     = "$x$";
ytitle      = "$\zeta^{\prime}$";
%--------------------------------------------------------------------------
fid1        = sprintf("Oc_P%d.pdf", p);
%--------------------------------------------------------------------------
P1          = cell(1, 2);
P2          = cell(2, 1);
P3          = cell(2, numel(nk));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~export(1, 1)
        fig1 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
    else
        fig1 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
    end
    hold on;
    %----------------------------------------------------------------------
    xline(0, '--', "$t_{o}^{\prime}$", ...
        'Color', 'k', ...
        'LineWidth', 0.5, ...
        'FontSize', 22.5, ...
        'LabelVerticalAlignment', 'bottom', ...
        'LabelOrientation', 'horizontal');
    %----------------------------------------------------------------------
    for i = 2:numel(nk)-1
        plot(xc, Fh1(nk(1, i), :)', '-', 'Color', lblue(11, :), 'LineWidth', LW);
    end
    P1{1, 1} = plot(xc, Fh1(1  , :)', M(1, 1), 'Color', Color(2, :), 'LineWidth', LW);
    P1{1, 2} = plot(xc, Fh1(end, :)', M(1, 2), 'Color', Color(2, :), 'LineWidth', LW);
    %----------------------------------------------------------------------
    PlotX1(...
        xlim, ylim, ...
        xtick, ytick, ...
        xticklabel, yticklabel, ...
        xtitle1, ytitle);
    PlotX2(...
        P1(1, :), legendlabel, 'northeast', numcols);
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
            fig2 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        else
            fig2 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
        end
        hold on;
        %------------------------------------------------------------------
        for j = 1:numel(nk)
            P3{i, j} = plot(xc, Fh1(nk(1, j), :)', '-', 'Color', lblue(11, :), 'LineWidth', LW, 'Visible', 'off');
        end
        P2{i, 1} = plot(xc, Fh1(1, :)', '-', 'Color', Color(2, :), 'LineWidth', LW);
        %------------------------------------------------------------------
        if i == 2
            set(P3{i, 1}, 'Visible', 'on');
        end
        %------------------------------------------------------------------
        PlotX1(...
            xlim, ylim, ...
            xtick, ytick, ...
            xticklabel, yticklabel, ...
            xtitle2, ytitle);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        video.path_mp4 = sprintf("%s/Oc%d_P%d", path2, i, p);
        video.path_png = sprintf("%s/FRAMES/", path2);
        video.fig      = fig2;
        video.P1       = P2(i, 1);
        video.P2       = P3(i, :);
        video.Fh1      = Fh1;
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
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end