function [] = PlotT09C(g, export)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color = linspecer(9, 'Qualitative');
nc    = 50;
ns    = 3;
nf    = g.data.nf+1;
R1    = linspace(1, Color(1, 1), nc);
G1    = linspace(1, Color(1, 2), nc);
B1    = linspace(1, Color(1, 3), nc);
R2    = linspace(1, Color(2, 1), nc);
G2    = linspace(1, Color(2, 2), nc);
B2    = linspace(1, Color(2, 3), nc);
lred  = [R1', G1', B1'];
lblue = [R2', G2', B2'];
M     = [":", "-", "-."];
LW    = 2.5;
flag  = [1, 1];
path1 = "Plot/T09 - Breakup of a Gaussian hump/DATA";
path2 = "Plot/T09 - Breakup of a Gaussian hump";
path3 = "Plot/T09 - Breakup of a Gaussian hump/";
%--------------------------------------------------------------------------
xm    =-100;
xM    = 100;
K     = 2000;
xv    = linspace(xm, xM, K+1)';
xc    = zeros(K, 1);
for i = 1:K
    xc(i, 1) = 1./2.*(xv(i, 1)+xv(i+1, 1));
end
%--------------------------------------------------------------------------
m     = 2;
p     = 3;
Fh1   = zeros(nf, K  , m);
Fh2   = zeros(ns, 1+K, m);
th1   = zeros(nf, m);
str1  = ["01", "20"];
for i = 1:m
    Fh1(:, :, i) = getdata1(sprintf("%s/%s/Bgh_P%d_%s_05_GN.mat", path1, str1(1, i), p, str1(1, i)), p, 1);
    Fh2(:, :, i) = getdata2(sprintf("%s/%s/Bgh_P%d_%s_05_GN.mat", path1, str1(1, i), p, str1(1, i)), p, 1, ns);
    th1(:, i)    = getdata3(sprintf("%s/%s/Bgh_P%d_%s_05_GN.mat", path1, str1(1, i), p, str1(1, i)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim         = [-100, 100];
xtick        =-100:25:100;
xticklabel   = ["-100", "-75", "-50", "-25", "0", "25", "50", "75", "100"];
%--------------------------------------------------------------------------
ylim         = [-0.25, 1.65];
ytick        = [-0.25, 0.00:0.25:1.50, 1.65];
yticklabel   = ["-0.25", "0.00", "0.25", "0.50", "0.75", "1.00", "1.25", "1.50", "1.65"];
%--------------------------------------------------------------------------
legendlabel1 = ["$\lambda = 1$", "$\lambda = 20$"];
legendlabel2 = ["$t^{\prime} = t_{0}^{\prime}$", "$t^{\prime} = t_{f}^{\prime}$", "$t^{\prime} = t_{r}^{\prime}$"];
%--------------------------------------------------------------------------
numcols      = 1;
xtitle1      = "$x$";
xtitle2      = "$x-ct^{\prime}$";
ytitle       = "$\zeta^{\prime}$";
%--------------------------------------------------------------------------
fid1         = [...
    "Bgh.pdf", ...
    sprintf("Bgh_P%d (01).pdf", p), ...
    sprintf("Bgh_P%d (20).pdf", p)];
%--------------------------------------------------------------------------
P1           = cell(ns, m);
P2_1         = cell(1 , m);
P2_2         = cell(1 , m);
P3           = cell(m , ns);
Z1           = @(x) exp(-x.^2);
Z2           = @(x) exp(-x.^2./20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~export(1, 1)
        fig11 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
    else
        fig11 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
    end
    hold on;
    %----------------------------------------------------------------------
    P1{1, 1} = fplot(@(x) Z1(x), '-', 'Color', Color(1, :), 'LineWidth', LW);
    P1{1, 2} = fplot(@(x) Z2(x), '-', 'Color', Color(2, :), 'LineWidth', LW);
    %----------------------------------------------------------------------
    PlotX1(...
        xlim, ylim, ...
        xtick, ytick, ...
        xticklabel, yticklabel, ...
        xtitle1, ytitle);
    PlotX2(...
        P1(1, :), legendlabel1, 'northeast', numcols);
    %----------------------------------------------------------------------
    if export(1, 1)
        exportgraphics(fig11, fid1(1, 1), 'ContentType', 'Vector', 'Resolution', 600);
        movefile      (fid1(1, 1), path3);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:m
        %------------------------------------------------------------------
        if ~export(1, 1)
            fig12 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        else
            fig12 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
        end
        hold on;
        %------------------------------------------------------------------
        P1{i, 1} = plot(xc, Fh2(1  , 2:K+1, i)', M(1, 1), 'Color', Color(i, :), 'LineWidth', LW);
        P1{i, 2} = plot(xc, Fh1(end,     :, i)', M(1, 2), 'Color', Color(i, :), 'LineWidth', LW);
        switch i
            case 1
                [~, j] = ...
                    max(Fh1(:, 1, 1)); % FIND WALL TIME
                P1{i, 3} = plot(xc, Fh1(j,     :, 1)', M(1, 3), 'Color', Color(i, :), 'LineWidth', LW);
            case 2
                P1{i, 3} = plot(xc, Fh2(3, 2:K+1, 2)', M(1, 3), 'Color', Color(i, :), 'LineWidth', LW);
            otherwise
                return
        end
        %------------------------------------------------------------------
        PlotX1(...
            xlim, ylim, ...
            xtick, ytick, ...
            xticklabel, yticklabel, ...
            xtitle2, ytitle);
        PlotX2(...
            P1(i, :), legendlabel2, 'northeast', numcols);
        %------------------------------------------------------------------
        if export(1, 1)
            exportgraphics(fig12, fid1(1, 1+i), 'ContentType', 'Vector', 'Resolution', 600);
            movefile      (fid1(1, 1+i), path3);
        end
        %------------------------------------------------------------------
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 2:2
        %------------------------------------------------------------------
        if ~export(1, 2)
            fig2 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        else
            fig2 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
        end
        hold on;
        %------------------------------------------------------------------
        for j = 1:ns
            if j ~= 3
                P3{1, j} = plot(xc, Fh2(j, 2:K+1, 1)', '-', 'Color', lred (11, :), 'LineWidth', LW, 'Visible', 'off');
            else
                [~, k] = ...
                    max(Fh1(:, 1, 1)); % FIND WALL TIME
                P3{1, j} = plot(xc, Fh1(k, :, 1)', '-', 'Color', lred (11, :), 'LineWidth', LW, 'Visible', 'off');
            end
            P3{2, j} = plot(xc, Fh2(j, 2:K+1, 2)', '-', 'Color', lblue(11, :), 'LineWidth', LW, 'Visible', 'off');
        end
        for j = 1:m
            P2_1{1, j} = plot(xc, Fh1(1, :, j)', '-', 'Color', Color(j, :), 'LineWidth', LW);
            P2_2{1, j} = plot(nan, nan, '-', 'Color', Color(j, :), 'LineWidth', LW);
        end
        %------------------------------------------------------------------
        if i == 2
            for j = 1:m
                set(P3{j, 1}, 'Visible', 'on');
            end
        end
        %------------------------------------------------------------------
        PlotX1(...
            xlim, ylim, ...
            xtick, ytick, ...
            xticklabel, yticklabel, ...
            xtitle1, ytitle);
        PlotX2(...
            P2_2, legendlabel1, 'northeast', numcols);
        %------------------------------------------------------------------
        video.path_mp4 = sprintf("%s/Bgh%d_P%d", path2, i, p);
        video.path_png = sprintf("%s/FRAMES/", path2);
        video.fig      = fig2;
        video.P1       = P2_1;
        video.P2       = P3;
        video.Fh1      = Fh1;
        video.nf       = nf;
        video.nk       = [sum(Fh2(:, 1, 1)' > th1(:, 1), 1); sum(Fh2(:, 1, 2)' > th1(:, 2), 1)];
        video.nk(1, 3) = k;
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
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end