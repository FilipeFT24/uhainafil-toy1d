function [] = PlotS2T1(g, export)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xx = 1;

aa = load("Plot/S2/T1 - Head-on-colission/DATA/0.21/Hoc_P2.mat");


bb = reshape(aa.U(200, :, :, 1), [500, 3]);

Color  = linspecer(9, 'Qualitative');
grey   = repmat(0.80, 1, 3);
nc     = 50;
ns     = 3;
nf     = g.data.nf+numel(g.data.tk);
R      = linspace(1, Color(2, 1), nc);
G      = linspace(1, Color(2, 2), nc);
B      = linspace(1, Color(2, 3), nc);
lblue  = [R', G', B'];
M      = [":", "-", "-."];
LW     = 2.5;
flag   = [0, 1, 0];
path11 = "Plot/T06 - Head-on-colission/DATA";
path12 = "Plot/T06 - Head-on-colission";
path13 = "Plot/T06 - Head-on-colission/";
path21 = "Plot/T07 - Wall/DATA";
path22 = "Plot/T07 - Wall";
path_1 = [path11; path21];
path_2 = [path12; path22];
%--------------------------------------------------------------------------
xm     =-100;
xM1    = 100;
xM2    = 0;
K1     = 2000;
K2     = 1000;
xv1    = linspace(xm, xM1, K1+1)';
xv2    = linspace(xm, xM2, K2+1)';
xc1    = zeros(K1, 1);
xc2    = zeros(K2, 1);
for i = 1:K1
    xc1(i, 1) = 1./2.*(xv1(i, 1)+xv1(i+1, 1));
end
for i = 1:K2
    xc2(i, 1) = 1./2.*(xv2(i, 1)+xv2(i+1, 1));
end
%--------------------------------------------------------------------------
m      = 2;
p      = 3;
Fh2C   = zeros(ns, 1+K1, m);
Fh1C   = zeros(nf, K1  , m);
th1C   = zeros(nf, m);
Fh2W   = zeros(ns, 1+K2, m);
Fh1W   = zeros(nf, K2  , m);
th1W   = zeros(nf, m);
str1   = ["96", "21"];
str2   = ["Hoc", "Wall"];
for i = 1:m
    Fh1C(:, :, i) = getdata1(sprintf("%s/0.%s/%s_P%d_%s_05_GN.mat", path11, str1(1, i), str2(1, 1), p, str1(1, i)), p, 1);
    Fh2C(:, :, i) = getdata2(sprintf("%s/0.%s/%s_P%d_%s_05_GN.mat", path11, str1(1, i), str2(1, 1), p, str1(1, i)), p, 1, ns);
    th1C(:, i)    = getdata3(sprintf("%s/0.%s/%s_P%d_%s_05_GN.mat", path11, str1(1, i), str2(1, 1), p, str1(1, i)));
    Fh1W(:, :, i) = getdata1(sprintf("%s/0.%s/%s_P%d_%s_05_GN.mat", path21, str1(1, i), str2(1, 2), p, str1(1, i)), p, 1);
    Fh2W(:, :, i) = getdata2(sprintf("%s/0.%s/%s_P%d_%s_05_GN.mat", path21, str1(1, i), str2(1, 2), p, str1(1, i)), p, 1, ns);
    th1W(:, i)    = getdata3(sprintf("%s/0.%s/%s_P%d_%s_05_GN.mat", path21, str1(1, i), str2(1, 2), p, str1(1, i)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1) || flag(1, 2)
    %----------------------------------------------------------------------
    xlim11       = [-100, 100];
    xtick11      =-100:25:100;
    xticklabel11 = ["-100", "-75", "-50", "-25", "0", "25", "50", "75", "100"];
    xlim12       = [-100, 2.50];
    xtick12      =-100:25:0;
    xticklabel12 = ["-100", "-75", "-50", "-25", "0"];
    %----------------------------------------------------------------------
    xlim1        = {xlim11; xlim12};
    xtick1       = {xtick11; xtick12};
    xticklabel1  = {xticklabel11; xticklabel12};
    %----------------------------------------------------------------------
    xtitle1      = "$x$";
    xtitle2      = "$x-ct^{\prime}$";
    ytitle1      = "$\zeta^{\prime}$";
    numcols1     = 1;
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P11           = cell(2, m);
    P12           = cell(m, 3);
    fid11         = [sprintf("Hoc_P%d.pdf", p); sprintf("Wall_P%d.pdf", p)];
    fid12         = [sprintf("HocWall1_P%d (0.%s).pdf", p, str1(1, 1)); sprintf("HocWall1_P%d (0.%s).pdf", p, str1(1, 2))];
    %----------------------------------------------------------------------
    ylim10        = [-0.05, 1.15];
    ytick10       = [-0.05, 0.00:0.25:1.00, 1.15];
    yticklabel10  = ["-0.05", "0.00", "0.25", "0.50", "0.75", "1.00", "1.15"];
    ylim11        = [-0.15, 3.00];
    ytick11       = [-0.15, 0.00:0.50:3.00];
    yticklabel11  = ["-0.15", "0.00", "0.50", "1.00", "1.50", "2.00", "2.50", "3.00"];
    ylim12        = [-0.02, 0.50];
    ytick12       = [-0.02, 0.00:0.10:0.50];
    yticklabel12  = ["-0.02", "0.00", "0.10", "0.20", "0.30", "0.40", "0.50"];
    %----------------------------------------------------------------------
    ylim1         = {ylim11; ylim12};
    ytick1        = {ytick11; ytick12};
    yticklabel1   = {yticklabel11; yticklabel12};
    %----------------------------------------------------------------------
    legendlabel11 = ["$\epsilon = 0.96$", "$\epsilon = 0.21$"];
    legendlabel12 = ["$t^{\prime} = t_{0}^{\prime}$", "$t^{\prime} = t_{f}^{\prime}$", "$t^{\prime} = t_{p}^{\prime}$"];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:2
        %------------------------------------------------------------------
        if ~export(1, 1)
            fig11 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        else
            fig11 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
        end
        hold on;
        %------------------------------------------------------------------
        switch i
            case 1
                P11{i, 1} = plot(xc1, Fh1C(1, :, 1)', '-', 'Color', Color(1, :), 'LineWidth', LW);
                P11{i, 2} = plot(xc1, Fh1C(1, :, 2)', '-', 'Color', Color(2, :), 'LineWidth', LW);
            case 2
                patch(...
                    [ 0.00, 0.00, 2.50,  2.50], ...
                    [-0.05, 1.15, 1.15, -0.05], grey, 'EdgeColor', grey);
                P11{i, 1} = plot(xc2, Fh1W(1, :, 1)', '-', 'Color', Color(1, :), 'LineWidth', LW);
                P11{i, 2} = plot(xc2, Fh1W(1, :, 2)', '-', 'Color', Color(2, :), 'LineWidth', LW);
            otherwise
                return
        end
        %------------------------------------------------------------------
        PlotX1(...
            xlim1{i, 1}, ylim10, ...
            xtick1{i, 1}, ytick10, ...
            xticklabel1{i, 1}, yticklabel10, ...
            xtitle1, ytitle1);
        PlotX2(...
            P11(i, :), legendlabel11, 'northeast', numcols1);
        %------------------------------------------------------------------
        if export(1, 1)
            exportgraphics(fig11, fid11(i, 1), 'ContentType', 'Vector', 'Resolution', 600);
            movefile      (fid11(i, 1), path13);
        end
        %------------------------------------------------------------------
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
        patch(...
            [ 0.00, 0.00, 5.00,  5.00], ...
            [-0.15, 3.00, 3.00, -0.15], grey, 'EdgeColor', grey);
        xline(0, '--', ...
            'Color', 'k');
        %------------------------------------------------------------------
        P12{i, 1} = plot(xc1, [Fh1W(  1,      :, i), Fh1C(  1, K2+1:end, i)]', M(1, 1), 'Color', Color(i, :), 'LineWidth', LW);
        P12{i, 2} = plot(xc1, [Fh1W(end,      :, i), Fh1C(end, K2+1:end, i)]', M(1, 2), 'Color', Color(i, :), 'LineWidth', LW);
        P12{i, 3} = plot(xc1, [Fh2W(  3, 2:K2+1, i), Fh2C(  3, K2+2:end, i)]', M(1, 3), 'Color', Color(i, :), 'LineWidth', LW);
        %------------------------------------------------------------------
        PlotX1(...
            xlim11, ylim1{i, 1}, ...
            xtick11, ytick1{i, 1}, ...
            xticklabel11, yticklabel1{i, 1}, ...
            xtitle2, ytitle1);
        PlotX2(...
            P12(i, :), legendlabel12, 'northeast', numcols1);
        %------------------------------------------------------------------
        if export(1, 1)
            exportgraphics(fig12, fid12(i, 1), 'ContentType', 'Vector', 'Resolution', 600);
            movefile      (fid12(i, 1), path13);
        end
        %------------------------------------------------------------------
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P2           = cell(2 , ns);
    P3           = cell(2 , ns);
    %----------------------------------------------------------------------
    ylim21       = [-0.50, 2.65];
    ytick21      = [-0.50:0.50:2.50, 2.65];
    yticklabel21 = ["-0.50", "0.00", "0.50", "1.00", "1.50", "2.00", "2.50", "2.65"];
    ylim22       = [-0.02, 0.50];
    ytick22      = [-0.02, 0.00:0.10:0.50];
    yticklabel22 = ["-0.02", "0.00", "0.10", "0.20", "0.30", "0.40", "0.50"];
    %----------------------------------------------------------------------
    ylim2        = {ylim21; ylim22};
    ytick2       = {ytick21; ytick22};
    yticklabel2  = {yticklabel21; yticklabel22};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 2:2
        for j = 1:m
            %--------------------------------------------------------------
            if ~export(1, 2)
                fig2 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
            else
                fig2 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
            end
            hold on;
            %--------------------------------------------------------------
            switch i
                case 1
                    P3{i, 1} = plot(xc1, Fh2C(1, 2:K1+1, j), '-', 'Color', lblue(11, :), 'LineWidth', LW, 'Visible', 'off');
                    P3{i, 2} = plot(xc1, Fh2C(2, 2:K1+1, j), '-', 'Color', lblue(11, :), 'LineWidth', LW, 'Visible', 'off');
                    P3{i, 3} = plot(xc1, Fh2C(3, 2:K1+1, j), '-', 'Color', lblue(11, :), 'LineWidth', LW, 'Visible', 'off');
                    P2{i, 1} = plot(xc1, Fh1C(1,      :, j), '-', 'Color', Color( 2, :), 'LineWidth', LW);
                case 2
                    P3{i, 1} = plot(xc2, Fh2W(1, 2:K2+1, j), '-', 'Color', lblue(11, :), 'LineWidth', LW, 'Visible', 'off');
                    P3{i, 2} = plot(xc2, Fh2W(2, 2:K2+1, j), '-', 'Color', lblue(11, :), 'LineWidth', LW, 'Visible', 'off');
                    P3{i, 3} = plot(xc2, Fh2W(3, 2:K2+1, j), '-', 'Color', lblue(11, :), 'LineWidth', LW, 'Visible', 'off');
                    P2{i, 1} = plot(xc2, Fh1W(1,      :, j), '-', 'Color', Color( 2, :), 'LineWidth', LW);
                    patch(...
                        [ 0.00, 0.00, 2.50,  2.50], ...
                        [-0.50, 2.65, 2.65, -0.50], grey, 'EdgeColor', grey);
                otherwise
                    return
            end
            %--------------------------------------------------------------
            PlotX1(...
                xlim1{i, 1}, ylim2{j, 1}, ...
                xtick1{i, 1}, ytick2{j, 1}, ...
                xticklabel1{i, 1}, yticklabel2{j, 1}, ...
                xtitle1, ytitle1);
            %--------------------------------------------------------------
            video.path_mp4 = sprintf("%s/%s1_P%d (0.%s)", path_2(i, 1), str2(1, i), p, str1(1, j));
            video.path_png = sprintf("%s/FRAMES/", path_2(i, 1));
            video.fig      = fig2;
            video.P1       = P2(i, 1);
            video.P2       = P3(i, :);
            video.nf       = nf;
            switch i
                case 1
                    video.Fh1 = Fh1C(:, :, j);
                    video.nk  = sum(Fh2C(:, 1, j)' > th1C(:, j), 1);
                case 2
                    video.Fh1 = Fh1W(:, :, j);
                    video.nk  = sum(Fh2W(:, 1, j)' > th1W(:, j), 1);
                otherwise
                    return
            end
            video.nk(1, 2) = nf+1;
            %--------------------------------------------------------------
            WriteVideo1(video, export(1, 2));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 3)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P4           = cell(2, 4);
    fid3         = [sprintf("Hoc_P%d (xp).pdf", p); sprintf("Wall_P%d (tp).pdf", p)];
    %----------------------------------------------------------------------
    Za           = {...
        za(1, 1, 0.96);
        za(1, 1, 0.21)};
    na           = 1000;
    Fa1          = zeros(na, 2, 2);
    Fa2          = zeros(na, 2, 2);
    h0           = 1;
    A            = [0.96; 0.21].*h0;
    c            = sqrt(G.*(h0+A));
    ta           = [...
        linspace(0, 100./c(1, 1), na)', ...
        linspace(0, 100./c(2, 1), na)'];
    for i = 1:2
        [ZavM, ZaiM] = ...
            max(Za{i, 1}(xc1', ta(:, i)), [], 2);
        Fa1(:, 1, i) = ta(:, i);
        Fa2(:, 1, i) = ta(:, i);
        Fa1(:, 2, i) = abs(xc1(ZaiM, 1));
        Fa2(:, 2, i) = ZavM;
    end
    %----------------------------------------------------------------------
    xlim3        = [0, 100];
    xtick3       = 0:10:100;
    xticklabel3  = ["0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"];
    %----------------------------------------------------------------------
    ylim31       = [-50, 50];
    ytick31      =-50:25:50;
    yticklabel31 = ["-50", "-25", "0", "25", "50"];
    ylim32       = [0.50, 2.65];
    ytick32      = [0.50, 0.65, 1.00:0.50:2.50, 2.65];
    yticklabel32 = ["0.50", "0.65", "1.00", "1.50", "2.00", "2.50", "2.65"];
    ylim33       = [0.15, 0.50];
    ytick33      = [0.15, 0.20:0.10:0.50];
    yticklabel33 = ["0.15", "0.20", "0.30", "0.40", "0.50"];
    %----------------------------------------------------------------------
    xtitle3      = "$t^{\prime}$";
    ytitle31     = "$x_{p}$";
    ytitle32     = "$\zeta_{p}^{\prime}$";
    legendlabel3 = ["$\epsilon = 0.96$", "$\epsilon = 0.21$", "$\textrm{Analytical}$", "$k = 3$"];
    numcols3     = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:2
        %------------------------------------------------------------------
        if ~export(1, 1)
            fig3 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        else
            fig3 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
        end
        hold on;
        %------------------------------------------------------------------
        for j = 1:m
            xline(50./c(j, 1), '--', ...
                'Color', Color(j, :));
            xline(Fh2C(3, 1, j), '-', "$t_{p}^{\prime}$", ...
                'Color', Color(j, :), ...
                'LineWidth', 0.5, ...
                'FontSize', 22.5, ...
                'LabelVerticalAlignment', 'bottom', ...
                'LabelOrientation', 'horizontal');
            xline(100./c(j, 1), '--', "$t_{f}^{\prime}$", ...
                'Color', Color(j, :), ...
                'LineWidth', 0.5, ...
                'FontSize', 22.5, ...
                'LabelVerticalAlignment', 'bottom', ...
                'LabelOrientation', 'horizontal');
        end
        %------------------------------------------------------------------
        [vmax1, imax1] = ...
            max(Fh1W(:, :, 1), [], 2);
        [vmax2, imax2] = ...
            max(Fh1W(:, :, 2), [], 2);
        X1 = xc2(imax1);
        X2 = xc2(imax2);
        V1 = vmax1;
        V2 = vmax2;
        switch i
            case 1
                %----------------------------------------------------------
                PlotX1(...
                    xlim3, ylim31, ...
                    xtick3, ytick31, ...
                    xticklabel3, yticklabel31, ...
                    xtitle3, ytitle31);
                %----------------------------------------------------------
                plot(Fa1 (:, 1, 1),  Fa1(:, 2, 1), M(1, 1), 'Color', Color(1, :), 'LineWidth', LW(1, 1));
                plot(Fa1 (:, 1, 1), -Fa1(:, 2, 1), M(1, 1), 'Color', Color(1, :), 'LineWidth', LW(1, 1));
                plot(Fa1 (:, 1, 2),  Fa1(:, 2, 2), M(1, 1), 'Color', Color(2, :), 'LineWidth', LW(1, 1));
                plot(Fa1 (:, 1, 2), -Fa1(:, 2, 2), M(1, 1), 'Color', Color(2, :), 'LineWidth', LW(1, 1));
                plot(th1W(:, 1)   ,  X1          , M(1, 2), 'Color', Color(1, :), 'LineWidth', LW(1, 1));
                plot(th1W(:, 1)   , -X1          , M(1, 2), 'Color', Color(1, :), 'LineWidth', LW(1, 1));
                plot(th1W(:, 2)   ,  X2          , M(1, 2), 'Color', Color(2, :), 'LineWidth', LW(1, 1));
                plot(th1W(:, 2)   , -X2          , M(1, 2), 'Color', Color(2, :), 'LineWidth', LW(1, 1));
                %----------------------------------------------------------
            case 2
                %----------------------------------------------------------
                PlotX1(...
                    xlim3, ylim32, ...
                    xtick3, ytick32, ...
                    xticklabel3, yticklabel32, ...
                    xtitle3, ytitle32);
                %----------------------------------------------------------
                plot(Fa2 (:, 1, 1), Fa2(:, 2, 1), M(1, 1), 'Color', Color(1, :), 'LineWidth', LW(1, 1));
                plot(th1W(:, 1)   , V1          , M(1, 2), 'Color', Color(1, :), 'LineWidth', LW(1, 1));
                %----------------------------------------------------------
                yyaxis right;
                ax = gca;
                ax.YAxis(1).Color = Color(1, :); ax.YAxis(1).Label.Color = 'k';
                ax.YAxis(2).Color = Color(2, :);
                %----------------------------------------------------------
                PlotX1(...
                    xlim3, ylim33, ...
                    xtick3, ytick33, ...
                    xticklabel3, yticklabel33, ...
                    xtitle3, "");
                %----------------------------------------------------------
                plot(Fa2 (:, 1, 2), Fa2(:, 2, 2), M(1, 1), 'Color', Color(2, :), 'LineWidth', LW(1, 1));
                plot(th1W(:, 2)   , V2          , M(1, 2), 'Color', Color(2, :), 'LineWidth', LW(1, 1));
                %----------------------------------------------------------
            otherwise
                return
        end
        %------------------------------------------------------------------
        P4{i, 1} = plot(nan, nan, '-', 'Color', Color(1, :), 'LineWidth', LW);
        P4{i, 2} = plot(nan, nan, '-', 'Color', Color(2, :), 'LineWidth', LW);
        P4{i, 3} = plot(nan, nan, ':', 'Color', 'k'        , 'LineWidth', LW);
        P4{i, 4} = plot(nan, nan, '-', 'Color', 'k'        , 'LineWidth', LW);
        PlotX2(...
            P4(i, :), legendlabel3, 'northeast', numcols3);
        %------------------------------------------------------------------
        if export(1, 1)
            exportgraphics(fig3, fid3(i, 1), 'ContentType', 'Vector', 'Resolution', 600);
            movefile      (fid3(i, 1), path13);
        end
        %------------------------------------------------------------------
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z] = za(G, h0, A)
syms x t;
c  = sqrt(G.*(h0+A));
k  = sqrt(3.*A)./(2.*h0.*sqrt(h0+A));
xs =-50;
z1 = A.*sech(k.*(x-xs-c.*t)).^2;
z2 = A.*sech(k.*(x+xs+c.*t)).^2;
z  = z1+z2;
Z  = matlabFunction(vpa(z, 16), 'Vars', [x, t]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%