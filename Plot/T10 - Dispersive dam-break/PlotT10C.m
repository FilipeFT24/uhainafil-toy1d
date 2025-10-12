function [] = PlotT10C(g, export)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color = linspecer(9, 'Qualitative');
M     = [":", "-", "-.", "--"];
LW    = 2.5;
flag  = [1, 0];
path1 = "Plot/T10 - Dispersive dam-break/DATA";
path2 = "Plot/T10 - Dispersive dam-break/";
%--------------------------------------------------------------------------
h0    = g.data.h0;
G     = g.data.G;
A     = g.data.A;
xm    =-300;
xM    = 300;
K     = 6000;
xv    = linspace(xm, xM, K+1)';
xc    = zeros(K, 1);
for i = 1:K
    xc(i, 1) = 1./2.*(xv(i, 1)+xv(i+1, 1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m   = 2;
n   = 3;
nf  = g.data.nf+numel(g.data.tk);
P1  = cell (m , 3+n);
P2  = cell (m , 3);
Fh1 = zeros(nf, K, n, m);
for i = 1:n
    [Fh1(:, :, i, 1), Fh1(:, :, i, 2)] = geth1(sprintf("%s/Ddb_P%d_01_GN.mat", path1, i), i, h0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid1          = ["Ddb_h.pdf"; "Ddb_u.pdf"];
%--------------------------------------------------------------------------
xlim          = [-300, 300];
xtick         =-300:100:300;
xticklabel    = ["-300", "-200", "-100", "0", "100", "200", "300"];
xtitle        = "$x$";
%--------------------------------------------------------------------------
ylim1         = [0.95, 2.00];
ytick1        = [0.95, 1.00:0.20:2.00];
yticklabel1   = ["\phantom{-}0.95", "1.00", "1.20", "1.40", "1.60", "1.80", "2.00"];
ytitle1       = "$h$";
legendlabel1  = ["$h_{0}$", "$h_{\phantom{0}}^{\ast}$", "$h_{R}+a^{+}$", "$k=1$", "$k=2$", "$k=3$"];
ylim2         = [-0.10, 2.25];
ytick2        = [-0.10, 0.00:0.50:2.00, 2.25];
yticklabel2   = ["-0.10",  "0.00", "0.50", "1.00", "1.50", "2.00", "2.25"];
ytitle2       = "$u$";
legendlabel2  = ["$u_{0}$", "$u_{\phantom{0}}^{\ast}$", "", "$k=1$", "$k=2$", "$k=3$"];
%--------------------------------------------------------------------------
xlim1z        = [45, 55];
xtick1z       = 45:5:55;
xticklabel1z  = ["45", "50", "55"];
xlim2z        = [170, 200];
xtick2z       = 170:10:200;
xticklabel2z  = ["170", "180", "190", "200"];
%--------------------------------------------------------------------------
ylim11z       = [1.35, 1.39];
ytick11z      = 1.35:0.01:1.39;
yticklabel11z = ["1.35", "1.36", "1.37", "1.38", "1.39"];
ylim12z       = [1.00, 1.80];
ytick12z      = 1.00:0.20:1.80;
yticklabel12z = ["1.00", "1.20", "1.40", "1.60", "1.80"];
ylim21z       = [1.04, 1.10];
ytick21z      = 1.04:0.02:1.10;
yticklabel21z = ["1.04", "1.06", "1.08", "1.10"];
ylim22z       = [0.00, 1.80];
ytick22z      = 0.00:0.60:1.80;
yticklabel22z = ["0.00", "0.60", "1.20", "1.80"];
%--------------------------------------------------------------------------
ylim          = {ylim1; ylim2};
ytick         = {ytick1; ytick2};
yticklabel    = {yticklabel1; yticklabel2};
ytitle        = {ytitle1; ytitle2};
legendlabel   = {legendlabel1; legendlabel2};
%--------------------------------------------------------------------------
xlimz         = {xlim1z, xlim2z; xlim1z, xlim2z};
xtickz        = {xtick1z, xtick2z; xtick1z, xtick2z};
xticklabelz   = {xticklabel1z, xticklabel2z; xticklabel1z, xticklabel2z};
ylimz         = {ylim11z, ylim12z; ylim21z, ylim22z};
ytickz        = {ytick11z, ytick12z; ytick21z, ytick22z};
yticklabelz   = {yticklabel11z, yticklabel12z; yticklabel21z, yticklabel22z};
%--------------------------------------------------------------------------
zf            = 1.025./2;
xz1_          = xlim1z.*zf;
xz2_          = xlim2z.*zf;
xr11          = [-250.0,  80.0].*zf;
xr12          = [-112.5, 112.5].*zf;
xr21          = [-112.5, 130.0].*zf;
xr22          = [-250.0,  80.0].*zf;
yz11          = ylim11z;
yz12          = ylim12z;
yr11          = [1.04, 1.32];
yr12          = [1.50, 1.78];
yz21          = ylim21z;
yz22          = ylim22z;
yr21          = [0.0875, 0.0875+(2.25./1.05).*0.28];
yr22          = [1.2500, 1.2500+(2.25./1.05).*0.35];
%--------------------------------------------------------------------------
xz            = {xz1_, xz2_; xz1_, xz2_};
xr            = {xr11, xr12; xr21, xr22};
yz            = {yz11, yz12; yz21, yz22};
yr            = {yr11, yr12; yr21, yr22};
%--------------------------------------------------------------------------
numcols       = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 1)
    for i = 1:m
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fig1 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        cp = get(gca, 'Position');
        hold on;
        %------------------------------------------------------------------
        P1{i, 1} = plot(nan, nan, M(1, 1), 'Color', Color(1, :), 'LineWidth', LW);
        P1{i, 2} = plot(nan, nan, M(1, 1), 'Color', Color(9, :), 'LineWidth', LW);
        switch i
            case 1
                P1{i, 3} = plot(nan, nan, M(1, 1), 'Color', Color(5, :), 'LineWidth', LW);
            case 2
                P1{i, 3} = plot(nan, nan, M(1, 1), 'Color', 'w', 'LineWidth', LW);
            otherwise
                return
        end
        for j = 4:6
            P1{i, j} = plot(nan, nan, M(1, j-2), 'Color', Color(j-2, :), 'LineWidth', LW);
        end
        %------------------------------------------------------------------
        switch i
            case 1
                fplot(@(x) g.data.H(x, 0), M(1, 1), 'Color', Color(1, :), 'LineWidth', LW);
            case 2
                fplot(@(x) g.data.U(x, 0), M(1, 1), 'Color', Color(1, :), 'LineWidth', LW);
            otherwise
                return
        end
        for j = 1:n
            plot(xc, Fh1(end, :, j, i)', M(1, 1+j), 'Color', Color(1+j, :), 'LineWidth', LW);
        end
        %------------------------------------------------------------------
        switch i
            case 1
                H1 = (sqrt(h0)+sqrt(h0+A)).^2./4;
                H2 = h0+A-1./12.*A.^2;
                plot([ 51.05, 300.00], [H1, H1], M(1, 1), 'Color', Color(9, :), 'LineWidth', LW);
                plot([150.00, 300.00], [H2, H2], M(1, 1), 'Color', Color(5, :), 'LineWidth', LW);
            case 2
                H1 = (sqrt(h0)+sqrt(h0+A)).^2./4;
                U1 = 2.*(sqrt(G.*H1)-sqrt(G.*1));
                plot([ 51.05, 300.00], [U1, U1], M(1, 1), 'Color', Color(9, :), 'LineWidth', LW);
            otherwise
                return
        end
        %------------------------------------------------------------------
        PlotX1(...
            xlim, ylim{i, 1}, ...
            xtick, ytick{i, 1}, ...
            xticklabel, yticklabel{i, 1}, ...
            xtitle, ytitle{i, 1});
        PlotX2(...
            P1(i, :), legendlabel{i, 1}, 'northEast', numcols);
        %------------------------------------------------------------------
        MagInset(gcf, -1, ...
            [xz{i, 1}, yz{i, 1}], ...
            [xr{i, 1}, yr{i, 1}], {'NW', 'NW'; 'NE', 'NE'}, cp, xlimz{i, 1}, ylimz{i, 1}, xtickz{i, 1}, ytickz{i, 1}, xticklabelz{i, 1}, yticklabelz{i, 1}, 10.5);
        MagInset(gcf, -1, ...
            [xz{i, 2}, yz{i, 2}], ...
            [xr{i, 2}, yr{i, 2}], {'SE', 'SE'; 'NE', 'NE'}, cp, xlimz{i, 2}, ylimz{i, 2}, xtickz{i, 2}, ytickz{i, 2}, xticklabelz{i, 2}, yticklabelz{i, 2}, 10.5);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if export(1, 1)
            exportgraphics(fig1, fid1(i, 1), 'ContentType', 'Vector', 'Resolution', 600);
            movefile      (fid1(i, 1), path2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag(1, 2)
    for i = 1:m
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~export(1, 2)
            fig2 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
        else
            fig2 = figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized', 'MenuBar', 'none', 'Toolbar', 'none');
        end
        hold on;
        %------------------------------------------------------------------
        %
        switch i
            case 1
                P2{i, 1} = fplot(@(x) g.data.H(x, 0), M(1, 1), 'Color', Color(1, :), 'LineWidth', LW);
            case 2
                P2{i, 1} = fplot(@(x) g.data.U(x, 0), M(1, 1), 'Color', Color(1, :), 'LineWidth', LW);
            otherwise
                return
        end
        %
        P3 = plot(xc, Fh1(1, :, 3, i)', '-', 'Color', Color(2, :), 'LineWidth', LW);
        %
        switch i
            case 1
                H1       = (sqrt(h0)+sqrt(h0+A)).^2./4;
                H2       = h0+A-1./12.*A.^2;
                P2{i, 2} = plot([ 51.05, 300.00], [H1, H1], M(1, 1), 'Color', Color(9, :), 'LineWidth', LW);
                P2{i, 3} = plot([150.00, 300.00], [H2, H2], M(1, 1), 'Color', Color(5, :), 'LineWidth', LW);
                PlotX2(...
                    P2(i, 1:3), legendlabel{i, 1}(1, 1:3), 'northEast', 1);
            case 2
                H1       = (sqrt(h0)+sqrt(h0+A)).^2./4;
                U1       = 2.*(sqrt(G.*H1)-sqrt(G.*1));
                P2{i, 2} = plot([ 51.05, 300.00], [U1, U1], M(1, 1), 'Color', Color(9, :), 'LineWidth', LW);
                PlotX2(...
                    P2(i, 1:2), legendlabel{i, 1}(1, 1:2), 'northEast', 1);
            otherwise
                return
        end
        %
        PlotX1(...
            xlim, ylim{i, 1}, ...
            xtick, ytick{i, 1}, ...
            xticklabel, yticklabel{i, 1}, ...
            xtitle, ytitle{i, 1});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        video.src      = sprintf("%s/0.20/VIDEO.mat", path1);
        switch i
            case 1
                path_mp4 = sprintf("%s/VIDEO2_P3 (h)", path2);
                path_png = "Plot/T10 - Dispersive dam-break/FRAMES/h (2)";
            case 2
                path_mp4 = sprintf("%s/VIDEO2_P3 (u)", path2);
                path_png = "Plot/T10 - Dispersive dam-break/FRAMES/u (2)";
            otherwise
                return
        end
        video.path_mp4 = path_mp4;
        video.path_png = path_png;
        video.fig      = fig2;
        video.P1       = P3;
        video.Fh1      = Fh1(:, :, 3, i);
        video.nf       = nf;
        %------------------------------------------------------------------
        WriteVideo1(video, export(1, 2));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fh1, Fh2] = geth1(fid, p, h0)
Fh1 = getdata1(fid, p, 1)+h0;
Fh2 = getdata1(fid, p, 2);
Fh2 = Fh2./Fh1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%