function [obj] = PlotA1(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Color  = linspecer(9, 'Qualitative');
grey   = repmat(0.80, 1, 3);
LW     = [1.0, 2.5];
%--------------------------------------------------------------------------
K      = g.numE;
p      = g.p;
h0     = g.data.h0;
flag   = 0; % Plot auxiliary vars.
test   = g.test;
wetdry = g.data.wetdry;
xv     = g.data.xv;
xc     = zeros(K, 1);
for i = 1:K
    xc(i, 1) = 1./2.*(xv(i, 1)+xv(i+1, 1));
end
aux    = PlotA2(g);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test == 11
    Xa = g.data.H(xc, 0);
    Xh = sum(pagemtimes(g.x(:, :, 1)-g.zbinit, g.Wbf'), 2)./sum(g.W, 2);
else
    h0 = g.data.h0;
    Xa = g.data.N(xc, 0)-h0;
    Xh = sum(pagemtimes(g.x(:, :, 1)-h0, g.Wbf'), 2)./sum(g.W, 2);
    Xa = Xa./h0;
    Xh = Xh./h0;
end
if test == 12
    xc = xc./h0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test == 2 && flag
    %----------------------------------------------------------------------
    ph = cell(4, 2);
    Ph = cell(4, 2);
    for i = 1:2
        figure('Color', 'w', 'Windowstate', 'maximized');
        for j = 1:2
            %--------------------------------------------------------------
            k = 2*(i-1)+j;
            %--------------------------------------------------------------
            subplot(1, 2, j);
            PlotX1(...
                aux.xlim, aux.ylim(k, :), ...
                aux.xtick, aux.ytick{k, 1}, ...
                aux.xticklabel, aux.yticklabel{k, 1}, ...
                xtitle, ytitle);
            if i ~= 1 || j ~= 1
                Xa = zeros(K, 1);
                Xh = zeros(K, 1);
            end
            ph{k, 1} = plot(nan, nan, 'Color', Color(1, :), 'LineWidth', LW(1, 1));
            Ph{k, 1} = plot(xc , Xa , 'Color', Color(1, :), 'LineWidth', LW(1, 1));
            ph{k, 2} = plot(nan, nan, 'Color', Color(2, :), 'LineWidth', LW(1, 2));
            Ph{k, 2} = plot(xc , Xh , 'Color', Color(2, :), 'LineWidth', LW(1, 2));
            PlotX2(...
                ph(k, :), aux.legendlabel(k, :), 'northeast', numcols);
            %--------------------------------------------------------------
            grid on;
            grid minor;
            %--------------------------------------------------------------
        end
    end
    %----------------------------------------------------------------------
else
    %----------------------------------------------------------------------
    ph = cell(1, 2);
    Ph = cell(1, 2);
    figure('Color', 'w', 'Renderer', 'painters', 'Windowstate', 'maximized');
    hold on;
    PlotX1(...
        aux.xlim, aux.ylim, ...
        aux.xtick, aux.ytick, ...
        aux.xticklabel, aux.yticklabel, ...
        aux.xtitle, aux.ytitle);
    if ~ismembc(test, [12, 13, 14, 15, 16, 17])
        ph{1, 1} = plot(nan, nan, 'Color', Color(1, :), 'LineWidth', LW(1, 1));
        Ph{1, 1} = plot(xc , Xa , 'Color', Color(1, :), 'LineWidth', LW(1, 1));
    end
    ph{1, 2} = plot(nan, nan, 'Color', Color(2, :), 'LineWidth', LW(1, 2));
    Ph{1, 2} = plot(xc , Xh , 'Color', Color(2, :), 'LineWidth', LW(1, 2));
    if ismembc(test, [1, 2, 7, 8, 9])
        PlotX2(...
            ph(1, :), aux.legendlabel, 'northeast', aux.numcols);
    end
    %----------------------------------------------------------------------
end
if test == 1 || wetdry
    xbathy = sum(g.coordV0T, 2)./2;
    ybathy = meanval(g, g.zbinit);
    patch(...
        xbathy, ybathy, grey, 'EdgeColor', grey);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = sprintf('%s_P%d_%d.mat', aux.path, p, K);
if isfile(fid)
    delete(fid);
end
obj = Output(g, fid, Ph);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end