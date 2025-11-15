function [obj] = PlotB(g, obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xc   = obj.xc;
t    = g.t;
tend = g.data.tend;
test = g.test;
if test == 11
    Xa = g.data.H(xc, t);
    Xh = sum(pagemtimes(g.x(:, :, 1)-g.zbinit, g.Wbf'), 2)./sum(g.W, 2);
else
    Xa = g.data.Z(xc, t);
    Xh = sum(pagemtimes(g.x(:, :, 1), g.Wbf'), 2)./sum(g.W, 2);
    h0 = g.data.h0;
    Xa = Xa./h0;
    Xh = Xh./h0;
end
%--------------------------------------------------------------------------
if mod(g.nit, 25) == 0
    if size(obj.Ph, 1) == 1
        if ismembc(test, [1, 2, 6, 7, 8, 9, 10])
            set(obj.Ph{1, 1}, 'YData', Xa);
            set(obj.Ph{1, 2}, 'YData', Xh);
        elseif ismembc(test, [3, 4, 5, 11, 12, 13, 14, 15, 16, 17])
            set(obj.Ph{1, 2}, 'YData', Xh);
        end
    else
        PHIa = g.data.PHI(xc, t);
        HYDa = g.data.d1Z(xc, t);
        HQ1a = g.data.HQ1(xc, g.t);
        PHIh = sum(pagemtimes(g.PHIN   , g.Wbf'), 2)./sum(g.W, 2);
        HYDh = sum(pagemtimes(g.GHd1ZNX, g.Wbf'), 2)./sum(g.W, 2);
        HQ1h = sum(pagemtimes(g.HQ1N   , g.Wbf'), 2)./sum(g.W, 2);
        set(obj.Ph{1, 1}, 'YData', Xa);
        set(obj.Ph{1, 2}, 'YData', Xh);
        set(obj.Ph{2, 1}, 'YData', PHIa);
        set(obj.Ph{2, 2}, 'YData', PHIh);
        set(obj.Ph{3, 1}, 'YData', HYDa);
        set(obj.Ph{3, 2}, 'YData', HYDh);
        set(obj.Ph{4, 1}, 'YData', HQ1a);
        set(obj.Ph{4, 2}, 'YData', HQ1h);
    end
    drawnow limitrate;
end
%--------------------------------------------------------------------------
if ismembc(test, [1, 2])
    if g.t+eps >= tend
        %------------------------------------------------------------------
        K  = g.numE;
        N  = g.N;
        D  = 1;
        V  = 1+D;
        bf = g.bf;
        if test == 2
            V = V+2;
        end
        ua = zeros(K, N, V);
        uh = zeros(K, N, V);
        ec = zeros(K, V);
        e2 = zeros(1, V);
        %------------------------------------------------------------------
        for i = 1:1+D
            uh(:, :, i) = g.x(:, :, i);
        end
        ua(:, :, 1) = inittype(g.itype, @(x) g.data.Z (x, g.t), g.xydc, g.xyqc, g.fc);
        ua(:, :, 2) = inittype(g.itype, @(x) g.data.HU(x, g.t), g.xydc, g.xyqc, g.fc);
        if test == 2
            uh(:, :, 3) = g.GHd1ZN;
            uh(:, :, 4) = g.HQ1N;
            ua(:, :, 3) = inittype(g.itype, @(x) g.data.HYD(x, g.t), g.xydc, g.xyqc, g.fc);
            ua(:, :, 4) = inittype(g.itype, @(x) g.data.HQ1(x, g.t), g.xydc, g.xyqc, g.fc);
        end
        for i = 1:V
            eN       = ua(:, :, i)-uh(:, :, i);
            eq       = eN*bf';
            ec(:, i) = eq*g.W';
            e2(1, i) = sqrt(sum(((eq.^2)*g.W').*g.detJ0T, 1));
        end
        save(obj.fid, 'ua', 'uh', 'ec');
    end
    %----------------------------------------------------------------------
else
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
end




% obj.Write1(g.t, g.x, 0);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if g.t > (obj.i-1)*obj.ndt || any(g.t == g.data.tk, 2)
%     %----------------------------------------------------------------------
%     obj.Write2(g.t, g.x, nan, 0);
%     %----------------------------------------------------------------------
%     if obj.i > obj.nf
%         if ismembc(test, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15])
%             tm    = obj.tm;
%             tM    = obj.tM;
%             vm    = obj.vm;
%             vM    = obj.vM;
%             Um    = obj.Um;
%             UM    = obj.UM;
%             time  = obj.time;
%             data  = obj.data;
%             datag = obj.datag;
%             if obj.n ~= 0
%                 save(obj.fid, 'time', 'data', 'datag', 'tm', 'vm', 'Um', 'tM', 'vM', 'UM');
%             else
%                 save(obj.fid, 'time', 'data', 'tm', 'vm', 'Um', 'tM', 'vM', 'UM');
%             end
%         elseif ismembc(test, [12, 13])
%             time  = obj.time;
%             data  = obj.data;
%             wd    = obj.wd;
%             save(obj.fid, 'time', 'data', 'wd');
%         end
%     end
%     %----------------------------------------------------------------------
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end