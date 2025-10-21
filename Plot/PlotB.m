function [obj] = PlotB(g, obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h0   = g.data.h0;
test = g.test;
Za   = g.data.Z(obj.xc, g.t);
H_   = sum(pagemtimes(g.x(:, :, 1)-g.zb, g.Wbf'), 2)./sum(g.W, 2);
Z_   = sum(pagemtimes(g.x(:, :, 1)     , g.Wbf'), 2)./sum(g.W, 2);
if test ~= 13 && test ~= 16
    Za = Za./h0;
    Z_ = Z_./h0;
end
%--------------------------------------------------------------------------
if mod(g.nit, 1) == 0
    if test ~= 7
        if ismembc(test, [6, 8, 12, 13])
            set(obj.Ph{1, 1}, 'YData', Za);
            set(obj.Ph{1, 2}, 'YData', Z_);
        elseif ismembc(test, [1, 2, 3, 4, 5, 9, 11, 15, 16])
            set(obj.Ph{1, 2}, 'YData', Z_);
        elseif test == 10
            set(obj.Ph{1, 2}, 'YData', H_);
        end
    else
        PHIa = g.data.PHI(obj.xc, g.t);
        HYDa = g.data.HYD(obj.xc, g.t);
        HQ1a = g.data.HQ1(obj.xc, g.t);
        PHIh = sum(pagemtimes(g.PHIN  , g.Wbf'), 2)./sum(g.W, 2);
        HYDh = sum(pagemtimes(g.GHd1ZN, g.Wbf'), 2)./sum(g.W, 2);
        HQ1h = sum(pagemtimes(g.HQ1N  , g.Wbf'), 2)./sum(g.W, 2);
        set(obj.Ph{1, 1}, 'YData', Za);
        set(obj.Ph{1, 2}, 'YData', Z_);
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
obj.Write1(g.t, g.x, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if g.t > (obj.i-1)*obj.ndt || any(g.t == g.data.tk, 2)
    %----------------------------------------------------------------------
    if test == 1 && g.data.opt == 3
        phi = g.PHIN;
    else
        phi = nan;
    end
    obj.Write2(g.t, g.x, phi, 0);
    %----------------------------------------------------------------------
    if obj.i > obj.nf
        if ismembc(test, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15])
            tm    = obj.tm;
            tM    = obj.tM;
            vm    = obj.vm;
            vM    = obj.vM;
            Um    = obj.Um;
            UM    = obj.UM;
            time  = obj.time;
            data  = obj.data;
            datag = obj.datag;
            if obj.n ~= 0
                save(obj.fid, 'time', 'data', 'datag', 'tm', 'vm', 'Um', 'tM', 'vM', 'UM');
            else
                save(obj.fid, 'time', 'data', 'tm', 'vm', 'Um', 'tM', 'vM', 'UM');
            end
        elseif ismembc(test, [12, 13])
            time  = obj.time;
            data  = obj.data;
            wd    = obj.wd;
            save(obj.fid, 'time', 'data', 'wd');
        end
    end
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end