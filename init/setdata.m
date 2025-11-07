function [data] = setdata(test, flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x t;
switch test
    case 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S0/T0 - Convergence (Green-Nagdhi)
        %------------------------------------------------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S0/T1 - Convergence (SW)
        %------------------------------------------------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S0/T2 - Wave generation
        %------------------------------------------------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S1/T1 - Lake at rest
        %------------------------------------------------------------------
        G        = 1;
        h0       = 0.25;
        nm       = 0;
        wetdry   = 0;
        option   = 0;
        %------------------------------------------------------------------
        xm       =-3;
        xM       = 3;
        K        = 1000;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        switch option
            case 0 % EXP
                zb1 =-0.50.*exp(-50.*(x+1).^2); % WET
                zb2 = 0.50.*exp(-25.*(x-1).^2); % DRY
            case 1 % LINEAR
                zb1 =-0.50.*((x+1.50).*heaviside(x+1.50).*heaviside(-(1.00+x))+(-0.50-x).*heaviside(-0.50-x).*heaviside(x+1.00));
                zb2 = 0.50.*((x-0.50).*heaviside(x-0.50).*heaviside( (1.00-x))+( 1.50-x).*heaviside( 1.50-x).*heaviside(x-1.00));
            case 2 % QUADRATIC
                zb1 =-0.50.*(1-((x+1)./0.5).^2).*heaviside(x+1.5).*heaviside(-0.5-x);
                zb2 = 0.50.*(1-((x-1)./0.5).^2).*heaviside(x-0.5).*heaviside( 1.5-x);
            otherwise
                return
        end
        zb       = zb1+zb2;
        z        = h0+(zb-h0).*(0.5.*(sign(zb2.*heaviside(x)-h0)+1));
        h        = z-zb;
        u        = 0;
        hu       = h.*u;
        tend     = 10e3;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 5
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S1/T2 - Drying of a lake
        %------------------------------------------------------------------
        G        = 1;
        h0       = 0.60;
        nm       = 0;
        wetdry   = 1;
        option   = 1;
        %------------------------------------------------------------------
        xm       =-0.50;
        xM       = 0.50;
        K        = 200;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        a        = 0.25;
        b        = 0.10;
        k        = pi./b;
        zb       = a.*(cos(k.*x)+1).*heaviside(x+b).*heaviside(b-x);
        z        = h0+(zb-h0).*heaviside(zb-h0);
        h        = z-zb;
        u        = 0;
        hu       = h.*u;
        tend     = 1000;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 6
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S1/T3 - Oscillating lake
        %------------------------------------------------------------------
        G        = 1;
        h0       = 0.1005.*G;
        nm       = 0;
        wetdry   = 1;
        %------------------------------------------------------------------
        xm       =-1.5;
        xM       = 1.5;
        K        = 3000;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        a0       = 0.1;
        omega    = sqrt(2.*h0);
        zb       = h0.*x.^2;
        h_aux    = h0+2.*h0.*a0.*cos(omega.*t).*(x-a0./2.*cos(omega.*t))-zb;
        h        = h_aux.*heaviside(h_aux);
        z        = h+zb;
        u        =-a0.*omega.*sin(omega.*t);
        hu       = h.*u;
        tend     = 2.*pi./omega;
        tk       = [0.25, 0.50, 0.75, 1.00].*tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 7
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S2/T1 - Head-on collision
        %------------------------------------------------------------------
        alpha    = 1;
        G        = 1;
        h0       = 1;
        nm       = 0;
        wetdry   = 0;
        option   = 1;
        %------------------------------------------------------------------
        xm       =-100;
        xM       = 100;
        K        = 2000;
        %{
        dxmin    = 2e-02;
        func     = @(y) (L./(2.*N)).*y./sinh(y./2)-dxmin;
        beta     = fzero(func, [1, 25]);
        xi       = linspace(0, 1, N+1)';
        xv       = xm+L./2.*(1+sinh(beta.*(xi-1./2))./sinh(beta./2));
        %}
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        switch option
            case 1
                A1 = 0.96.*h0;
                A2 = A1;
            case 2
                h0 = 1;
                A1 = 0.21.*h0;
                A2 = A1;
            otherwise
                return
        end
        c1       = sqrt(G.*(h0+A1));
        c2       = sqrt(G.*(h0+A2));
        k1       = sqrt(3.*A1)./(2.*h0.*sqrt(h0+A1));
        k2       = sqrt(3.*A2)./(2.*h0.*sqrt(h0+A2));
        xs1      =-50;
        xs2      = 50;
        zb       =-h0;
        z1       = A1.*sech(k1.*(x-xs1-c1.*t)).^2;
        z2       = A2.*sech(k2.*(x-xs2+c2.*t)).^2;
        z        = z1+z2;
        h        = z-zb;
        u1       = c1.*(1-h0./(h0+z1));
        u2       =-c2.*(1-h0./(h0+z2));
        u        = u1+u2;
        hu       = h.*u;
        tend     = 100./c2;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 8
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S2/T2 - Wall
        %------------------------------------------------------------------
        alpha    = 1;
        G        = 1;
        h0       = 1;
        nm       = 0;
        wetdry   = 0;
        option   = 2;
        %------------------------------------------------------------------
        xm       =-100;
        xM       = 0;
        K        = 1000;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        switch option
            case 1
                A1 = 0.96.*h0;
            case 2
                A1 = 0.21.*h0;
            otherwise
                return
        end
        c        = sqrt(G.*(h0+A1));
        k        = sqrt(3.*A1)./(2.*h0.*sqrt(h0+A1));
        xs       =-100;
        zb       =-h0;
        z        = A1.*sech(k.*(x-xs-c.*t)).^2;
        h        = z-zb;
        u        = c.*(1-h0./(h0+z));
        hu       = h.*u;
        tend     = 100./c;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 9
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S2/T3 - Overtaking collision
        %------------------------------------------------------------------
        alpha    = 1;
        G        = 1;
        h0       = 1;
        nm       = 0;
        wetdry   = 0;
        %------------------------------------------------------------------
        xm       =-200;
        xM       = 200;
        K        = 4000;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        A1       = 0.96.*h0;
        A2       = 0.44.*h0;
        c1       = sqrt(G.*(h0+A1));
        c2       = sqrt(G.*(h0+A2));
        k1       = sqrt(3.*A1)./(2.*h0.*sqrt(h0+A1));
        k2       = sqrt(3.*A2)./(2.*h0.*sqrt(h0+A2));
        xs2      =-150;
        xs1      = xs2.*c1./c2;
        zb       =-h0;
        z1       = A1.*sech(k1.*(x-xs1-c1.*t)).^2;
        z2       = A2.*sech(k2.*(x-xs2-c2.*t)).^2;
        z        = z1+z2;
        h        = z-zb;
        u1       = c1.*(1-h0./(h0+z1));
        u2       = c2.*(1-h0./(h0+z2));
        u        = u1+u2;
        hu       = h.*u;
        tend     =-2.*xs1./c1;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 10
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S2/T4 - Breakup of a Gaussian hump
        %------------------------------------------------------------------
        alpha    = 1;
        G        = 1;
        h0       = 1;
        nm       = 0;
        wetdry   = 0;
        option   = 2;
        %------------------------------------------------------------------
        xm       =-100;
        xM       = 100;
        K        = 2000;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        switch option
            case 1
                A1 = 1;
            case 2
                A1 = 20;
            otherwise
                return
        end
        zb       =-h0;
        z        = exp(-x.^2./A1);
        h        = z-zb;
        u        = 0;
        hu       = h.*u;
        tend     = 100;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 11
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S2/T5 - Dispersive dam-break
        %------------------------------------------------------------------
        alpha    = 1;
        G        = 9.81;
        h0       = 1;
        nm       = 0;
        wetdry   = 0;
        %------------------------------------------------------------------
        xm       =-300;
        xM       = 300;
        K        = 6000;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        A1       = 0.8;
        zb       =-h0;
        z        = 1./2.*A1.*(1-tanh(x./0.4));
        h        = z-zb;
        u        = eps.*x;
        hu       = h.*u;
        tend     = 47.434;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 12
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S3/T1 - Grilli
        %------------------------------------------------------------------
        alpha    = 1;
        %        = 1.159;
        G        = 9.81;
        h0       = 0.44;
        nm       = 0;
        wetdry   = 0;
        option   = 3;
        %
        bathy    = 0;
        slope    = 1./34.70;
        %------------------------------------------------------------------
        xm       =-52.32.*h0;
        xM       = h0./slope;
        K        = 383;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        switch option
            case 1
                A1   = 0.10.*h0;
                tend = 52.52.*(sqrt(h0./G));
            case 2
                A1   = 0.15.*h0;
                tend = 48.52.*(sqrt(h0./G));
            case 3 % THIS
                A1   = 0.20.*h0;
                tend = 44.52.*(sqrt(h0./G));
            case 4
                A1   = 0.25.*h0;
                tend = 42.52.*(sqrt(h0./G));
            case 5
                A1   = 0.30.*h0;
                tend = 40.52.*(sqrt(h0./G));
            case 6
                A1   = 0.40.*h0;
                tend = 38.52.*(sqrt(h0./G));
            otherwise
                return
        end
        c        = sqrt(G.*(h0+A1));
        k        = sqrt(3.*A1)./(2.*h0.*sqrt(h0+A1));
        xs       =-5.2.*h0-(c.*13.56.*(sqrt(h0./G)));
        ds       = 0.5;
        switch bathy
            case 0
                zb =-h0+heaviside(x).*slope.*x;
            case 1
                f1 = bathyC1(slope, ds, 1);
                zb =-h0+...
                    (heaviside(x+ds)-heaviside(x-ds)).*f1(x)+...
                    (heaviside(x-ds)).*slope.*x;
            case 2
                f2 = bathyC2(slope, ds, 1);
                zb =-h0+...
                    (heaviside(x+ds)-heaviside(x-ds)).*f2(x)+...
                    (heaviside(x-ds)).*slope.*x;
            case 3
                f3 = bathyC3(slope, ds, 1);
                zb =-h0+...
                    (heaviside(x+ds)-heaviside(x-ds)).*f3(x)+...
                    (heaviside(x-ds)).*slope.*x;
            otherwise
                return
        end
        z        = A1.*sech(k.*(x-xs-c.*t)).^2;
        h        = z-zb;
        u        = c.*(1-h0./(h0+z));
        hu       = h.*u;
        xg       = [-5.00, 20.96, 22.55, 23.68, 24.68, 25.91]'.*h0;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 13
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S3/T2 - Reversibility check (shelf)
        %------------------------------------------------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 14
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S3/T3 - Reflection of shoaling waves (M. Walkley)
        %------------------------------------------------------------------
        alpha    = 1;
        %        = 1.159;
        G        = 9.81;
        h0       = 0.7;
        nm       = 0;
        wetdry   = 0;
        option   = 1;
        %
        bathy    = 0;
        slope    = 1./50;
        %------------------------------------------------------------------
        xm       =-55;
        xM       = 20;
        K        = 750;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        switch option
            case 1
                A1 = 0.1000.*h0;
            case 2
                A1 = 0.1834.*h0;
                %  = 0.1710.*h0;
            otherwise
                return
        end
        c        = sqrt(G.*(h0+A1));
        k        = sqrt(3.*A1)./(2.*h0.*sqrt(h0+A1));
        xs       =-30;
        z        = A1.*sech(k.*(x-xs-c.*t)).^2;
        zb       =-h0+heaviside(x).*slope.*(x);
        h        = z-zb;
        u        = c.*(1-h0./(h0+z));
        hu       = h.*u;
        xg       = [0.00, 16.25, 17.75]';
        tend     = 30;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 15
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S3/T4 - Synolakis et al.
        %------------------------------------------------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 16
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S3/T5 - Submerged bar
        %------------------------------------------------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 17
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S3/T6 - Wave overtopping over a seawall
        %------------------------------------------------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
        return
end





switch test
    case 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 0) T00 - Convergence
        % NOTE: No need to make 3 steps for the RK3 scheme (dt = constant!)
        %------------------------------------------------------------------
    case 1
    case 3
    case 4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4) T04 - Reflection of shoaling waves

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 5
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 5) T05 - Submerged bar
        %---------------------------------------------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 6
    case 7
    case 8
    case 9

    case 10
    case 11
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 11) T11 - Reversibility check
        %------------------------------------------------------------------
        wetdry   = 0;
        h0       = 1;
        bathy    = 3;
        slope    = 1./20;
        nm       = 0;
        %------------------------------------------------------------------
        abslayer = 0;
        alpha    = 1;
        G        = 9.81;
        A1       = 0.20.*h0;
        xm       = 0;
        xM       = 240;
        K        = 2400;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        c        = sqrt(G.*(h0+A1));
        k        = sqrt(3.*A1)./(2.*h0.*sqrt(h0+A1));
        xs       = 80;
        z        = A1.*sech(k.*(x-xs-c.*t)).^2;
        switch bathy
            case 0
                zb =-h0+...
                    slope.*(x-130).*(heaviside(x-130)-heaviside(x-140))+...
                    slope.*10.*heaviside(x-140);
            case 3
                ds = 0.5;
                f3 = bathyC3(slope, ds, 1);
                g3 = bathyC3(slope, ds, 3);
                zb =-h0+...
                    (heaviside(x-(130-ds))-heaviside(x-(130+ds))).*f3(x-130)+...
                    (heaviside(x-(130+ds))-heaviside(x-(140-ds))).*slope.*(x-130)+...
                    (heaviside(x-(140-ds))-heaviside(x-(140+ds))).*g3(x-140)+...
                    slope.*10.*heaviside(x-(140-ds));
            otherwise
                return
        end
        h        = z-zb;
        u        = c.*(1-h0./(h0+z));
        hu       = h.*u;
        tend     = 50;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 12
    case 13
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 15
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 15) T15 - Synolakis
        %------------------------------------------------------------------
        %
        wetdry   = 1;
        h0       = 1;
        bathy    = 0;
        slope    = 1./19.85;
        nm       = 0;
        option   = 2;
        %------------------------------------------------------------------
        abslayer = 0;
        alpha    = 1;
        %        = 1.159;
        G        = 9.81;
        sqrth0_G = sqrt(h0./G);
        switch option
            case 1
                A1   = 0.019.*h0;
                tend = 70.*sqrth0_G;
                tk   = (25:5:70).*sqrth0_G;
            case 2
                A1   = 0.040.*h0;
                tend = 100.*sqrth0_G;
                tk   = (20:6:62).*sqrth0_G;
            otherwise
                return
        end
        gamma    = sqrt(3.*A1./(4.*h0));
        xs       = 1./slope+1./gamma.*acosh(sqrt(1./0.05));
        xm       =-20;
        xM       = 80;
        K        = 300;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        ds       = 0.5;
        switch bathy
            case 0
                zb =-h0+...
                    heaviside(-x+1./slope).*slope.*(-x+1./slope);
            case 3
                f3 = bathyC3(slope, ds, 2);
                zb =-h0+...
                    slope.*(-x+1./slope).*heaviside(-x+1./slope-ds)+...
                    f3(x-1./slope).*(heaviside(x-(1./slope-ds))-heaviside(x-(1./slope+ds)));
            otherwise
                return
        end
        c        =-sqrt(G.*(h0+A1));
        k        = sqrt(3.*A1)./(2.*h0.*sqrt(h0+A1));
        z        = A1.*sech(k.*(x-xs-c.*t)).^2;
        z        = z.*heaviside(x-eps)+heaviside(-(x+eps)).*zb;
        h        = z-zb;
        u        = c.*(1-h0./(h0+z));
        hu       = h.*u;
        %{
        %------------------------------------------------------------------
        h0       = 1;
        hi       = 10;
        option   = 1;
        wetdry   = 1;
        nm       = 0;
        %------------------------------------------------------------------
        abslayer = 0;
        alpha    = 1;
        G        = 10;
        xm       = 0;
        xM       = 12;
        K        = 1000;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        zb       = hi.*heaviside(x-10);
        hu       = hi.*heaviside(x+10)+(  -hi).*(heaviside(x-10)-heaviside(x-(12.5)));
        z        = h0.*heaviside(x+10)+(hi-h0).*(heaviside(x-10)-heaviside(x-(12.5)));
        h        = z-zb;
        u        = hu./h;
        tend     = 100;
        tk       = tend;
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 16

    case 17
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 17) T17 - Wave generation
        %------------------------------------------------------------------
        wetdry   = 0;
        h0       = 0.5;
        nm       = 0;
        %------------------------------------------------------------------
        abslayer = 0;
        alpha    = 1;
        G        = 10;
        sqrth0_G = sqrt(h0./G);
        xm       = 0;
        xM       = 400;
        K        = 2000;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        % WAVE GEN.
        a1       = -(alpha-1)./3;
        a        = -alpha./3;
        T        = 4;
        omega    = 2.*pi./T;
        wd       = omega.*h0;
        wd2      = wd.^2;
        k        = roots([G.*h0.^3.*a1 0 -a.*wd2-G.*h0 0 omega.^2]);
        %        = omega.*sqrt(1./(G.*h0-1./3.*wd2));
        kd       = k(1, 1).*h0;
        lambda   = 2.*pi./k(1, 1); %#ok<NASGU> 
        A        = 0.1.*h0;
        %------------------------------------------------------------------
        zb       =-h0;
        c        = sqrt(G.*(h0+A));
        z        = A.*sech(k(1, 1).*(x-0-c.*t)).^2;
        h        = z-zb;
        
        hu       = c.*(1-h0./(h0+z));
        u        = hu./h;


        ur1      = @(x, t) A.*sech(k(1, 1).*(x+0-c.*t)).^2;%A.*cos(omega.*t-k(1, 1).*x);
        ur2      = @(x, t) c.*(1-h0./(h0+ur1(x, t)));%h0.*ur1(x, t).*(omega./(kd));
        data.Ur1 = ur1;
        data.Ur2 = ur2;
        %{
        figure;
        hold on;
        plot(xv, ur1(xv, -5, 0), '-b');
        plot(xv, ur2(xv, -5, 0), '-r');
        %}
        tend     = 1000.*sqrth0_G;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
otherwise
    return;
end
%--------------------------------------------------------------------------
if test ~= 12
    data.Z    = mymatlabFunction(vpa(z , 16), [x, t]); % sol1
    data.HU   = mymatlabFunction(vpa(hu, 16), [x, t]); % sol2
    data.H    = mymatlabFunction(vpa(h , 16), [x, t]);
    data.U    = mymatlabFunction(vpa(u , 16), [x, t]);
end
data.Zb       = mymatlabFunction(vpa(zb, 16), [x, t]);
%--------------------------------------------------------------------------
data.wetdry   = wetdry;
data.abslayer = abslayer;
data.nm       = nm;
data.G        = G;
data.nf       = 1800;
data.tend     = tend;
data.tk       = tk;
data.xv       = xv;
data.dx       = dx;
if ismembc(test, [1, 2, 3, 4, 5, 6, 7, 9, 11, 12, 15, 16])
    if ismembc(test, [1, 11])
        data.bathy = bathy;
    end
    if ismembc(test, [1, 2, 3, 4, 5])
        data.xg = xg;
    end
    if test ~= 11
        data.opt = option;
    end
end
if ismembc(test, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 17])
    data.alpha = alpha;
end
if test == 10
    data.h0       = h0;
    data.A        = A1;
elseif test == 12
    data.A        = A;
    data.x0       = x0;
    data.L        = L;
    if option == 2
        data.e    = 0.1;
    end
    data.aL       = slope.*L;
    data.sqrtGaL  = sqrt(G.*slope.*L);
    data.sqrtL_aG = sqrt(L./(slope.*L));
else
    data.h0       = h0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismembc(test, [0, 1, 6, 7, 8, 14, 15, 17])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if test == 3 || test == 7
        z1 = z;
        z2 = A1.*sech(k.*(x-(2.*xM-xs)+c.*t)).^2;
        z  = z1+z2;
        if test == 3
            zb = zb+data.Zb(2.*xM-x, t);
        end
        h  = z-zb;
        u1 = u;
        u2 =-c.*(1-h0./(h0+z2));
        u  = u1+u2;
        hu = h.*u;
    end
    %----------------------------------------------------------------------
    u          = hu./h;
    d1h        = diff(h , x, 1);
    d1u        = diff(u , x, 1);
    d1hu       = diff(hu, x, 1);
    d1z        = diff(z , x, 1);
    d1b        = diff(zb, x, 1);
    d2h        = diff(h , x, 2);
    d2u        = diff(u , x, 2);
    d2hu       = diff(hu, x, 2);
    d2zb       = diff(zb, x, 2);
    d3zb       = diff(zb, x, 3);
    %----------------------------------------------------------------------
    ghd1z      = G.*h.*d1z;
    q1         = 2.*h.*(d1h+d1b./2).*d1u.^2+4./3.*h.^2.*d1u.*d2u+h.*d2zb.*u.*d1u+(d1h.*d2zb+h./2.*d3zb).*u.^2;
    hq1        = h.*q1;
    huu        = hu.*u;
    tp         = diff(huu, x, 1);
    phi        = diff(hu, t, 1)+tp+ghd1z;
    psi        = phi-ghd1z./alpha;
    psi_h      = psi./h;
    %----------------------------------------------------------------------
    div1       = diff(psi_h, x, 1);
    div2       = d1b.*psi_h;
    lhs        = psi-alpha.*(1./3.*diff(h.^3.*div1, x, 1)+1./2.*diff(h.^2.*div2, x, 1)-1./2.*h.^2.*div1.*d1b+h.*div2.*d1b);
    rhs        =-ghd1z./alpha-hq1;
    frict      =-G.*nm.^2.*h.^(-1/3).*(u.^2.*heaviside(u)-u.^2.*heaviside(-u)).*heaviside(h);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data.d1H   = mymatlabFunction(vpa(d1h  , 16), [x, t]);
    data.d1Hu  = mymatlabFunction(vpa(d1hu , 16), [x, t]);
    data.d1U   = mymatlabFunction(vpa(d1u  , 16), [x, t]);
    data.d1Z   = mymatlabFunction(vpa(d1z  , 16), [x, t]);
    data.d1Zb  = mymatlabFunction(vpa(d1b , 16), [x, t]);
    data.d2H   = mymatlabFunction(vpa(d2h  , 16), [x, t]);
    data.d2U   = mymatlabFunction(vpa(d2u  , 16), [x, t]);
    data.d2Hu  = mymatlabFunction(vpa(d2hu , 16), [x, t]);
    data.d2Zb  = mymatlabFunction(vpa(d2zb , 16), [x, t]);
    data.d3Zb  = mymatlabFunction(vpa(d3zb , 16), [x, t]);
    %----------------------------------------------------------------------
    data.HYD   = mymatlabFunction(vpa(ghd1z, 16), [x, t]);
    data.HQ1   = mymatlabFunction(vpa(hq1  , 16), [x, t]);
    data.PHI   = mymatlabFunction(vpa(phi  , 16), [x, t]);
    data.PSI   = mymatlabFunction(vpa(psi  , 16), [x, t]);
    data.PSI_H = mymatlabFunction(vpa(psi_h, 16), [x, t]);
    data.LHS   = mymatlabFunction(vpa(lhs  , 16), [x, t]);
    data.RHS   = mymatlabFunction(vpa(rhs  , 16), [x, t]);
    data.FRICT = mymatlabFunction(vpa(frict, 16), [x, t]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag && test ~= 12
    %----------------------------------------------------------------------
    figure('Color', 'w', 'Renderer', 'painters');
    subplot(1, 2, 1);
    hold on;
    p1 = plot(xv, data.Z (xv, 0), '-*r');
    p2 = plot(xv, data.Zb(xv, 0), '-og');
    p3 = plot(xv, data.H (xv, 0), '-sm'); legend([p1, p2, p3], '$\zeta$', '$z_{b}$', '$h$');
    subplot(1, 2, 2);
    p3 = plot(xv, data.HU(xv, 0), '-or'); legend(p3, '$hu$');
    %----------------------------------------------------------------------
    %{
    if ismembc(test, [0, 3, 6, 7, 8])
        time = 0;
        hq1t = data.HQ1(xv, time);
        lhst = data.LHS(xv, time);
        hydt = data.HYD(xv, time);
        phit = data.PHI(xv, time);
        psit = data.PSI(xv, time);
        figure('Color', 'w', 'Renderer', 'painters');
        hold on;
        p4 = plot(xv,-hq1t     , '-og');
        p5 = plot(xv, lhst+hydt, '-*r'); legend([p4, p5], 'Q1', 'LHS-RHS');
        figure('Color', 'w', 'Renderer', 'painters');
        hold on;
        p6 = plot(xv, phit     , '-*r');
        p7 = plot(xv, psit+hydt, '-ob'); legend([p6, p7], 'PHI', 'PSI+HYD');
    end
    %}
    %----------------------------------------------------------------------
    if test == 7
        nt   = 1000;
        t_   = linspace(0, tend, nt);
        xv_  = [xv; xM+(xv-xm)];
        var_ = data.PSI_H(xv_', t_');
        figure('Color', 'w', 'Renderer', 'painters');
        p8 = plot(xv_, var_(1, :), '-b');
        set(gca, ...
            'Box', 'on', ...
            'Clipping', 'on', ...
            'Layer', 'top', ...
            'TickLabelInterpreter', 'latex', ...
            'XMinorTick', 'on', ...
            'YMinorTick', 'on');
        for i = 2:nt
            set(p8, 'YData', var_(i, :));
            drawnow;
            disp(var_(i, K+1));
            pause(0.025);
        end
    end
    %----------------------------------------------------------------------
    close all;
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g] = mymatlabFunction(f, vars)
if ~has(f, vars(1, 1))
    c = matlabFunction(f, 'Vars', vars);
    g = @(x, varargin) c(varargin{:}).*ones(size(x));
else
    g = matlabFunction(f, 'Vars', vars);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%