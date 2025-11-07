function [data] = setdata(test, flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x t;
switch test
    case 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S0/T0 - Convergence (Green-Nagdhi)
        % NOTE: No need to make 3 steps for the RK3 scheme (dt = constant!)
        %------------------------------------------------------------------
        alpha    = 1;
        dispers  = 1;
        G        = 1;
        h0       = 1;
        nm       = 0;
        wetdry   = 0;
        %------------------------------------------------------------------
        xm       = 0;
        xM       = 200;
        K        = 2048;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        A1       = 0.2.*h0;
        c        = sqrt(G.*(h0+A1));
        k        = sqrt(3.*A1)./(2.*h0.*sqrt(h0+A1));
        xs       = 50;
        b        =-h0;
        z        = A1.*sech(k.*(x-xs-c.*t)).^2;
        h        = z-b;
        u        = c.*(1-h0./(h0+z));
        hu       = h.*u;
        tend     = 100./c;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S0/T1 - Convergence (SW)
        %------------------------------------------------------------------
        alpha    = 1;
        dispers  = 0;
        G        = 1;
        h0       = 1;
        nm       = 0;
        wetdry   = 0;
        %------------------------------------------------------------------
        xm       = 0;
        xM       = 200;
        K        = 512;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        h        = 5+exp(cos(2.*pi.*x));
        hu       = sin(cos(2.*pi.*x));
        b        = sin(pi.*x).^2;
        z        = h+b;
        u        = hu./h;
        tend     = 10;
        tk       = tend;
        %------------------------------------------------------------------
        % SOURCE:
        dt_z     = diff(z , t, 1);
        dthu     = diff(hu, t, 1);
        dxhu     = diff(hu, x, 1);
        dxzb     = diff(b , x, 1);
        hu2      = hu.^2./h;
        f1       = hu2+G./2.*(z.^2-2.*z.*b);
        dxf1     = diff(f1, x, 1);
        f2       = G.*z.*dxzb;
        s1       = dt_z+dxhu;
        s2       = dthu+dxf1+f2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S0/T2 - Wave generation
        %------------------------------------------------------------------
        alpha    = 1;
        %        = 1.159;
        dispers  = 0;
        G        = 9.81;
        h0       = 1;
        sqrth0_G = sqrt(h0./G);
        nm       = 0;
        wetdry   = 0;
        %------------------------------------------------------------------
        xm       = 0;
        xM       = 100;
        K        = 500;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        A1       = 0.20.*h0;
        c        = sqrt(G.*(h0+A1));
        k        = sqrt(3.*A1)./(2.*h0.*sqrt(h0+A1));
        xs       = 0;
        b        =-h0;
        z        = A1.*sech(k.*(x-xs-c.*t)).^2;
        h        = z-b;
        u        = c.*(1-h0./(h0+z));
        hu       = h.*u;
        %{
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
        %}
        tend     = 100.*sqrth0_G;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S1/T1 - Lake at rest
        %------------------------------------------------------------------
        alpha    = 1;
        dispers  = 0;
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
                b1 =-0.50.*exp(-50.*(x+1).^2); % WET
                b2 = 0.50.*exp(-25.*(x-1).^2); % DRY
            case 1 % LINEAR
                b1 =-0.50.*((x+1.50).*heaviside(x+1.50).*heaviside(-(1.00+x))+(-0.50-x).*heaviside(-0.50-x).*heaviside(x+1.00));
                b2 = 0.50.*((x-0.50).*heaviside(x-0.50).*heaviside( (1.00-x))+( 1.50-x).*heaviside( 1.50-x).*heaviside(x-1.00));
            case 2 % QUADRATIC
                b1 =-0.50.*(1-((x+1)./0.5).^2).*heaviside(x+1.5).*heaviside(-0.5-x);
                b2 = 0.50.*(1-((x-1)./0.5).^2).*heaviside(x-0.5).*heaviside( 1.5-x);
            otherwise
                return
        end
        b        = b1+b2;
        z        = h0+(b-h0).*(0.5.*(sign(b2.*heaviside(x)-h0)+1));
        h        = z-b;
        u        = 0;
        hu       = h.*u;
        tend     = 10e3;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 5
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S1/T2 - Drying of a lake
        %------------------------------------------------------------------
        alpha    = 1;
        dispers  = 0;
        G        = 1;
        h0       = 0.60;
        nm       = 0;
        wetdry   = 1;
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
        c1       = 0.25;
        c2       = 0.10;
        k        = pi./c2;
        b        = c1.*(cos(k.*x)+1).*heaviside(x+c2).*heaviside(c2-x);
        z        = h0+(b-h0).*heaviside(b-h0);
        h        = z-b;
        u        = 0;
        hu       = h.*u;
        tend     = 1000;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 6
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S1/T3 - Oscillating lake
        %------------------------------------------------------------------
        alpha    = 1;
        dispers  = 0;
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
        b        = h0.*x.^2;
        h_aux    = h0+2.*h0.*a0.*cos(omega.*t).*(x-a0./2.*cos(omega.*t))-b;
        h        = h_aux.*heaviside(h_aux);
        z        = h+b;
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
        dispers  = 1;
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
        b        =-h0;
        z1       = A1.*sech(k1.*(x-xs1-c1.*t)).^2;
        z2       = A2.*sech(k2.*(x-xs2+c2.*t)).^2;
        z        = z1+z2;
        h        = z-b;
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
        dispers  = 1;
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
        b        =-h0;
        z        = A1.*sech(k.*(x-xs-c.*t)).^2;
        h        = z-b;
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
        dispers  = 1;
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
        b        =-h0;
        z1       = A1.*sech(k1.*(x-xs1-c1.*t)).^2;
        z2       = A2.*sech(k2.*(x-xs2-c2.*t)).^2;
        z        = z1+z2;
        h        = z-b;
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
        dispers  = 1;
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
        b        =-h0;
        z        = exp(-x.^2./A1);
        h        = z-b;
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
        dispers  = 1;
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
        b        =-h0;
        z        = 1./2.*A1.*(1-tanh(x./0.4));
        h        = z-b;
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
        dispers  = 1;
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
                b =-h0+heaviside(x).*slope.*x;
            case 1
                f1 = bathyC1(slope, ds, 1);
                b  =-h0+...
                    (heaviside(x+ds)-heaviside(x-ds)).*f1(x)+...
                    (heaviside(x-ds)).*slope.*x;
            case 2
                f2 = bathyC2(slope, ds, 1);
                b  =-h0+...
                    (heaviside(x+ds)-heaviside(x-ds)).*f2(x)+...
                    (heaviside(x-ds)).*slope.*x;
            case 3
                f3 = bathyC3(slope, ds, 1);
                b  =-h0+...
                    (heaviside(x+ds)-heaviside(x-ds)).*f3(x)+...
                    (heaviside(x-ds)).*slope.*x;
            otherwise
                return
        end
        z        = A1.*sech(k.*(x-xs-c.*t)).^2;
        h        = z-b;
        u        = c.*(1-h0./(h0+z));
        hu       = h.*u;
        xg       = [-5.00, 20.96, 22.55, 23.68, 24.68, 25.91]'.*h0;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 13
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S3/T2 - Reversibility check (shelf)
        %------------------------------------------------------------------
        alpha    = 1;
        %        = 1.159;
        dispers  = 1;
        G        = 9.81;
        h0       = 1;
        nm       = 0;
        wetdry   = 0;
        option   = 3;
        %
        bathy    = 3;
        slope    = 1./20;
        %------------------------------------------------------------------
        xm       = 0;
        xM       = 240;
        K        = 2400;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        A1       = 0.20.*h0;
        c        = sqrt(G.*(h0+A1));
        k        = sqrt(3.*A1)./(2.*h0.*sqrt(h0+A1));
        xs       = 80;
        switch bathy
            case 0
                b  =-h0+...
                    slope.*(x-130).*(heaviside(x-130)-heaviside(x-140))+...
                    slope.*10.*heaviside(x-140);
            case 3
                ds = 0.5;
                f3 = bathyC3(slope, ds, 1);
                g3 = bathyC3(slope, ds, 3);
                b  =-h0+...
                    (heaviside(x-(130-ds))-heaviside(x-(130+ds))).*f3(x-130)+...
                    (heaviside(x-(130+ds))-heaviside(x-(140-ds))).*slope.*(x-130)+...
                    (heaviside(x-(140-ds))-heaviside(x-(140+ds))).*g3(x-140)+...
                    slope.*10.*heaviside(x-(140-ds));
            otherwise
                return
        end
        z        = A1.*sech(k.*(x-xs-c.*t)).^2;
        h        = z-b;
        u        = c.*(1-h0./(h0+z));
        hu       = h.*u;
        tend     = 50;
        tk       = tend;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 14
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % S3/T3 - Reflection of shoaling waves (M. Walkley)
        %------------------------------------------------------------------
        alpha    = 1;
        %        = 1.159;
        dispers  = 1;
        G        = 9.81;
        h0       = 0.7;
        nm       = 0;
        wetdry   = 0;
        option   = 1;
        %
        bathy    = 0; %#ok<NASGU> 
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
        b        =-h0+heaviside(x).*slope.*(x);
        h        = z-b;
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
        alpha    = 1;
        %        = 1.159;
        dispers  = 1;
        G        = 9.81;
        h0       = 1;
        sqrth0_G = sqrt(h0./G);
        nm       = 0;
        wetdry   = 1;
        option   = 2;
        %
        bathy    = 0;
        slope    = 1./19.85;     
        %------------------------------------------------------------------
        xm       =-20;
        xM       = 80;
        K        = 300;
        xv       = linspace(xm, xM, K+1)';
        dx       = zeros(K, 1);
        for i = 1:K
            dx(i, 1) = xv(i+1, 1)-xv(i, 1);
        end
        %------------------------------------------------------------------
        switch option
            case 1
                A1   = 0.019.*h0;
                tend = 100.*sqrth0_G;
                tk   = (25:5:70).*sqrth0_G;
            case 2
                A1   = 0.040.*h0;
                tend = 100.*sqrth0_G;
                tk   = (20:6:62).*sqrth0_G;
            otherwise
                return
        end
        c        =-sqrt(G.*(h0+A1));
        k        = sqrt(3.*A1)./(2.*h0.*sqrt(h0+A1));
        gamma    = sqrt(3.*A1./(4.*h0));
        xs       = 1./slope+1./gamma.*acosh(sqrt(1./0.05));
        ds       = 0.5;
        switch bathy
            case 0
                b  =-h0+...
                    heaviside(-x+1./slope).*slope.*(-x+1./slope);
            case 3
                f3 = bathyC3(slope, ds, 2);
                b  =-h0+...
                    slope.*(-x+1./slope).*heaviside(-x+1./slope-ds)+...
                    f3(x-1./slope).*(heaviside(x-(1./slope-ds))-heaviside(x-(1./slope+ds)));
            otherwise
                return
        end
        z        = A1.*sech(k.*(x-xs-c.*t)).^2;
        z        = z.*heaviside(x-eps)+heaviside(-(x+eps)).*b; %?
        h        = z-b;
        u        = c.*(1-h0./(h0+z));
        hu       = h.*u;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.Z       = mymatlabFunction(vpa(z , 16), [x, t]); % sol1
data.HU      = mymatlabFunction(vpa(hu, 16), [x, t]); % sol2
data.B       = mymatlabFunction(vpa(b , 16), [x, t]);
data.H       = mymatlabFunction(vpa(h , 16), [x, t]); % NOT NECESSARY
data.U       = mymatlabFunction(vpa(u , 16), [x, t]); % NOT NECESSARY
if ismembc(test, 2)
    data.S{1, 1} = mymatlabFunction(vpa(s1, 16), [x, t]); % TAKE CONVERGENCE
    data.S{1, 2} = mymatlabFunction(vpa(s2, 16), [x, t]);
end
d1b          = diff(b , x, 1);
d1h          = diff(h , x, 1);
d1u          = diff(u , x, 1);
d1z          = diff(z , x, 1);
d1hu         = diff(hu, x, 1);
d2b          = diff(b , x, 2);
d2h          = diff(h , x, 2);
d2u          = diff(u , x, 2);
d2hu         = diff(hu, x, 2);
d3b          = diff(b , x, 3);
data.d1B     = mymatlabFunction(vpa(d1b , 16), [x, t]);
data.d1H     = mymatlabFunction(vpa(d1h , 16), [x, t]);
data.d1U     = mymatlabFunction(vpa(d1u , 16), [x, t]);
data.d1Z     = mymatlabFunction(vpa(d1z , 16), [x, t]);
data.d1Hu    = mymatlabFunction(vpa(d1hu, 16), [x, t]);
data.d2B     = mymatlabFunction(vpa(d2b , 16), [x, t]);
data.d2H     = mymatlabFunction(vpa(d2h , 16), [x, t]);
data.d2U     = mymatlabFunction(vpa(d2u , 16), [x, t]);
data.d2Hu    = mymatlabFunction(vpa(d2hu, 16), [x, t]);
data.d3B     = mymatlabFunction(vpa(d3b , 16), [x, t]);
%--------------------------------------------------------------------------
data.alpha   = alpha;
data.dispers = dispers;
data.G       = G;
data.h0      = h0;
data.nm      = nm;
data.wetdry  = wetdry;
data.tk      = tk;
data.tend    = tend;
%--------------------------------------------------------------------------
data.xv      = xv;
data.nf      = 1800;
%--------------------------------------------------------------------------
if ismembc(test, [1, 7, 8])
    ghd1z      = G.*h.*d1z;
    q1         = 2.*h.*(d1h+d1b./2).*d1u.^2+4./3.*h.^2.*d1u.*d2u+h.*d2b.*u.*d1u+(d1h.*d2b+h./2.*d3b).*u.^2;
    hq1        = h.*q1;
    huu        = hu.*u;
    tp         = diff(huu, x, 1);
    phi        = diff(hu, t, 1)+tp+ghd1z;
    psi        = phi-ghd1z./alpha;
    psi_h      = psi./h;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag
    f1 = figure('Color', 'w', 'Renderer', 'painters');
    subplot(1, 2, 1);
    hold on;
    p1 = plot(xv, data.Z (xv, 0), '-*r');
    p2 = plot(xv, data.B (xv, 0), '-og');
    p3 = plot(xv, data.H (xv, 0), '-sm'); legend([p1, p2, p3], '$\zeta$', '$b$', '$h$');
    subplot(1, 2, 2);
    hold on;
    p3 = plot(xv, data.HU(xv, 0), '-or'); legend(p3, '$hu$');
    close(f1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ismembc(test, [0, 1, 6, 7, 8, 14, 15, 17])
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if test == 3 || test == 7
%         z1 = z;
%         z2 = A1.*sech(k.*(x-(2.*xM-xs)+c.*t)).^2;
%         z  = z1+z2;
%         if test == 3
%             zb = zb+data.Zb(2.*xM-x, t);
%         end
%         h  = z-zb;
%         u1 = u;
%         u2 =-c.*(1-h0./(h0+z2));
%         u  = u1+u2;
%         hu = h.*u;
%     end
%     %----------------------------------------------------------------------
% 
%     %----------------------------------------------------------------------
%     div1       = diff(psi_h, x, 1);
%     div2       = d1b.*psi_h;
%     lhs        = psi-alpha.*(1./3.*diff(h.^3.*div1, x, 1)+1./2.*diff(h.^2.*div2, x, 1)-1./2.*h.^2.*div1.*d1b+h.*div2.*d1b);
%     rhs        =-ghd1z./alpha-hq1;
%     frict      =-G.*nm.^2.*h.^(-1/3).*(u.^2.*heaviside(u)-u.^2.*heaviside(-u)).*heaviside(h);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     %----------------------------------------------------------------------
%     data.HYD   = mymatlabFunction(vpa(ghd1z, 16), [x, t]);
%     data.HQ1   = mymatlabFunction(vpa(hq1  , 16), [x, t]);
%     data.PHI   = mymatlabFunction(vpa(phi  , 16), [x, t]);
%     data.PSI   = mymatlabFunction(vpa(psi  , 16), [x, t]);
%     data.PSI_H = mymatlabFunction(vpa(psi_h, 16), [x, t]);
%     data.LHS   = mymatlabFunction(vpa(lhs  , 16), [x, t]);
%     data.RHS   = mymatlabFunction(vpa(rhs  , 16), [x, t]);
%     data.FRICT = mymatlabFunction(vpa(frict, 16), [x, t]);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if flag && test ~= 12
%     %----------------------------------------------------------------------
% 
%     %----------------------------------------------------------------------
%     %{
%     if ismembc(test, [0, 3, 6, 7, 8])
%         time = 0;
%         hq1t = data.HQ1(xv, time);
%         lhst = data.LHS(xv, time);
%         hydt = data.HYD(xv, time);
%         phit = data.PHI(xv, time);
%         psit = data.PSI(xv, time);
%         figure('Color', 'w', 'Renderer', 'painters');
%         hold on;
%         p4 = plot(xv,-hq1t     , '-og');
%         p5 = plot(xv, lhst+hydt, '-*r'); legend([p4, p5], 'Q1', 'LHS-RHS');
%         figure('Color', 'w', 'Renderer', 'painters');
%         hold on;
%         p6 = plot(xv, phit     , '-*r');
%         p7 = plot(xv, psit+hydt, '-ob'); legend([p6, p7], 'PHI', 'PSI+HYD');
%     end
%     %}
%     %----------------------------------------------------------------------
%     if test == 7
%         nt   = 1000;
%         t_   = linspace(0, tend, nt);
%         xv_  = [xv; xM+(xv-xm)];
%         var_ = data.PSI_H(xv_', t_');
%         figure('Color', 'w', 'Renderer', 'painters');
%         p8 = plot(xv_, var_(1, :), '-b');
%         set(gca, ...
%             'Box', 'on', ...
%             'Clipping', 'on', ...
%             'Layer', 'top', ...
%             'TickLabelInterpreter', 'latex', ...
%             'XMinorTick', 'on', ...
%             'YMinorTick', 'on');
%         for i = 2:nt
%             set(p8, 'YData', var_(i, :));
%             drawnow;
%             disp(var_(i, K+1));
%             pause(0.025);
%         end
%     end
%     %----------------------------------------------------------------------
%     close all;
%     %----------------------------------------------------------------------
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end