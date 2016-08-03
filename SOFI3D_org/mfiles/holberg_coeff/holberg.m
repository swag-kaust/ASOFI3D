% *************************************************************************
% Calculation of Holberg coefficients
% for standard staggered grid FD operators
%
% Holberg, O., 1987, Computational aspects of the choice of operator and
% sampling interval for numerical differentiation in large-scale simulation
% of wave phenomena: Geophys. Prosp., 35, 629-655
%
% Kindelan, M., 1990,  On the construction and efficiency of staggered
% numerical differentiators for the wave equation: Geophysics, 55, 107-110
% *************************************************************************


clear all, close all, clc, pause(0.1);

dx     =  1.00;    % spacing in m
Ls     =  12;      % order L (for comparison of operators using Taylor and holberg coefficients)
E      =  0.01;   % maximum relative error (for comparison...)
method =  2;       % 1: usage of b-coeff. in table 1 in Kindelan et al. (order L=2,..,8)
                   % 2: solving Holberg's equations (attention: solutions may be incorrect)
dim    =  2;       % dimensionality of the problem (for calculation of numerical cost)
coeff  =  7;       % plot coefficient di (i=1,3,...) vs. E

Ev = logspace(-3,0,100)';     % several errors for plots y vs. E
k  = logspace(-2,1,10000)';   % vector of wave numbers k in 1/m
Lv = 2:2:12;                  % comparison of different orders

titel    = 0;        % plot title? 
elimplot = [-40 10]; % error axis limitation (where rel. error is ordinate)

% old FD coefficients from Taylor series
f(2,:)  = [ 1              0   0             0  0              0   0             0  0             0   0          ];
f(4,:)  = [ 9/8            0  -1/24          0  0              0   0             0  0             0   0          ];
f(6,:)  = [ 75/64          0  -25/384        0  3/640          0   0             0  0             0   0          ];
f(8,:)  = [ 1225/1024      0  -245/3072      0  49/5120        0  -5/7168        0  0             0   0          ];
f(10,:) = [ 19845/16384    0  -735/8192      0  567/40960      0  -405/229376    0  35/294912     0   0          ];
f(12,:) = [ 160083/131072  0  -12705/131072  0  22869/1310720  0  -5445/1835008  0  847/2359296   0  -63/2883584 ];


% table of b coefficients to calculate Holberg coefficients
b(:,:,2) = [1.0 1.0 0.0 0.0 0.0 0.0;
            0.0 0.0 0.0 0.0 0.0 0.0;
            0.0 0.0 0.0 0.0 0.0 0.0;
            0.0 0.0 0.0 0.0 0.0 0.0];

b(:,:,4) = [ 1.125    0.4330 -0.4583  0.2566 0.0 0.0;
            -0.04167 -0.1443 -0.1806 -0.0855 0.0 0.0;
             0.0      0.0     0.0     0.0    0.0 0.0;
             0.0      0.0     0.0     0.0    0.0 0.0];

b(:,:,6) = [ 1.172     0.2742 -0.3006 0.2637 -0.2391 0.0;
            -0.06510  -0.1371 -0.0100 0.1077 -0.0086 0.0;
             0.004688  0.0274  0.0661 0.0826  0.0530 0.0;
             0.0       0.0     0.0    0.0     0.0    0.0];
  
b(:,:,8) = [ 1.1963     0.1987  -0.2195   0.2096  -0.2038   0.1779;
            -0.07975   -0.1192   0.03929  0.04049 -0.05841  0.04389;
             0.009570   0.03975  0.04853 -0.00967 -0.05507 -0.01117;
            -0.0006975 -0.00568 -0.02014 -0.04039 -0.04937 -0.03625];
  

% ==================================
iv = find(Lv == Ls);
idxL = iv(1);


% ==================================
if (method == 1) && (Lv(end)>8)
    Lv = 2:2:8;
end
if (method == 1) && (Ls>8)
    error('L > 8 not supported by current method!');
end


% ==================================
% Holberg coeff. for one max. rel. error E
for L = Lv
    switch method
        case 1
            rv = 0:5;
            for l = 1:L/2
                btmp = b(l,:,L);
                d(2*l-1,L) = sum(btmp.*E.^(2*rv/L));
            end
        case 2
            if L == Ls
                [d(1:2:L-1,L),hcerr(L)] = gethbcoeff(L,dx,E,1);
            else
                [d(1:2:L-1,L),hcerr(L)] = gethbcoeff(L,dx,E,1);
            end
    end

    Dk(:,L) = 0*k;
    for l = 1:L/2
        Dk(:,L) = Dk(:,L) + d(2*l-1,L)*sin((2*l-1)*k*dx/2)/(dx/2);
    end

    dDk(:,L) = 0*k;
    for l = 1:L/2
        dDk(:,L) = dDk(:,L) + d(2*l-1,L)*cos((2*l-1)*k*dx/2)*(2*l-1);
    end

    eDk(:,L) = dDk(:,L) - 1;

    il = find(abs(eDk(:,L)) > 1.2*E);
    il = il(1);
    while abs(eDk(il,L)) > E
        il = il(1) - 1;
        if il < 1
            il = 1;
            break
        end
    end
    kcH(L) = k(il+1);
end


% ==================================
% rel. error for old FD coeff.
for L = Lv
    Fk(:,L) = 0*k;
    for l = 1:L/2
        Fk(:,L) = Fk(:,L) + f(L,2*l-1)*sin((2*l-1)*k*dx/2)/(dx/2);
    end

    dFk(:,L) = 0*k;
    for l = 1:L/2
        dFk(:,L) = dFk(:,L) + f(L,2*l-1)*cos((2*l-1)*k*dx/2)*(2*l-1);
    end

    eFk(:,L) = dFk(:,L) - 1;

    il = find(abs(eFk(:,L)) > 1.2*E);
    il = il(1);
    while abs(eFk(il,L)) > E
        il = il(1) - 1;
        if il < 1
            il = 1;
            break
        end
    end
    kcT(L) = k(il+1);
end


% ==================================
col1 = [0 0 0.95];
col2 = [0.8 0 0];


% ==================================
% Plot: rel error vs. k for one E and L
figure(1);
hold on;
box on;
Eiv = E*ones(length(k),1);
plot(k,-100*Eiv);
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineWidth', 1,...
            'color',[0.3 0.3 0.3]);
plot(k,100*Eiv);
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineWidth', 1,...
            'color',[0.3 0.3 0.3]);
ph1 = plot(k,100*eFk(:,Lv(idxL)));
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineWidth', 1,...
            'color',col1);
ph2 = plot(k,100*eDk(:,Lv(idxL)));
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineWidth', 1,...
            'color',col2);
legend([ph1 ph2],'Taylor','Holberg',3);
xlabel('Wave number k in m^{-1}','fontsize',11);
ylabel('Relative error \epsilon in %','fontsize',11);
t1 = sprintf('Relative error for an operator of order L = %d,',Ls);
t2 = sprintf('a maximum relative error E = %.2G and ', E);
t3 = sprintf(' = %.2G m',dx);
title([sprintf([t1 '\n' t2]) '\Deltax' sprintf(t3)], 'fontsize', 12);
set(gca,'xlim',[k(1) 1.1*kcH(Ls)],...
        'ylim',elimplot);
ax = axis;
line([kcH(Ls) kcH(Ls)],elimplot);
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineStyle', '--',...
            'LineWidth', 1,...
            'color',col2);
line([kcT(Ls) kcT(Ls)],elimplot);
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineStyle', '--',...
            'LineWidth', 1,...
            'color',col1);
th = text(kcH(Ls),(ax(3)+(ax(4)-ax(3))/10),sprintf('k_{c,H} = %.2G m^{-1}  ',kcH(Ls)));
set(th, 'color',col2,...
        'fontsize',12,...
        'horizontalalignment','right');
th = text(kcT(Ls),(ax(3)+(ax(4)-ax(3))/10),sprintf('k_{c,T} = %.2G m^{-1}  ',kcT(Ls)));
set(th, 'color',col1,...  
        'fontsize',12,...
        'horizontalalignment','right');
tx0 = ax(1) + (ax(2)-ax(1))/30;
th = text(tx0,180*E,'E');
set(th, 'color',[0 0 0 ],...
        'fontsize',12,...
        'horizontalalignment','left');
th = text(0.8*tx0,-200*E,'-E');
set(th, 'color',[0 0 0 ],...
        'fontsize',12,...
        'horizontalalignment','left');


% ==================================
% Plot: rel. error vs. number of grid points N for one E and L
figure(7);
hold on;
box on;
Eiv = E*ones(length(k),1);
Nv = 2*pi./(k*dx);
NcH = 2*pi./(kcH*dx);
NcH1 = NcH(Ls);
NcT = 2*pi./(kcT*dx);
NcT1 = NcT(Ls);
plot(Nv,-100*Eiv);
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineWidth', 1,...
            'color',[0.3 0.3 0.3]);
plot(Nv,100*Eiv);
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineWidth', 1,...
            'color',[0.3 0.3 0.3]);
ph1 = plot(Nv,100*eFk(:,Lv(idxL)));
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineWidth', 1,...
            'color',col1);
ph2 = plot(Nv,100*eDk(:,Lv(idxL)));
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineWidth', 1,...
            'color',col2);
legend([ph1 ph2],'Taylor','Holberg',4);
xlabel('Grid points per shortest wavelength N','fontsize',11);
ylabel('Relative error \epsilon in %','fontsize',11);
t1 = sprintf('Relative error for an operator of order L = %d,',Ls);
t2 = sprintf('a maximum relative error E = %g and ', E);
t3 = sprintf(' = %g m',dx);
title([sprintf([t1 '\n' t2]) '\Deltax' sprintf(t3)], 'fontsize', 12);
set(gca,'xlim',[2 5],...
        'ylim',elimplot,...
        'xscale','log');
ax = axis;
xtck = ceil(ax(1)):floor(ax(2));
set(gca,'xtickmode','manual',...
        'xtick',xtck);
line([NcH1 NcH1],elimplot);
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineStyle', '--',...
            'LineWidth', 1,...
            'color',col2);
line([NcT1 NcT1],elimplot);
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineStyle', '--',...
            'LineWidth', 1,...
            'color',col1);
th = text(NcH1,(ax(3)+(ax(4)-ax(3))/10),sprintf('  N_{c,H} = %.1f',NcH1));
set(th, 'color',col2,...
        'fontsize',12,...
        'horizontalalignment','left');
th = text(NcT1,(ax(3)+(ax(4)-ax(3))/10),sprintf('  N_{c,T} = %.1f',NcT1));
set(th, 'color',col1,...
        'fontsize',12,...
        'horizontalalignment','left');
th = text(1.05*ax(1),180*E,' E');
set(th, 'color',[0 0 0 ],...
        'fontsize',12,...
        'horizontalalignment','left');
th = text(1.05*ax(1),-200*E,'-E');
set(th, 'color',[0 0 0 ],...
        'fontsize',12,...
        'horizontalalignment','left');


% ==================================
% Plot: rel. error vs. k for one E and all L
figure(2);
hold on;
box on;
Eiv = E*ones(length(k),1);
plot(k,-100*Eiv);
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineWidth', 1,...
            'color',[0.3 0.3 0.3]);
plot(k,100*Eiv);
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineWidth', 1,...
            'color',[0.3 0.3 0.3]);

p = 0;
for L = [2 4 8 12]
    p = p + 1;
    plot(k,100*eFk(:,L),'--');
    set(findobj(gca,'Type','line','Color',[0 0 1]),...
                'LineWidth', 1,...
                'color',[1-(L-2)/(Lv(end)-2) 0 (L-2)/(Lv(end)-2)]);
    ph(p) = plot(k,100*eDk(:,L));
    set(findobj(gca,'Type','line','Color',[0 0 1]),...
                'LineWidth', 1,...
                'color',[1-(L-2)/(Lv(end)-2) 0 (L-2)/(Lv(end)-2)]);
end
ph2(1) = plot(nan,nan,'color',[0 0 0]);
ph2(2) = plot(nan,nan,'color',[0 0 0],'linestyle','--');
legend([ph ph2], 'L = 2','L = 4','L = 8','L = 12','Holberg','Taylor');
xlabel('Wave number k in m^{-1}','fontsize',11);
ylabel('Relative error \epsilon in %','fontsize',11);
set(gca,'xlim',[k(1) 1.05*max(kcH)],...
        'ylim',[-400*E,300*E]);
t1 = sprintf('Relative error for operators of order L = 2 ... 12');
t2 = sprintf('a maximum relative error E = %g and ', E);
t3 = sprintf(' = %g m',dx);
title([sprintf([t1 '\n' t2]) '\Deltax' sprintf(t3)], 'fontsize', 12);
drawnow;
pause(0.1);


% ==================================
% Plot: rel. error vs. N for one E and all L
figure(8);
hold on;
box on;
Eiv = E*ones(length(k),1);
plot(Nv,-100*Eiv);
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineWidth', 1,...
            'color',[0.3 0.3 0.3]);
plot(Nv,100*Eiv);
set(findobj(gca,'Type','line','Color',[0 0 1]),...
            'LineWidth', 1,...
            'color',[0.3 0.3 0.3]);

p = 0;
ph = [];
for L = [2 4 8 12]
    p = p + 1;
    plot(Nv,100*eFk(:,L),'--');
    set(findobj(gca,'Type','line','Color',[0 0 1]),...
                'LineWidth', 1,...
                'color',[1-(L-2)/(Lv(end)-2) 0 (L-2)/(Lv(end)-2)]);
    ph(p) = plot(Nv,100*eDk(:,L));
    set(findobj(gca,'Type','line','Color',[0 0 1]),...
                'LineWidth', 1,...
                'color',[1-(L-2)/(Lv(end)-2) 0 (L-2)/(Lv(end)-2)]);
end
ph2(1) = plot(nan,nan,'color',[0 0 0]);
ph2(2) = plot(nan,nan,'color',[0 0 0],'linestyle','--');
legend([ph ph2], 'L = 2','L = 4','L = 8','L = 12','Holberg','Taylor');
xlabel('Grid points per shortest wavelength N','fontsize',11);
ylabel('Relative error \epsilon in %','fontsize',11);
xl = [2 7];
set(gca,'xscale','log',...
        'xlim',xl,...
        'ylim',[-40 10],...
        'xtickmode','manual',...
        'xtick',xl(1):xl(2),...
        'xticklabel',xl(1):xl(2));
t1 = sprintf('Relative error for operators of order L = 2 ... 12');
t2 = sprintf('a maximum relative error E = %g and ', E);
t3 = sprintf(' = %g m',dx);
title([sprintf([t1 '\n' t2]) '\Deltax' sprintf(t3)], 'fontsize', 12);
drawnow;
pause(0.1);


% ==================================
% Test
% figure(100);
% hold on;
% box on;
% Eiv = E*ones(length(k),1);
% plot(Nv,-100*Eiv);
% set(findobj(gca,'Type','line','Color',[0 0 1]),...
%             'LineWidth', 1,...
%             'color',[0.3 0.3 0.3]);
% plot(Nv,100*Eiv);
% set(findobj(gca,'Type','line','Color',[0 0 1]),...
%             'LineWidth', 1,...
%             'color',[0.3 0.3 0.3]);
% 
% p = 0;
% ph = [];
% for L = [2 4 6 8 10 12]
%     p = p + 1;
%     plot(Nv,100*eFk(:,L),'--');
%     set(findobj(gca,'Type','line','Color',[0 0 1]),...
%                 'LineWidth', 1,...
%                 'color',[1-(L-2)/10 0 (L-2)/10]);
%     ph(p) = plot(Nv,100*eDk(:,L));
%     set(findobj(gca,'Type','line','Color',[0 0 1]),...
%                 'LineWidth', 1,...
%                 'color',[1-(L-2)/10 0 (L-2)/10]);
% end
% legend(ph, 'L = 2','L = 4','L = 6','L = 8','L = 10','L = 12');
% xlabel('Grid points per shortest wavelength N','fontsize',11);
% ylabel('Relative error \epsilon in %','fontsize',11);
% xl = [2 10];
% set(gca,'xscale','log',...
%         'xlim',xl,...
%         'ylim',[-40 10],...
%         'xtickmode','manual',...
%         'xtick',xl(1):xl(2),...
%         'xticklabel',xl(1):xl(2));
% t1 = sprintf('Relative error for operators of order L = 2 ... 12');
% t2 = sprintf('a maximum relative error E = %g and ', E);
% t3 = sprintf(' = %g m',dx);
% title([sprintf([t1 '\n' t2]) '\Deltax' sprintf(t3)], 'fontsize', 12);
% drawnow;
% pause(0.1);


% ==================================
% calculate Holberg coeff for several E
Lvold = Lv;
Lv = 2:2:8;
method = 1;

eDk = [];
kc = [];
for n = 1:length(Ev)
    for L = Lv
        switch method
            case 1
                rv = 0:5;
                for l = 1:L/2
                    btmp = b(l,:,L);
                    dtmp(l) = sum(btmp.*Ev(n).^(2*rv/L));
                end
            case 2
                warning off all;
                [dtmp(1:L/2), hcerr] = gethbcoeff(L,dx,Ev(n),0);
                if hcerr
                    dtmp(1:L/2) = nan;
                end
        end

        Dtmp = 0*k;
        for l = 1:L/2
            Dtmp = Dtmp + dtmp(l)*sin((2*l-1)*k*dx/2)/(dx/2);
        end

        dDktmp = 0*k;
        for l = 1:L/2
            dDktmp = dDktmp + dtmp(l)*cos((2*l-1)*k*dx/2)/(dx/2)*(2*l-1)*dx/2;
        end

        eDk(:,n,L) = dDktmp - 1;

        il = find(abs(eDk(:,n,L)) > 1.2*Ev(n));
        if ~isempty(il)
            il = il(1);
            while abs(eDk(:,n,L)) > Ev(n)
                il = il(1) - 1;
                if il < 1
                    il = 1;
                    break
                end
            end
            kcE(n,L) = k(il);
        else
            kcE(n,L) = nan;
        end
        
        if (coeff+1)/2 <= length(dtmp)
            dsav(n,L/2) = dtmp((coeff+1)/2);
        else
            dsav(n,L/2) = NaN;
        end
    end
end


% ==================================
% Plot: grid points per shortest wavelength vs. E
for L = Lv
    Nc(:,L) = 2*pi./kcE(:,L)/dx;
end
figure(3);
hold on;
box on;
for L = Lv
    plot(Ev,Nc(:,L));
    set(findobj(gca,'Type','line','Color',[0 0 1]),...
                'LineWidth', 1,...
                'color',[1-(L-2)/(Lv(end)-2) 0 (L-2)/(Lv(end)-2)]);
end
set(gca,'xscale','log',...
        'yscale','log',...
        'yticklabelmode', 'manual',...
        'ytick', [1 2 5 10 20 50 100],...
        'yticklabel', [1 2 5 10 20 50 100]);
legend('L = 2','L = 4','L = 6','L = 8');
xlabel('Maximum relative error E','fontsize',11);
ylabel('Grid points per shortest wavelength N_c','fontsize',11);
t1 = 'Number of gridpoints per shortest wavelength';
t2 = 'needed by Holberg''s differentiator';
title(sprintf([t1 '\n' t2]), 'fontsize', 12);


% ==================================
% Plot: numerical cost
for L = Lv
    nop = 3*L/2 - 1;
    C(:,L) = nop*Nc(:,L).^dim;
end
figure(4);
hold on;
box on;
for L = Lv
    plot(Ev,C(:,L));
    set(findobj(gca,'Type','line','Color',[0 0 1]),...
                'LineWidth', 1,...
                'color',[1-(L-2)/(Lv(end)-2) 0 (L-2)/(Lv(end)-2)]);
    set(gca,'xscale','log','yscale','log');
end
legend('L = 2','L = 4','L = 6','L = 8');
xlabel('Maximum relative error E','fontsize',11);
ylabel('Relative cost C','fontsize',11);
title('Cost of numerical differentiation', 'fontsize', 12);


d = d';


% ==================================
figure(5);
hold on;
box on;
n = 0;
for L = Lv
    n = n + 1;
    ph1(n) = plot(Lv-1,d(L,Lv-1));
    set(findobj(gca,'Type','line','Color',[0 0 1]),...
                'LineWidth', 1,...
                'color',[1-(L-2)/(Lv(end)-2) 0 (L-2)/(Lv(end)-2)]);
    plot(Lv-1,f(L,Lv-1),'--');
    set(findobj(gca,'Type','line','Color',[0 0 1]),...
                'LineWidth', 1,...
                'color',[1-(L-2)/(Lv(end)-2) 0 (L-2)/(Lv(end)-2)]);
end
xlabel('l','fontsize',11);
ylabel('Coefficient','fontsize',11);
ph2(1) = plot(nan,nan,'color',[0 0 0]);
ph2(2) = plot(nan,nan,'color',[0 0 0],'linestyle','--');
legend([ph1 ph2(1) ph2(2)],'L = 2','L = 4','L = 6','L = 8',...
        'Holberg', 'Taylor');


% ==================================
% Plot: selected coeff. vs. E
figure(6);
hold on;
box on;
legstr = {};
for L = Lv/2
    if isnan(sum(dsav(:,L)))
        continue;
    end
    ph1 = plot(Ev,dsav(:,L));
    cst = (L-1)*0.99/((Lv(end)/2-1));
    set(findobj(gca,'Type','line','Color',[0 0 1]),...
                'LineWidth', 1,...
                'color',[1-cst 0 cst]);
    legstr = [legstr {sprintf('L = %d',L*2)}];
end
legend(legstr,2);
xlabel('Maximum relative error E','fontsize',11);
ylabel(sprintf('Coefficient d_{%d}',coeff),'fontsize',11);
set(gca,'xscale','log',...
        'ylim',[min(min(dsav)) max(max(dsav))]);


% ==================================
% format figures
scs = get(0,'screensize');
for n = 1:8
    figure(n);
    set(gcf,'position',[10 10 0.6*scs(3) 0.39*scs(4)]);
    cax = gca;
    axprop = get(cax);
    
    if ~titel
        set(axprop.Title,'Visible','off');
        set(cax,'position', [0.1 0.11 0.8 0.85]);
    else
        set(cax,'position', [0.1 0.12 0.8 0.63]);
    end
    
    set(axprop.XLabel,'fontweight','bold',...
                      'fontsize',12);
    set(axprop.YLabel,'fontweight','bold',...
                      'fontsize',12);
    
    set(cax,'xgrid','off',...
            'ygrid','off');
    
    if n==1 || n==2 || n==7 || n==8
        set(cax,'ylim',elimplot,...
                'yticklabelmode','manual',...
                'ytick',[-40 -30 -20 -10 0 10],...
                'yticklabel',[-40 -30 -20 -10 0 10]);
    end
    
    drawnow;
end
pause(0.1);
figure(1);


% ==================================
% Output to command window
fprintf('\n\nHolberg coefficients for E = %.2G (row = order 2,4,...):\n',E);
format short g;
d2 = d(2:2:end,1:2:end);
disp(d2);
if exist('hcerr')
    badhc = find(hcerr);
    if ~isempty(badhc)
        fprintf('\nCoefficients for order %d might be incorrect!\n',badhc);
    end
end

if exist('sym')
    syms dsym dpr;
    dd = round(10000*d2)/10000;
    dsym = sym(dd);
%     dpr = pretty(dsym);
    fprintf('\n\n');
    disp(dsym);
end

fprintf('\n\nMinimum number of grid points per shortest wavelength:\n\tHolberg\t\t\tTaylor\n');
disp([NcH(Lvold)' NcT(Lvold)']);
fprintf('==============================================================================\n\n');
format long;
