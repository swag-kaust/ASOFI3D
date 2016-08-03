function [d,hcerr] = gethbcoeff(L,dx,E,outp,cntlim)

global LL EE;
LL = L;
EE = E;

if nargin < 4
    outp = 1;
end
if nargin < 5
    cntlim = 10;
end

% Startwerte f�r Wellenzahl k f�r dx=2 -> dx/2=1 -> Normierung
k0(2,:)  = [ 0.0  0.0   0.0   0.0   0.0   0.0  ];    % L = 2
k0(4,:)  = [ 0.0  0.26  0.0   0.0   0.0   0.0  ];    % L = 4
k0(6,:)  = [ 0.0  0.32  0.56  0.0   0.0   0.0  ];    % L = 6
k0(8,:)  = [ 0.0  0.32  0.58  0.77  0.0   0.0  ];    % L = 8
k0(10,:) = [ 0.0  0.34  0.58  0.78  1.09  0.0  ];    % L = 10
k0(12,:) = [ 0.0  0.32  0.60  0.79  0.98  1.17 ];    % L = 12

% Startwerte f�r Koeff. d sind Taylor-Koeff. (Spalten: 1,3,5,...)
d0(2,:)  = [ 1               0             0              0             0             0          ];
d0(4,:)  = [ 9/8            -1/24          0              0             0             0          ];
d0(6,:)  = [ 75/64          -25/384        3/640          0             0             0          ];
d0(8,:)  = [ 1225/1024      -245/3072      49/5120       -5/7168        0             0          ];
d0(10,:) = [ 19845/16384    -735/8192      567/40960     -405/229376    35/294912     0          ];
d0(12,:) = [ 160083/131072  -12705/131072  22869/131072  -5445/1835008  847/2359296  -63/2883584 ];


options = optimset('levenbergmarquardt','on',...
                   'maxfunevals', 1000,...
                   'maxiter', 1000,...
                   'tolfun',1e-8,...
                   'tolx', 1e-8);
%                    'typicalx',[1;-0.1;0;1];...

% Startvektor
dk0 = [d0(LL,1:LL/2) k0(LL,1:LL/2)]';

% Gleichungssystem l�sen
warning off all;
cnt = 0;
while (cnt < cntlim)
    [x,y,exitflag,output] = fsolve(@tgs2,dk0,options);
    if sum(abs(y)) <= 1e-8
        break
    end
    dk0(1:LL/2) = dk0(1:LL/2)/1.05;
    dk0(LL/2+1:end) = dk0(LL/2+1:end)*1.02;
    cnt = cnt + 1;
end
warning on all;

lx = length(x)/2 + 1;
% Holberg-Koeffizienten
d = x(1:lx-1);
% Gl.-system berechnet sin(c*k), tats�chlich sin(c*k*dx/2), deswegen:
k = x(lx:end)/(dx/2);  % Wellenzahl k (Norm. im GS rückgängig machen)
k0 = k0/(dx/2);        % dito


% Vgl. zu Startwerten
qd = d./(d0(LL,1:LL/2)');
qk = k./(k0(LL,1:LL/2)');

% Ausgabe
if outp
    fprintf('==============================================================================\n');
    fprintf('Optimization results\n');
    fprintf('Holberg-coefficients for order L = %d and max. error E = %g %%:\n',LL,100*EE);
    for n = 1:lx-1
        fprintf('d(%d) = %.6G \t\t = %.2f*d0(%d)\n',2*n-1,d(n),qd(n),2*n-1);
    end

    fprintf('\nExtrema at wavenumber\n');
    for n = 1:lx-1
        fprintf('k(%d) = %.6G 1/m \t\t = %.2f*k0(%d)\n',n,k(n),qk(n),n);
    end

    fprintf('\nverify results (should be 0):\n');
    for n = 1:lx-1
        fprintf('y[d(%d)] = %.6g\n',2*n-1,y(n));
    end
    for n = 1:lx-1
        fprintf('y[k(%d)] = %.6g\n',n,y(n));
    end
    fprintf('==============================================================================\n');
end

hcerr = 0;
if sum(abs(y)) > 1e-8
    warning('At least 1 wrong result!');
    hcerr = 1;
end
