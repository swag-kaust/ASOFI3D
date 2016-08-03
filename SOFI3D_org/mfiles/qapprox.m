close all,clear all,clc;
% Q-approximation using an improved tau-method (Blanch et al., 1995)
% otimization routine: Marquardt-Levenberg
global L w Qf1 Qf2

%-------------------------INPUT-PARAMETERS----------------------------

Q0  = 125.0;                    % constant Q to be approximated
fp1 = 0.5; fp2=20; df=0.01;    % within frequency range fp1,..., fp2
          
L = 3;                     % number of relaxation mechanisms      

fl_st=[0.2 2 20];      % L starting values for the relaxation frequencies

t=2/Q0;                 % starting value for tau

mode = 2;       % old or new least squares function

%-------------------------END: INPUT-PARAMETER-----------------------

%%

f=fp1:df:fp2;
w=2*pi*f;

Qf1=Q0+f*0;           % here constant Q-approximation, arbitrary
                      % frequency dependency of Q possible, define
                      % Q as function of frequency here

% Qf1=Q0+sin(2*pi*f/fp2)*Q0;


x=[fl_st t];

switch mode
    % old
    case 1
        % if optimization toolbox is installed use 'leastq'
        % defining option for otimization (see 'help foptions')
        options(1)=1;
        options(2)=0.1;
        options(3)=0.1;
        options(5)=1;     % =0 : Levenberg-Marquardt Method, =1: Gauss-Newton Method
        options(14)=1000;
        [x,options]=leastsq('qflt',x,options);

    % new
    case 2
        newoptions = optimset('TolFun',1e-6, 'TolX', 1e-6, 'MaxIter',1000,...
                              'LevenbergMarquardt','on', 'Display','iter-detailed',...
                              'TolCon',1e-6);
        [x] = lsqnonlin(@qflt,x,[],[],newoptions);
end


%% output of results:
RELAXATIONSFEQUENZ=x(1:L),
TAU=x(L+1),
sigma=sqrt(sum((Qf2-Qf1).*(Qf2-Qf1))/size(Qf2,2))*100/Q0;
GEW_REL_STANDARDABW=sigma,



% plot Q as function of frequency:
plot(f,Qf1,f,Qf2,'r--');
xlabel('frequency [hz]');
ylabel('quality factor (Q)');
%axis([fp1 fp2 10 100]);

