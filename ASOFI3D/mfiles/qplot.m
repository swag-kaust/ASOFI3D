%Plotting of Q(f,L,f_l,tau) for a generalized SLS

% for plotting
Q0=50.0;                    
fp1=1.0;fp2=200;df=0.1;               

%------------------------------------------------------------------
% Relaxation frequencies and tau od

L=1;
tau1=2.0/Q0;
fl1=[70.0];

%----------------------------------------------------------
f=fp1:df:fp2;
w=2*pi*f;


ts=1.0./(2.0*pi*fl1);
te=ts*(tau1+1.0);
Qft1=qgsls(te,ts,1,w);



;
p_h1=plot(f,Qft1,'k-');
xlabel('frequency [hz]');
ylabel('quality factor (Q)');
axis([fp1 fp2 10 200]);




