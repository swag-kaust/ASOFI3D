% Load ASCII-data
load /seismics/data1/SUPGE053/MILET/C/tr94_24c.asc
trace=tr94_24c;
clear tr94_24c;

% Old sample intervall
dt=0.001;
n=size(trace,1);
t=dt:dt:n*dt;
trace=trace/max(trace);
plot(t,trace,'+')
%pause

% new sampling rate (corresponds to FD time step):
dtnew=0.00025;
tnew=dtnew:dtnew:n*dt;
% spline interpolation
trace_new=spline(t,trace,tnew);
hold on
plot(tnew,trace_new,'b')

% write interpolated data into new file:
fid=fopen('/seismics/data1/SUPGE053/FDVISKO/MILET2D/signals/tr94_24c_25.asc','w');
n_new=size(trace_new,2);
for i=1:n_new, fprintf(fid,'%10.5e\n',trace_new(i)); end
fclose(fid);
