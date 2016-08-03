% Plot der Benchmark Ergebnisse
% 
% Daniel Koehn
% Freiberg, den 16.4.2007

clear all
close all

data=load('test.log_2nd_chevron.timings');

figure;
plot(data(:,1),data(:,2),'b-',data(:,1),data(:,3),'c-',data(:,1),data(:,4),'y-',data(:,1),data(:,5),'g-',data(:,1),data(:,6),'r-');

% Achsen und Beschriftungen
af=xlabel('timestep');
bf=ylabel('calculation time [s]');
df=title('Benchmark Results (CHIC), acoustic with malloc, 1 shot 32 CPU, Chevron model - 25400 secs');
ef=legend('V-Update','S-Update','V-Exchange','S-Exchange','Total',2);

% Anpassen der Fonts
set(gca,'FontSize',15,'FontName','Times');
set(af,'FontSize',15,'FontName','Times');
set(bf,'FontSize',15,'FontName','Times');
set(df,'FontSize',15,'FontName','Times');
set(ef,'FontSize',15,'FontName','Times');

sum=0.0;
for i=1:6000
sum = sum + data(i,6);
end

sum
