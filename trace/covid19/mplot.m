clear all;
close all;
clc;

data = load('trace.dat');

h1 = figure(1);
% set(h1,'Visible', 'off'); 
hold on
plot(data(:,1), data(:,3)*100,'b','LineWidth',2)
plot(data(:,1), data(:,4)*100,'m','LineWidth',2)
xlabel('time (days)')
ylabel('population (percent of total population)')
grid on;
legend('Exposed','Infectious','Location', 'South')
title('Exposed and Infectious population from 2020 May 11 in France')
print -djpeg EI

h2 = figure(2);
% set(h2,'Visible', 'off'); 
hold on
plot(data(:,1), data(:,2)*100,'k','LineWidth',2)
plot(data(:,1), data(:,5)*100,'r','LineWidth',2)
xlabel('time (days)')
ylabel('population (percent of total population)')
grid on;
legend('Susceptible','Removed','Location', 'South')
title('Susceptible and Removed population from 2020 May 11 in France')
print -djpeg SR

h3 = figure(3);
% set(h3,'Visible', 'off'); 
hold on
R0 = 3.4;
plot(data(:,1), (1-data(:,10))*R0,'b','LineWidth',2)
plot(data(:,1), ones(length(data(:,1)),1)*R0,'--b','LineWidth',2)
xlabel('time (days)')
ylabel('controlled contagiousness = R_0(1-u)')
title('controlled contagiousness = R_0(1-u)')
grid on;
print -djpeg Rt

