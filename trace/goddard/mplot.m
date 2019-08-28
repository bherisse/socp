clear all;
close all;
clc;

data = load('trace.dat');

h1 = figure(1);
set(h1,'Visible', 'off'); 
plot(data(:,2), data(:,3),'LineWidth',2)
xlabel('x')
ylabel('y')
grid on;
title('xy trajectory')
print -depsc xy_trajectory

h2 = figure(2);
% set(h2,'Visible', 'off'); 
plot(data(:,2), data(:,4),'LineWidth',2)
xlabel('x')
ylabel('z')
grid on;
title('xz trajectory')
print -depsc xz_trajectory

h3 = figure(3);
set(h3,'Visible', 'off'); 
hold on
plot(data(:,1), data(:,5),'LineWidth',2)
plot(data(:,1), data(:,6),'k','LineWidth',2)
plot(data(:,1), data(:,7),'r','LineWidth',2)
xlabel('t')
ylabel('speed')
grid on;
title('Velocity')
ylim([-1 1])
legend('vx', 'vy', 'vz')
print -depsc velocity

h4 = figure(4);
set(h4,'Visible', 'off'); 
plot(data(:,1), data(:,8),'LineWidth',2)
xlabel('t')
ylabel('mass')
ylim([0 1])
title('mass')
grid on;
print -depsc mass

h5 = figure(5);
set(h5,'Visible', 'off'); 
hold on
plot(data(:,1), data(:,16), 'LineWidth',2)
plot(data(:,1), data(:,17),'k','LineWidth',2)
plot(data(:,1), data(:,18),'r','LineWidth',2)
plot(data(:,1), 1*ones(length(data(:,1)),1),'--','LineWidth',2)
plot(data(:,1), -1*ones(length(data(:,1)),1),'--','LineWidth',2)
xlabel('t')
ylabel('u')
grid on;
ylim([-1.5 1.5])
title('Control')
legend('ux', 'uy', 'uz')
print -depsc control

h6 = figure(6);
set(h6,'Visible', 'off'); 
hold on
plot(data(:,1), sqrt(data(:,5).^2+data(:,6).^2+data(:,7).^2),'LineWidth',2);
xlabel('t')
ylabel('speed norm')
grid on;
title('Velocity norm')
print -depsc velocityNorm

h7 = figure(7);
set(h7,'Visible', 'off'); 
plot(data(:,1), data(:,19),'LineWidth',2)
xlabel('t')
ylabel('hamiltonian')
grid on;
title('Hamiltonian')
print -depsc hamiltonian

h8 = figure(8);
% set(h8,'Visible', 'off'); 
hold on
plot(data(:,1), sqrt(data(:,16).^2 + data(:,17).^2 + data(:,18).^2),'LineWidth',2)
plot(data(:,1), 1*ones(length(data(:,1)),1),'--','LineWidth',2)
xlabel('t')
ylabel('u')
grid on;
ylim([0 1.5])
title('Norm of the control')
print -depsc controlNorm
