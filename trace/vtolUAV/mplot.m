clear all;
close all;
clc;

%% run algorithm
modeMPP = 0; % 0: (MPP)_{1,sigma} ; 1: (MPP)_{1,0} ; 2: (MPP)_{1,inf}
sigma = 2;
mu = 1;
param = [modeMPP, sigma, mu];
args = [];
for k=1:3
    args = [args num2str(param(k),'%10.5e\n') ' '];
end

cd '../../binaries/Release'
command = ['testVtolUAV.exe' ' ' args];
system(command);
cd '../../trace/vtolUAV'

%% disp parameters for nominal
file = 'trace.dat';
color = [0,0,1];
lineStyle = '-';
islegend = 0;

%% load data
data = load(file);
dataWP = fopen('../../data/vtolUAV/waypoints');
dataObs = fopen('../../data/vtolUAV/obstacles');

%% waypoints parameters
fgetl(dataWP); % get "n_wp:"
nWP = str2num(fgetl(dataWP));
fgetl(dataWP); % get "n_wp:"
posWP = zeros(nWP,3);
for k=1:nWP
	posWP(k,:) = str2num(fgetl(dataWP));
end

%% obstacle parameters
fgetl(dataObs); % get "n:"
nObs = str2num(fgetl(dataObs));
typeObs = zeros(nObs,1);
posObs = zeros(nObs,3);
radiusObs = zeros(nObs,3);
fgetl(dataObs); % get "type:"
for k=1:nObs
    typeObs(k) = str2num(fgetl(dataObs));	
end
fgetl(dataObs); % get "position:"
for k=1:nObs
    posObs(k,:) = str2num(fgetl(dataObs));
end
fgetl(dataObs); % get "radius:"
for k=1:nObs
    radiusObs(k,:) = str2num(fgetl(dataObs));
end
% include the window 
posWindow = (posObs(nObs,:)+posObs(nObs-1,:))/2;
radiusWindow = [radiusObs(nObs,1), radiusObs(nObs,2), 6];
posObs(nObs+1,:) = posWindow;
radiusObs(nObs+1,:) = radiusWindow;

%% figures
h1 = figure(1);
% set(h1,'Visible', 'off'); 
hold on
for k=1:nObs-2
	%rectangle('Position', [posObs(k,1)-radiusObs(k,1), posObs(k,2)-radiusObs(k,2), 2*radiusObs(k,1), 2*radiusObs(k,2)], 'Curvature', [1-typeObs(k) 1-typeObs(k)], 'FaceColor', 'k', 'EdgeColor', 'none')
    xr = [posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)-radiusObs(k,1)];
    yr = [posObs(k,2)-radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2)];
    %zr = [posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3)];
    patch(xr,yr,'k', 'EdgeColor', 'none') 
end
k = nObs+1;
xr = [posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)-radiusObs(k,1)];
yr = [posObs(k,2)-radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2)];
patch(xr,yr,'c', 'EdgeColor', 'none') 
plot(data(:,2), data(:,3), 'LineWidth', 2)
% plot(posWP(:,1), posWP(:,2), 'LineWidth', 2)
plot(posWP(1,1), posWP(1,2), '.g', 'markersize', 20);
plot(posWP(nWP,1), posWP(nWP,2), '.r', 'markersize', 20);
for k=2:nWP-1
	plot(posWP(k,1), posWP(k,2), '.k', 'markersize', 20);
end
% alpha(0.2)
xlabel('x (m)','FontSize',20,'Interpreter','latex')
ylabel('y (m)','FontSize',20,'Interpreter','latex')
axis equal
ylim([-10 110])
xlim([4 100])
box on
% title('xy trajectory')
print -depsc xy_trajectory

h2 = figure(2);
set(h2,'Visible', 'off'); 
plot(data(:,2), data(:,4), 'LineWidth', 2)
xlabel('x')
ylabel('z')
grid on;
title('xz trajectory')
print -depsc xz_trajectory

h3 = figure(3);
set(h3,'Visible', 'off'); 
hold on
plot(data(:,1), data(:,5), 'LineWidth', 2)
plot(data(:,1), data(:,6),'k', 'LineWidth', 2)
plot(data(:,1), data(:,7),'r', 'LineWidth', 2)
xlabel('t')
ylabel('speed')
ylim([-2 2])
grid on
title('Velocity')
legend('vx', 'vy', 'vz')
print -depsc velocity

h4 = figure(4);
set(h4,'Visible', 'off'); 
hold on
plot(data(:,1), data(:,14), 'LineWidth', 2)
plot(data(:,1), data(:,15),'k', 'LineWidth', 2)
plot(data(:,1), data(:,16),'r', 'LineWidth', 2)
plot(data(:,1), 1*ones(length(data(:,1)),1),'--', 'LineWidth', 2)
plot(data(:,1), -1*ones(length(data(:,1)),1),'--', 'LineWidth', 2)
xlabel('t')
ylabel('u')
grid on
ylim([-1.5 1.5])
title('Control')
legend('ux', 'uy', 'uz')
print -depsc control

h5 = figure(5);
% set(h5,'Visible', 'off'); 
hold on
plot(data(:,1)/data(end,1), sqrt(data(:,5).^2+data(:,6).^2+data(:,7).^2), 'LineStyle', lineStyle, 'LineWidth', 2, 'color', color)
% plot(data(:,1), 1*ones(length(data(:,1)),1),'--', 'LineWidth', 2)
xlabel('time (normalized)','FontSize',20,'Interpreter','latex')
ylabel('$\left\|v\right\|$ (m/s)','FontSize',20,'Interpreter','latex')
grid on
ylim([0 2])
% title('Velocity norm')
if islegend
    legend(legendLabels,'Interpreter','latex','Location','South','FontSize',20)
end
print -depsc velocityNorm

h6 = figure(6);
% set(h6,'Visible', 'off'); 
hold on
plot(data(:,1)/data(end,1), sqrt(data(:,14).^2+data(:,15).^2+data(:,16).^2), 'LineStyle', lineStyle, 'LineWidth', 2, 'color', color)
xlabel('t')
xlabel('time (normalized)','FontSize',20,'Interpreter','latex')
ylabel('$\left\|u + mge_3\right\|$ (normalized)','FontSize',20,'Interpreter','latex')
grid on
ylim([0 1.1])
if islegend
    legend(legendLabels,'Interpreter','latex','Location','South','FontSize',20)
    plot(data(:,1)/data(end,1), 1*ones(length(data(:,1)),1),'--r','LineWidth',2)
end
% title('Norm of the control')
print -depsc controlNorm

h7 = figure(7);
set(h7,'Visible', 'off'); 
plot(data(:,1), data(:,17),'LineWidth',2)
xlabel('t')
ylabel('hamiltonian')
grid on;
title('Hamiltonian')
print -depsc hamiltonian

h10 = figure(10);
% set(h10,'Visible', 'off'); 
% set(h10,'Renderer','Painters')
hold on
for k=1:nObs-2
	%rectangle('Position', [posObs(k,1)-radiusObs(k,1), posObs(k,2)-radiusObs(k,2), 2*radiusObs(k,1), 2*radiusObs(k,2)], 'Curvature', [1-typeObs(k) 1-typeObs(k)], 'FaceColor', 'k', 'EdgeColor', 'none')
    xr = [posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)-radiusObs(k,1)];
    yr = [posObs(k,2)-radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2)];
    %zr = [posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3)];
    patch(xr,yr,'k', 'EdgeColor', 'none') 
end
k = nObs+1;
xr = [posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)-radiusObs(k,1)];
yr = [posObs(k,2)-radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2)];
patch(xr,yr,'c', 'EdgeColor', 'none') 
plot(posWP(:,1), posWP(:,2), 'LineWidth', 2)
plot(posWP(1,1), posWP(1,2), '.g', 'markersize', 20);
plot(posWP(nWP,1), posWP(nWP,2), '.r', 'markersize', 20);
for k=2:nWP-1
	plot(posWP(k,1), posWP(k,2), '.k', 'markersize', 20);
end
% alpha(0.2)
xlabel('x (m)','FontSize',20,'Interpreter','latex')
ylabel('y (m)','FontSize',20,'Interpreter','latex')
axis equal
ylim([-10 110])
xlim([4 100])
box on
% title('xy path')
% print -depsc2 -painters xy_path
print -depsc2 xy_path
% print -dpdf xy_path


h11 = figure(11);
% set(h11,'Visible', 'off'); 
hold on
plot3(data(:,2), data(:,3), data(:,4), 'LineWidth', 2)
for k=1:nObs+1
    xr = [posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)-radiusObs(k,1)];
    yr = [posObs(k,2)-radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2)];
    zr = [posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3)];
    if (k< nObs+1)
        patch(xr,yr,zr,'k', 'EdgeColor', 'none')    
    else
        patch(xr,yr,zr,'c', 'EdgeColor', 'none')    
    end
    xr = [posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1)];
    yr = [posObs(k,2)-radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2)];
    zr = [posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)-radiusObs(k,3)];
    if (k< nObs+1)
        patch(xr,yr,zr,'k', 'EdgeColor', 'none')    
    else
        patch(xr,yr,zr,'c', 'EdgeColor', 'none')    
    end   
    xr = [posObs(k,1)-radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1)];
    yr = [posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2)];
    zr = [posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)-radiusObs(k,3)];
    if (k< nObs+1)
        patch(xr,yr,zr,'k', 'EdgeColor', 'none')    
    else
        patch(xr,yr,zr,'c', 'EdgeColor', 'none')    
    end
    xr = [posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)-radiusObs(k,1)];
    yr = [posObs(k,2)-radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2)];
    zr = [posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3)];
    if (k< nObs+1)
        patch(xr,yr,zr,'k', 'EdgeColor', 'none')    
    else
        patch(xr,yr,zr,'c', 'EdgeColor', 'none')    
    end
    xr = [posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1)];
    yr = [posObs(k,2)-radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2)];
    zr = [posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)-radiusObs(k,3)];
    if (k< nObs+1)
        patch(xr,yr,zr,'k', 'EdgeColor', 'none')    
    else
        patch(xr,yr,zr,'c', 'EdgeColor', 'none')    
    end
    xr = [posObs(k,1)-radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1)];
    yr = [posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2)];
    zr = [posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)-radiusObs(k,3)];
    if (k< nObs+1)
        patch(xr,yr,zr,'k', 'EdgeColor', 'none')    
    else
        patch(xr,yr,zr,'c', 'EdgeColor', 'none')    
    end
end
plot3(posWP(1,1), posWP(1,2), posWP(1,3), '.g', 'markersize', 20);
plot3(posWP(nWP,1), posWP(nWP,2), posWP(nWP,3), '.r', 'markersize', 20);
for k=2:nWP-1
	plot3(posWP(k,1), posWP(k,2), posWP(k,3), '.k', 'markersize', 20);
end
alpha(0.2)                % set all patches transparency to 0.3
xlabel('x (m)','Interpreter','latex')
ylabel('y (m)','Interpreter','latex')
zlabel('z (m)','Interpreter','latex')
axis equal
ylim([-10 110])
xlim([4 100])
zlim([0 40])
box on
% title('xyz trajectory')
print -depsc xyz_trajectory
%alpha(1.0)                % set all patches transparency to 0.3

h12 = figure(12);
% set(h12,'Visible', 'off'); 
hold on
plot3(posWP(:,1), posWP(:,2), posWP(:,3), 'LineWidth', 2)
for k=1:nObs+1
    xr = [posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)-radiusObs(k,1)];
    yr = [posObs(k,2)-radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2)];
    zr = [posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3)];
    if (k< nObs+1)
        patch(xr,yr,zr,'k', 'EdgeColor', 'none')    
    else
        patch(xr,yr,zr,'c', 'EdgeColor', 'none')    
    end
    xr = [posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1)];
    yr = [posObs(k,2)-radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2)];
    zr = [posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)-radiusObs(k,3)];
    if (k< nObs+1)
        patch(xr,yr,zr,'k', 'EdgeColor', 'none')    
    else
        patch(xr,yr,zr,'c', 'EdgeColor', 'none')    
    end   
    xr = [posObs(k,1)-radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1)];
    yr = [posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2)];
    zr = [posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)-radiusObs(k,3)];
    if (k< nObs+1)
        patch(xr,yr,zr,'k', 'EdgeColor', 'none')    
    else
        patch(xr,yr,zr,'c', 'EdgeColor', 'none')    
    end
    xr = [posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)-radiusObs(k,1)];
    yr = [posObs(k,2)-radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2)];
    zr = [posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3)];
    if (k< nObs+1)
        patch(xr,yr,zr,'k', 'EdgeColor', 'none')    
    else
        patch(xr,yr,zr,'c', 'EdgeColor', 'none')    
    end
    xr = [posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1)];
    yr = [posObs(k,2)-radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)-radiusObs(k,2), posObs(k,2)-radiusObs(k,2)];
    zr = [posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)-radiusObs(k,3)];
    if (k< nObs+1)
        patch(xr,yr,zr,'k', 'EdgeColor', 'none')    
    else
        patch(xr,yr,zr,'c', 'EdgeColor', 'none')    
    end
    xr = [posObs(k,1)-radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)+radiusObs(k,1), posObs(k,1)-radiusObs(k,1), posObs(k,1)-radiusObs(k,1)];
    yr = [posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2), posObs(k,2)+radiusObs(k,2)];
    zr = [posObs(k,3)-radiusObs(k,3), posObs(k,3)-radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)+radiusObs(k,3), posObs(k,3)-radiusObs(k,3)];
    if (k< nObs+1)
        patch(xr,yr,zr,'k', 'EdgeColor', 'none')    
    else
        patch(xr,yr,zr,'c', 'EdgeColor', 'none')    
    end
end
plot3(posWP(1,1), posWP(1,2), posWP(1,3), '.g', 'markersize', 20);
plot3(posWP(nWP,1), posWP(nWP,2), posWP(nWP,3), '.r', 'markersize', 20);
for k=2:nWP-1
	plot3(posWP(k,1), posWP(k,2), posWP(k,3), '.k', 'markersize', 20);
end
alpha(0.2)                % set all patches transparency to 0.3
xlabel('x (m)','Interpreter','latex')
ylabel('y (m)','Interpreter','latex')
zlabel('z (m)','Interpreter','latex')
axis equal
ylim([-10 110])
xlim([4 100])
zlim([0 40])
box on
% title('xyz trajectory')
print -depsc xyz_path
%alpha(1.0)                % set all patches transparency to 0.3

%% close data
fclose(dataWP);
fclose(dataObs);

