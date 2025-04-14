data = load('trace.dat');

h1 = figure(1);
set(h1,'Visible', 'off'); 
plot(data(:,2), data(:,3))
xlabel('x')
ylabel('y')
title('xy trajectory')
print -depsc xy_trajectory

h2 = figure(2);
set(h2,'Visible', 'off'); 
plot(data(:,2), data(:,4))
xlabel('x')
ylabel('z')
title('xz trajectory')
print -depsc xz_trajectory

h3 = figure(3);
set(h3,'Visible', 'off'); 
hold on
plot(data(:,1), data(:,5))
plot(data(:,1), data(:,6),'k')
plot(data(:,1), data(:,7),'r')
xlabel('t')
ylabel('speed')
ylim([-2 2])
title('Velocity')
legend('vx', 'vy', 'vz')
print -depsc velocity

h4 = figure(4);
set(h4,'Visible', 'off'); 
hold on
plot(data(:,1), data(:,14))
plot(data(:,1), data(:,15),'k')
plot(data(:,1), data(:,16),'r')
plot(data(:,1), 1*ones(length(data(:,1)),1),'-.')
plot(data(:,1), -1*ones(length(data(:,1)),1),'-.')
xlabel('t')
ylabel('u')
ylim([-1.5 1.5])
title('Control')
legend('ux', 'uy', 'uz')
print -depsc control

h5 = figure(5);
set(h5,'Visible', 'off'); 
plot(data(:,1), data(:,17))
xlabel('time')
ylabel('Hamiltonian')
title('Hamiltonian')
print -depsc hamiltonian

