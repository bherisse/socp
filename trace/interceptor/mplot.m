clear all
close all

R_Earth = 6378145;

color = ['b'; 'g'; 'r'];

for k=1:3

    file = ['trace_S' num2str(k) '.dat'];
    data = load(file);

    h1 = figure(1);
    set(h1,'Visible', 'off'); 
    hold on
    plot(data(:,6)*R_Earth/1000, data(:,2)/1000, 'LineWidth', 2, 'Color', color(k))
    xlabel('r_T \cdot L (km)', 'FontSize', 15)
    ylabel('altitude (km)', 'FontSize', 15)
    title('Latitude-altitude trajectories', 'FontSize', 15)
    grid on
    print -depsc Lh_trajectory

    h2 = figure(2);
    set(h2,'Visible', 'off');
    hold on
    plot(data(:,6)*R_Earth/1000, data(:,7)*R_Earth/1000, 'LineWidth', 2, 'Color', color(k))
    xlabel('r_T \cdot L (km)', 'FontSize', 15)
    ylabel('r_T \cdot l (km)', 'FontSize', 15)
    title('Latitude-longitude trajectories', 'FontSize', 15)
    grid on
    print -depsc Ll_trajectory

    h3 = figure(3);
    set(h3,'Visible', 'off'); 
    hold on
    plot(data(:,1), abs(data(:,14)), 'LineWidth', 2, 'Color', color(k))
    plot(data(:,1), 1*ones(length(data(:,1)),1),'-.')
    plot(data(:,1), -1*ones(length(data(:,1)),1),'-.')
    xlabel('t (s)', 'FontSize', 15)
    ylabel('|u|', 'FontSize', 15)
    ylim([-1.5 1.5])
    title('Norm of the control', 'FontSize', 15)
    grid on
    print -depsc control_u

    h4 = figure(4);
    set(h4,'Visible', 'off'); 
    hold on
    plot(data(:,1), data(:,15), 'LineWidth', 2, 'Color', color(k))
    xlabel('t (s)', 'FontSize', 15)
    ylabel('beta (rad)', 'FontSize', 15)
    title('Bank angle', 'FontSize', 15)
    grid on
    print -depsc control_beta

    h5 = figure(5);
    set(h5,'Visible', 'off');
    hold on
    plot(data(:,1), data(:,3), 'LineWidth', 2, 'Color', color(k))
    xlabel('t (s)', 'FontSize', 15)
    ylabel('speed (m/s)', 'FontSize', 15)
    grid on
    title('Velocities', 'FontSize', 15)
    print -depsc velocity

    h6 = figure(6);
    %set(h6,'Visible', 'off'); 
    hold on
    plot3(data(:,6)*R_Earth/1000, data(:,7)*R_Earth/1000, data(:,2)/1000, 'LineWidth', 2, 'Color', color(k))
    xlabel('r_T \cdot L (km)', 'FontSize', 15)
    ylabel('r_T \cdot l (km)', 'FontSize', 15)
    zlabel('altitude(km)', 'FontSize', 15)
    grid on
    title('Latitude-longitude-altitude trajectories', 'FontSize', 15)
    % print -depsc Llh_trajectory

    h7 = figure(7);
    set(h7,'Visible', 'off'); 
    hold on
    plot(data(:,1), data(:,14).*cos(data(:,15)), 'LineWidth', 2, 'Color', color(k))
    plot(data(:,1), 1*ones(length(data(:,1)),1),'-.')
    plot(data(:,1), -1*ones(length(data(:,1)),1),'-.')
    xlabel('t (s)', 'FontSize', 15)
    ylabel('u_1', 'FontSize', 15)
    ylim([-1.5 1.5])
    title('Vertical control', 'FontSize', 15)
    grid on
    print -depsc control_u1

    h8 = figure(8);
    set(h8,'Visible', 'off'); 
    hold on
    plot(data(:,1), data(:,14).*sin(data(:,15)), 'LineWidth', 2, 'Color', color(k))
    plot(data(:,1), 1*ones(length(data(:,1)),1),'-.')
    plot(data(:,1), -1*ones(length(data(:,1)),1),'-.')
    xlabel('t (s)', 'FontSize', 15)
    ylabel('u_2', 'FontSize', 15)
    ylim([-1.5 1.5])
    title('Horizontal control', 'FontSize', 15)
    grid on
    print -depsc control_u2
    
    h9 = figure(9);
    set(h9,'Visible', 'off'); 
    hold on
    plot(data(:,1), data(:,16), 'LineWidth', 2, 'Color', color(k))
    xlabel('t (s)', 'FontSize', 15)
    ylabel('H', 'FontSize', 15)
    title('Hamiltonian', 'FontSize', 15)
    grid on
    print -depsc hamiltonian

end