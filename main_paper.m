clc;
clear;
close all;


global cost ap_matrix scenario_num t_span;

R_0 = 2500;
theta_0 = 0*pi/180;
x_p0 = 0;
y_p0 = 0;
x_t0 = R_0*cos(theta_0);
y_t0 = R_0*sin(theta_0);

%% TIME CONDITIONS
t_step = 0.1;
t_end = 500;
t_span = 0:t_step:t_end;
t_terminate = 5;
options = odeset('Events', @(t, y) event_terminal(t, y));
%% PPN, 2pPPN and BPPN
%N = 3;
alpha_P_df = deg2rad(-150);
V_T = 0;
V_P = 50;
alpha_P0 = pi/4;
alpha_T0 = 0;
V_R0 = V_T*cos(alpha_T0 - theta_0) - V_P*cos(alpha_P0 - theta_0);
V_theta_0 = V_T*sin(alpha_T0 - theta_0) - V_P*sin(alpha_P0 - theta_0);
y0 = [R_0, theta_0, V_theta_0, V_R0, alpha_P0, alpha_T0, x_t0, y_t0, x_p0, y_p0];
%[t,y] = ode45(@(t,y) PPN_paper(t, y, V_P, V_T, alpha_P_df), t_span, y0, options);
%[t,y] = ode45(@(t,y) two_p_PPN_without_mu(t, y, V_P, V_T, alpha_P_df, alpha_P0, theta_0), t_span, y0, options);
%[t,y] = ode45(@(t,y) BPPN(t, y, V_P, V_T, alpha_P_df, alpha_P0, theta_0, V_theta_0, R_0), t_span, y0, options);

%theta_d_values = linspace(alpha_P0, alpha_P_df, 10); 
theta_d_values = linspace(deg2rad(-90),deg2rad(-180), 10);  
results = cell(1, length(theta_d_values));   
costs = zeros(1, length(theta_d_values));        


N0 = (alpha_P_df - alpha_P0)/(alpha_P_df - theta_0)

N_ori = linspace(0,N0,14);
scenario_count = length(N_ori); 
ap_matrix = nan(length(t_span),scenario_count);

for scenario_num = 1:1:1 %length(theta_d_values)
    theta_d = theta_d_values(scenario_num);   
    %reset_flag = (i == 1);
    [t, y] = ode45(@(t, y) two_p_PPN_without_mu2(t, y, V_P, V_T, alpha_P_df, alpha_P0, theta_0, N_ori(scenario_num)), t_span, y0, options);
    results{scenario_num} = [t, y];
    disp(cost);
    costs(scenario_num) = cost;
    cost = 0;
    % Extract cost after solving the ODE
    %[~, costs(i)] = two_p_PPN_without_mu(0, y0, V_P, V_T, alpha_P_df, alpha_P0, theta_0, theta_d, reset_flag);  % Reset cost
end


% Plot cost vs theta_d
figure;
plot(rad2deg(theta_d_values), costs, '-o', 'LineWidth', 1.5);
xlabel('\theta_d (deg)');
ylabel('Cost');
title('Cost vs \theta_d');
grid on;





%% PLOTS AND ANIMATION

% figure;
% plot(t, y(:, 1));
% xlabel('Time (s)');
% ylabel('R');
% 
% title('R over time');
% grid on;
% 
% figure;
% plot(t, (180/pi).*y(:, 2));
% xlabel('Time (s)');
% ylabel('\theta');
% title('\theta over time');
% grid on;
% 
% figure;
% plot(t, y(:, 3));
% xlabel('Time (s)');
% ylabel('V_{\theta}');
% title('V_{\theta} over time');
% grid on;
% 
% figure;
% plot(t, y(:, 4));
% xlabel('Time (s)');
% ylabel('V_R');
% title('V_R over time');
% grid on;
% 
% figure;
% plot(t, y(:, 5));
% xlabel('Time (s)');
% ylabel('\alpha_P');
% title('\alpha_P over time');
% grid on;
% 
% 
% figure;
% plot(y(:, 3), y(:, 4));
% xlabel('V_{\theta}');
% ylabel('V_R');
% title('V_R vs V_{\theta}');
% grid on;
% 
% % N = (alpha_P_df - y(:, 5))./(alpha_P_df - y(:, 2));
% % aP = (V_P/0.01).*diff(y(: , 5)); % PPN
% % 
% % theta_dot = y(:, 3)./y(:, 1);
% % figure;
% % plot(t, (180/pi).*theta_dot)
% % xlabel('t');
% % ylabel('$\dot{\theta}$', 'Interpreter', 'latex');
% % title('$\dot{\theta}$ vs t', 'Interpreter', 'latex');
% % grid on;
% % 
% % 
% % figure;
% % plot(t(1:end-1), aP)
% % xlabel('t');
% % ylabel('a_p');
% % title('a_p vs t');
% % grid on;
% % 
% % figure;
% % mu = (180/pi).*(y(:, 5) - y(:, 2));
% % plot(t, mu)
% % xlabel('t');
% % ylabel('\mu');
% % title('\mu vs t');
% % grid on;
% 
%  
% % Extract the trajectories
% x_T = y(:, 7); 
% y_T = y(:, 8); 
% x_P = y(:, 9); 
% y_P = y(:, 10); 

% %Plotting the initial positions
% 
% figure;
% hT = plot(x_T, y_T, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Target point
% hold on;
% hP = plot(x_P(1), y_P(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Pursuer point
% hTrajT = plot(x_T(1), y_T(1), 'b-', 'LineWidth', 1.5); % Target trajectory
% hTrajP = plot(x_P(1), y_P(1), 'r-', 'LineWidth', 1.5); % Pursuer trajectory
% xlabel('x');
% ylabel('y');
% title('Trajectories of Target (T) and Pursuer (P)');
% %legend('Target (T)', 'Pursuer (P)');
% grid on;
% axis equal;
% hold off;

% Animation loop
% for i = 1:length(t)
%     Update target position
%     set(hT, 'XData', x_T(i), 'YData', y_T(i));
%     Update pursuer position
%     set(hP, 'XData', x_P(i), 'YData', y_P(i));
%     Update trajectories
%     set(hTrajT, 'XData', x_T(1:i), 'YData', y_T(1:i));
%     set(hTrajP, 'XData', x_P(1:i), 'YData', y_P(1:i));
%     
%     Pause to control animation speed
%     pause(0.01);
% end
% 
% Extract the final positions and trajectories
% x_T_final = x_T(end); % Final x-position of the target
% y_T_final = y_T(end); % Final y-position of the target
% x_P_final = x_P(end); % Final x-position of the pursuer
% y_P_final = y_P(end); % Final y-position of the pursuer
% 
% Plot the initial and final positions
% figure;
% Target's initial position (stationary target)
% plot(x_T(1), y_T(1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); 
% hold on;
% Pursuer's initial position
% plot(x_P(1), y_P(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); 
% Target's final position (stationary)
% plot(x_T_final, y_T_final, 'b*', 'MarkerSize', 8); 
% Pursuer's final position
% plot(x_P_final, y_P_final, 'r*', 'MarkerSize', 8); 
% 
% Plot the trajectories
% plot(x_T, y_T, 'b-', 'LineWidth', 1.5); % Target trajectory (should be stationary)
% plot(x_P, y_P, 'r-', 'LineWidth', 1.5); % Pursuer trajectory
% 
% xlabel('x');
% ylabel('y');
% title('Final Positions and Trajectories of Target (T) and Pursuer (P)');
% grid on;
% axis equal;
% hold off;
