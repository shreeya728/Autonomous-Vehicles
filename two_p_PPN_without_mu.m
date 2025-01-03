function [dy_dt, cost] = two_p_PPN_without_mu(t, y, V_P, V_T, alpha_P_df, alpha_P0, theta_0, theta_d, reset_flag)
    %N = 3;
    R = y(1);
    theta = y(2);
    V_theta = y(3);
    V_R = y(4);
    alpha_P = y(5);
    alpha_T = y(6);
    x_T = y(7);
    y_T = y(8);
    x_P = y(9);
    y_P = y(10);

    persistent t_prev cost_acc
    if reset_flag
        t_prev = 0;
        cost_acc = 0;
    elseif isempty(t_prev)
        t_prev = 0;
        cost_acc = 0;
    end

    mu_th = pi/4;
    mu = alpha_P - theta;


    %N = (alpha_P_df - alpha_P)/(alpha_P_df - theta);
    N = 2;
    if (alpha_P_df - alpha_P)/(alpha_P_df - theta) < 2
        %N = (2/pi)*(alpha_P0 - theta_d);
        %N = (theta_d - alpha_P0)/(theta_d - (pi/2));
        N = 2*(theta_d + theta_0 - alpha_P0)/(theta_d - theta_0);
    end

    dy_dt = zeros(10, 1);

    dy_dt(1) = V_T*cos(alpha_T - theta) - V_P*cos(mu);  
    dy_dt(2) =  (V_T*sin(alpha_T - theta) - V_P*sin(mu))/R;       
    aP = N * V_P * V_theta/R;

    cost_acc = cost_acc + aP^2 * (t - t_prev);
    t_prev = t;
    cost = cost_acc;
 
    dy_dt(5) = aP/V_P;                          
    dy_dt(6) = 0;  

    dy_dt(4) = -V_T*sin(alpha_T - theta)*(dy_dt(6) - V_theta/R) + V_P*sin(alpha_P - theta)*(dy_dt(5) - V_theta/R);      
    dy_dt(3) = V_T*cos(alpha_T - theta)*(dy_dt(6) - V_theta/R) - V_P*cos(alpha_P - theta)*(dy_dt(5) - V_theta/R);  

    dy_dt(7) = V_T*cos(alpha_T);                
    dy_dt(8) = V_T*sin(alpha_T);                
    dy_dt(9) = V_P*cos(alpha_P);                
    dy_dt(10) = V_P*sin(alpha_P);                 
end
