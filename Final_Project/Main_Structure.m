close all;clear;clc
% Physical Parameters
R_i = 10; R_o = 20;
rou = 0.01; E = 100; nu = 0.3;
lambda = E*nu/((1 + nu)*(1 - 2*nu));
mu = E/(2*(1 + nu));
% Boundary Condition
P0 = 1.0;
% Spices (2 Node Linear)
x_iso = @(x_e,xi) [0.5*(1-xi), 0.5*(1+xi)]*x_e;  % Isoparametric Mapping, plug this as the last argument in B
N = @(xi) [0.5*(1-xi) 0.5*(1+xi)];
D = [lambda + 2*mu, lambda, lambda; lambda, lambda + 2*mu, lambda; lambda, lambda, lambda + 2*mu];
B = @(h,xi,x_iso) [-1/h, 1/h; (1-xi)/(2*x_iso), (1+xi)/(2*x_iso); (1-xi)/(2*x_iso), (1+xi)/(2*x_iso)];
% % Two Point Gauss Quadrature
% xi_1 = -sqrt(1/3); xi_2 = sqrt(1/3); xi = [xi_1,xi_2];
% w1 = 1; w2 = 1; w = [w1,w2];
% Three Poing Gauss Quadrature
xi = [-sqrt(3/5), 0, sqrt(3/5)];
w = [5/9, 8/9, 5/9];
% Parameters of NewMark Method
beta_all = [0, 0.25];
gamma_all = [0.5, 0.5];
% Initialize a Cell to store all the results
% First Row: Central Difference with Lumped Mass
% Second Row: Average Acceletration with Consistent Mass
sol_disp = cell(2,4);
sol_sigma_rr = cell(2,4);
sol_sigma_pp = cell(2,4);
nof_element_all = [10,20,40,80];
% Location of points for plots
location = [13,15,17];
dt_all = zeros(1,4);

for method = 1 : 2
    beta = beta_all(method);
    gamma = gamma_all(method);
    for nelem = 1 : 4
        %% Discretization
        nof_element = nof_element_all(nelem);
        x = linspace(R_i,R_o, nof_element+1).';
        h = (R_o - R_i)/nof_element;
       
        %% Stiffness Matrix and Mass Matrix
        M = zeros(nof_element+1, nof_element+1);
        K = zeros(nof_element+1, nof_element+1);
        for i = 1 : nof_element  % Scanning through all elements
            x_e = x(i:i+1);
            M_e = zeros(2,2);
            K_e = zeros(2,2);
            % Gauss Integration
            for gp = 1 : 3
                N_e = N(xi(gp));
                x_e_iso = x_iso(x_e, xi(gp));
                M_e = M_e + w(gp)*2*pi*h*rou*N_e.'*N_e*x_e_iso^2;  % Element Mass Matrix
                B_e = B(h, xi(gp),x_e_iso);
                K_e = K_e + w(gp)*2*pi*h*B_e.'*D*B_e*x_e_iso^2;  % Element Stiffness Matrix
            end
            % Global Assembly
            M(i:i+1,i:i+1) = M(i:i+1,i:i+1) + M_e;
            K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + K_e;
        end
        if method == 1
            % Convert Consistent Mass to Lumped Mass by Row Sum Method
            M = sum(M,2);
            M = diag(M);
        end
        
        % Critical Time Step, only compute this when using central difference
        if method == 1
            dt_critical = 2/sqrt(max(eig(K,M)));
            dt = 0.5*dt_critical
            dt_all(nelem) = dt;
        end
        
        if method ==2
            dt = 0.001
        end
        
        % Temportal Discretization
        t_total = 5;
        t = 0:dt:t_total;
        nof_time_step = length(t)-1;
  
        %% Prelocations
        disp_mat = zeros(nof_element+1, nof_time_step+1); % Store d_n at all location at all time
        sigma_rr = zeros(3, nof_time_step+1); % Store radial stress at three different positions
        sigma_pp = zeros(3, nof_time_step+1); 
        
        %% Initial Condition
        v_n = zeros(nof_element+1,1);
        a_n = zeros(nof_element+1,1);
        d_n = zeros(nof_element+1,1);
        disp_mat(:,1) = d_n;

        %% Time Marching
        for j = 1 : nof_time_step
            f = zeros(nof_element+1, 1);
            R_current = R_i + d_n(1);
            
            f(1) = 4*pi*R_current^2*P0*(R_i/R_current)^3.75;
            % Displacement, Velocity and Acceletration at new time step; Which serves as d_n, v_n, a_n in the next loop
            M_star = M + beta*dt^2*K;
            f_star = f - K*(d_n + dt*v_n + 0.5*(1-2*beta)*dt^2*a_n);
            a_previous = a_n;
            
            a_n = M_star\f_star;
            d_n = d_n + dt*v_n + 0.5*dt^2*((1-2*beta)*a_previous + 2*beta*a_n);
            v_n = v_n + dt*((1-gamma)*a_previous + gamma*a_n);
            
            disp_mat(:,j+1) = d_n;
            
            % Do interpolation to get stressess at three location
            for s = 1 : 3
                index = find(x==location(s));
                x_e = x(index:index+1);
                x_e_iso = x_iso(x_e, -1);
                B_e = B(h, -1, x_e_iso)*d_n(index:index+1);
                sigma_interp = D*B_e;  % 3 by 1 vector which gives three components of stressess
                sigma_rr(s, j) = sigma_interp(1);
                sigma_pp(s,j) = sigma_interp(2);
            end              
        end
        sol_disp{method, nelem} = disp_mat;
        sol_sigma_rr{method, nelem} = sigma_rr;
        sol_sigma_pp{method, nelem} = sigma_pp;
    end
end

%%Plots for Central Differnece with Lumped Mass
disp_string = {'$u_{r} \ at \ r = 13$','$u_{r} \ at \ r = 15$','$u_{r} \ at \ r = 17$'};
sigma_rr_string = {'$\sigma_{rr} \ at \ r = 13$','$\sigma_{rr} \ at \ r = 15$','$\sigma_{rr} \ at \ r = 17$'};
sigma_pp_string = {'$\sigma_{\phi\phi} \ at \ r = 13$','$\sigma_{\phi\phi} \ at \ r = 15$','$\sigma_{\phi\phi} \ at \ r = 17$'};
linestyle_string = {'k:','g-.','r--','b-'};
fig_num = reshape(1:1:9,3,3).';
for loc = 1 : 3
    figure(fig_num(loc,1)); hold on; xlabel('time'); ylabel('Displacement');
    title(disp_string{loc}, 'interpreter','latex','fontweight','bold');
    
    figure(fig_num(loc,2)); hold on; xlabel('time'); ylabel('Radial Stress')
    title(sigma_rr_string{loc}, 'interpreter','latex','fontweight','bold');
    
    figure(fig_num(loc,3)); hold on; xlabel('time'); ylabel('Radial Stress')
    title(sigma_pp_string{loc}, 'interpreter','latex','fontweight','bold');
    for p = 1 : 4
        x = linspace(R_i,R_o,nof_element_all(p)+1);
        t = 0:dt_all(1,p):5;
        size(t)
        
        index = find(x==location(loc));
        figure(fig_num(loc,1));
        plot(t,sol_disp{1,p}(index, :), linestyle_string{p});
        
        
        figure(fig_num(loc,2));
        plot(t, sol_sigma_rr{1,p}(1,:), linestyle_string{p});
        
        figure(fig_num(loc,3));
        plot(t, sol_sigma_pp{1,p}(3,:), linestyle_string{p});
    end
end


%%Plots for average acceleration
disp_string = {'$u_{r} \ at \ r = 13$','$u_{r} \ at \ r = 15$','$u_{r} \ at \ r = 17$'};
sigma_rr_string = {'$\sigma_{rr} \ at \ r = 13$','$\sigma_{rr} \ at \ r = 15$','$\sigma_{rr} \ at \ r = 17$'};
sigma_pp_string = {'$\sigma_{\phi\phi} \ at \ r = 13$','$\sigma_{\phi\phi} \ at \ r = 15$','$\sigma_{\phi\phi} \ at \ r = 17$'};
linestyle_string = {'k:','g-.','r--','b-'};
fig_num = reshape(10:1:18,3,3).';
for loc = 1 : 3
    figure(fig_num(loc,1)); hold on; xlabel('time'); ylabel('Displacement');
    title(disp_string{loc}, 'interpreter','latex','fontweight','bold');
    
    figure(fig_num(loc,2)); hold on; xlabel('time'); ylabel('Radial Stress')
    title(sigma_rr_string{loc}, 'interpreter','latex','fontweight','bold');
    
    figure(fig_num(loc,3)); hold on; xlabel('time'); ylabel('Radial Stress')
    title(sigma_pp_string{loc}, 'interpreter','latex','fontweight','bold');
    for p = 1 : 4
        x = linspace(R_i,R_o,nof_element_all(p)+1);
        t = 0:0.001:5;
        size(t);
        
        index = find(x==location(loc));
        figure(fig_num(loc,1));
        plot(t,sol_disp{2,p}(index, :), linestyle_string{p});
        
        
        figure(fig_num(loc,2));
        plot(t, sol_sigma_rr{2,p}(1,:), linestyle_string{p});
        
        figure(fig_num(loc,3));
        plot(t, sol_sigma_pp{2,p}(3,:), linestyle_string{p});
    end
end


