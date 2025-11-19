clear; clc; 
%{
Aidan Moriarty
CFD Project 2: Numerical Advection
Program for part 3 - shock discontinuity 
%}

% -- Parameters -- 
M = 32; 
L = 1; 
a = 1; 
dx = L/M; 
T_values = [0.5, 1.0]; 
nu_values = [1, 0.75, 0.5, 0.25]; 

schemes = {'Lax-Friedrichs', 'First-Order Upwind', 'Lax-Wendroff', 'Minimum Dispersion'}; 
x = linspace(0, L, M+1); 
IC = zeros(M+1, 1); 
IC(1:M/4 + 1) = 1; % shock discontinuity initial condition 
L1_errors = zeros(length(schemes), length(nu_values), length(T_values)); 

% -- Main Function -- 
for t = 1:length(T_values)
    T = T_values(t); 
    u_exact = zeros(M+1, 1); 
    x_exact_shock = L/4 + a*T;
    % initial condition 
    for j = 1:M+1
        if x(j) <= x_exact_shock
            u_exact(j) = 1; 
        end
    end 
    for i = 1:length(schemes)
        scheme = schemes{i}; 
        figure; 
        hold on; 
        plot(x, u_exact, 'LineWidth', 2, 'DisplayName', 'Exact Solution'); 
        for j = 1:length(nu_values)
            nu = nu_values(j); 
            u_final = solve_advection_shock (M, dx, T, a, scheme, nu, IC); 
            P = M+1; 
            L1_error_norm = (1/P) * sum(abs(u_final - u_exact)); 
            L1_errors(i,j,t) = L1_error_norm; 
            plot(x, u_final, 'LineWidth', 2, 'DisplayName', sprintf('|\\nu| = %.2f', nu)); 
        end 
        xlabel('x'); 
        ylabel('u(x,T)'); 
        legend show; 
        grid on; 
    end 
end 

% -- Outputs -- 
fprintf('\n## L1 Error Norm Matrix for Shock Problem ##\n');
nu_string_labels = arrayfun(@(n) sprintf('nu_%.2f', n), nu_values, 'UniformOutput', false);
 for t = 1:length(T_values)
    T = T_values(t);
    error_table = array2table(L1_errors(:,:,t), 'RowNames', schemes, ...
        'VariableNames', nu_string_labels); 
    disp(error_table); 
    
    figure; 
    hold on; 
    xlabel('\nu'); 
    ylabel('L1 Error Norm'); 
    grid on; 

    for i = 1:length(schemes)
        plot(nu_values, L1_errors(i,:,t), 'LineWidth', 2, 'DisplayName', schemes{i}); 
    end 
    legend show; 
 end 

%{ 
Calculates the linear advection given the shock discontinuity 
%}
function u_final = solve_advection_shock(M, dx, T, a, scheme, nu, IC)
    dt = nu*dx; 
    nt = floor(T / dt); 
    u = zeros(M+1, nt+1); 
    u(:,1) = IC; 
    q = get_q(scheme, nu); 
    for n = 1:nt
        u(1, n+1) = 1; % Dirichlet condition
        for j = 2:M
            jp1 = j+1; 
            jm1 = j-1; 
            
            F_right = calc_flux(u(j,n), u(jp1,n), a, dx, dt, q); 
            F_left = calc_flux(u(jm1,n), u(j,n), a, dx, dt, q); 
            
            u(j,n+1) = u(j,n) - (dt/dx) * (F_right - F_left);
        end 
        u(M+1, n+1) = u(M, n+1); % Neumann condition
    end 
    % returns u(x,T)
    u_final = u(:, nt+1); 
end 

function q = get_q(scheme, nu)
    switch lower(scheme)
        case 'lax-friedrichs'
            q = 1; 
        case 'first-order upwind'
            q = abs(nu); 
        case 'lax-wendroff'
            q = nu^2; 
        case 'minimum dispersion'
            q = 1/3 + (2/3)*nu^2; 
        otherwise
            error('invalid scheme'); 
    end 
end 

function F = calc_flux(u_L, u_R, a, dx, dt, q)
    f_L = a*u_L; 
    f_R = a*u_R; 

    F = 0.5*(f_L + f_R) - (q*dx)/(2*dt)*(u_R - u_L);  
end 






            