clear; clc; 
%{
Aidan Moriarty 
CFD Project 2: Numerical Advection
Program for part 2
%}
% -- Parameters -- 
M = 32; 
L = 1; 
a = 1; 
T = 3/2; 
total_shift = 2*pi*T; 
nu_values = [1, 0.75, 0.5, 0.25]; 
x = linspace(0, L, M+1);
u_exact = sin(2*pi*(x-a*T))'; 

schemes = {'Lax-Friedrichs', 'First-Order Upwind', 'Lax-Wendroff', 'Minimum Dispersion'}; 
L1_errors = zeros(length(schemes), length(nu_values)); 
phase_errors = zeros(length(schemes), length(nu_values)); 

% -- Main Function -- 
for i = 1:length(schemes)
    scheme = schemes{i}; 
    figure;
    hold on; 
    plot(x, u_exact, 'LineWidth', 2, 'DisplayName', 'Exact Solution'); 
    for j = 1:length(nu_values)
        nu = nu_values(j);
        u_final = solve_advection(a, M, x, T, L, scheme, nu); 
        P = M + 1; 
        L1_error_norm = (1/P)* sum(abs(u_final - u_exact)); 
        L1_errors(i,j) = L1_error_norm; 

        x_exact_zero = find_zero(x, u_exact, L/2); 
        x_num_zero = find_zero(x, u_final, L/2); 

        rad_error = 2*pi*(x_num_zero - x_exact_zero); 
        percent_error = (abs(rad_error) / total_shift) * 100; 
        if nu == 1
            phase_errors(i,j) = 0; 
        else
  
            phase_errors(i,j) = percent_error; 
        end 
        
        plot(x,u_final, 'LineWidth', 2, 'DisplayName', sprintf('|\\nu| = %.4f', nu)); 
        
    end 
    xlabel('x'); 
    ylabel('u(x,T)'); 
    legend show; 
    grid on; 
end 

% -- Outputs -- 
nu_string_labels = arrayfun(@(n) sprintf('nu_%.2f', n), nu_values, 'UniformOutput', false);
error_table = array2table(L1_errors, 'RowNames', schemes, 'VariableNames', nu_string_labels); 
fprintf('\n### L1 Error Norms ###\n');
disp(error_table); 

phase_table = array2table(phase_errors, 'RowNames', schemes, 'VariableNames', ...
    nu_string_labels); 
fprintf('\n### Phase Errors ###\n')
disp(phase_table); 

figure; 
hold on; 
for i = 1:length(schemes)
    plot(nu_values, L1_errors(i,:), '-o', 'LineWidth', 2, 'DisplayName', schemes{i}); 
end 
xlabel('\nu'); 
ylabel('L1 Error Norm'); 
legend show; 
grid on; 

%{
Solves the linear advection equation for part 2. 
%}
function u_final = solve_advection(a, M, x, T, L, scheme, nu)

    dx = L/M; 
    nt = round(3 / (2*nu*dx)); 
    dt = T / nt;  
    
    u = zeros(M+1, nt+1); 
    u(:,1) = sin(2*pi*x)'; 
    q = get_q(scheme, nu); 
    % Periodic boundary condition u(0,t) = u(1,t)
    for n=1:nt
        for j=1:M
            if j == M
                jp1 = 1; 
            else
                jp1 = j + 1; 
            end 
            if j == 1
                jm1 = M; 
            else
                jm1 = j - 1; 
            end 
            % Calcuation of flux interface
            F_right = calc_flux(u(j,n), u(jp1,n), a, dx, dt, q); 
            F_left = calc_flux(u(jm1,n), u(j,n), a, dx, dt, q); 
            
            u(j,n+1) = u(j,n) - (dt/dx) * (F_right - F_left); 
        end
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

%{
Function to linearlly interpolate to find where u(x,T) crosses 0. 
%}
function x0 = find_zero(x, u, center)

    M = length(x); 
    min_dist = inf; 
    x0 = NaN; 
    
    for j = 1:(M-1)
        if u(j) * u(j+1) < 0 
            x1 = x(j); 
            x2 = x(j+1); 
            u1 = u(j); 
            u2 = u(j+1); 

            dx = x2 - x1; 
            du = u2 - u1; 

            current_x0 = x1 - u1 * (dx / du); 

            current_dist = abs(current_x0 - center); 
            if current_dist < min_dist
                min_dist = current_dist; 
                x0 = current_x0; 
            end 
        end 
    end
end


        