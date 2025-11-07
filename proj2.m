clear; clc; 
M = 32; 
L = 1; 
a = 1; 
T = 3/2; 

nu_values = [1, 0.75, 0.5, 0.25]; 
schemes = {'Lax-Friedrichs', 'First-Order Upwind', 'Lax-Wendroff', 'Minimum Dispersion'}; 

x = linspace(0, L, M+1);
u_exact = sin(2*pi*(x-a*T))'; 

L1_errors = zeros(length(schemes), length(nu_values)); 

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
        plot(x,u_final, 'LineWidth', 2, 'DisplayName', sprintf('|\\nu| = %.4f', nu)); 
        
    end 
    xlabel('x'); 
    ylabel('u(x,T)'); 
    legend show; 
    grid on; 
end 

nu_string_labels = arrayfun(@(n) sprintf('nu_%.2f', n), nu_values, 'UniformOutput', false);
error_table = array2table(L1_errors, 'RowNames', schemes, 'VariableNames', nu_string_labels); 
fprintf('\n### L1 Error Norms ###\n');
disp(error_table); 

figure; 
hold on; 
for i = 1:length(schemes)
    plot(nu_values, L1_errors(i,:), '-o', 'LineWidth', 2, 'DisplayName', schemes{i}); 
end 
xlabel('\nu'); 
ylabel('L1 Error Norm'); 
legend show; 
grid on; 

function u_final = solve_advection(a, M, x, T, L, scheme, nu)

    dx = L/M; 
    nt = round(3 / (2*nu*dx)); 
    dt = T / nt;  
    
    u = zeros(M+1, nt+1); 
    u(:,1) = sin(2*pi*x)'; 
    q = get_q(scheme, nu); 

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

            F_right = calc_flux(u(j,n), u(jp1,n), a, dx, dt, q); 
            F_left = calc_flux(u(jm1,n), u(j,n), a, dx, dt, q); 
            
            u(j,n+1) = u(j,n) - (dt/dx) * (F_right - F_left); 
        end 
        u(M+1, n+1) = u(1, n+1); 
    end 
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




        