clear; clc; 
%{
Aidan Moriarty 
CFD Project 2: Numerical Advection
Program for part 3
%}

% -- Parameters -- 
M = 32; 
L = 1; 
a = 1; 
dx = L/M; 
T = 0.5; 
nu_values = [1, 0.75, 0.5, 0.25]; 

schemes = {'Lax-Friedrichs', 'First-Order Upwind', 'Lax-Wendroff', 'Minimum Dispersion'}; 

x = linspace(0, L, M+1); 
IC = zeros(M+1, 1); 
IC(1:M/4 + 1) = 1; 

L1_errors = zeros(1,length(nu_values)); 
 
u_exact = zeros(M+1, 1); 
x_exact_shock = L/4 + a*T; 
% -- Main Function -- 
for j = 1:M+1
    if x(j) <= x_exact_shock
        u_exact(j) = 1; 
    end 
end 
figure; 
hold on; 
plot(x, u_exact, 'LineWidth', 2, 'DisplayName', 'Exact Solution'); 
for j = 1:length(nu_values)
    nu = nu_values(j); 
    u_final = solve_advection_shock_flux (M, dx, T, a, nu, IC); 

    P = M+1; 
    L1_error_norm = (1/P) * sum(abs(u_final - u_exact)); 
    L1_errors(j) = L1_error_norm; 
    plot(x, u_final, 'LineWidth', 2, 'DisplayName', sprintf('|\\nu| = %.2f', nu, L1_error_norm)); 
end 
xlabel('x'); 
ylabel('u(x,T)'); 
legend show; 
grid on; 

% -- Outputs -- 
fprintf('\n## L1 Error Norm Matrix for Shock Problem ##\n');
nu_string_labels = arrayfun(@(n) sprintf('nu_%.2f', n), nu_values, 'UniformOutput', false);

error_table = array2table(L1_errors, ...
    'VariableNames', nu_string_labels); 
disp(error_table); 

figure; 
hold on; 
xlabel('\nu'); 
ylabel('L1 Error Norm'); 
grid on; 


plot(nu_values, L1_errors(1,:), 'LineWidth', 2, 'DisplayName', 'Modified FLux'); 
legend show; 

% Calculates the linear advection given the new flux formula
function u_final = solve_advection_shock_flux(M, dx, T, a, nu, IC)
    dt = nu*dx; 
    nt = floor(T / dt); 
    u = zeros(M+1, nt+1); 
    u(:,1) = IC; 
    for n = 1:nt
        u(1, n+1) = 1; 
        for j = 2:M
            F_right = calc_F_right(u, j, n, a, nu); 
            F_left = calc_F_left(u, j, n, a, nu); 
            
            u(j,n+1) = u(j,n) - (dt/dx) * (F_right - F_left);
        end 

        % Neumann BC at x = 1 
        u(M+1, n+1) = u(M, n+1); 

    end 
    % outputs u(x,T)
    u_final = u(:, nt+1); 
end 

function F_right = calc_F_right(u, j, n, a, nu)
    ujm1 = u(j-1, n); 
    uj = u(j,n); 
    ujp1 = u(j+1, n); 

    F_right = calc_flux(uj, ujp1, ujm1, a, nu); 
end

function F_left = calc_F_left(u, j, n, a, nu)
    if j == 2
        ujm1 = u(1,n); 
    else
        ujm1 = u(j-2,n); 
    end 
    uj = u(j-1,n); 
    ujp1 = u(j,n);

    F_left = calc_flux(uj, ujp1, ujm1, a, nu); 
end

function F = calc_flux(uj, ujp1, ujm1, a, nu)
    x = ujp1 - uj; 
    y = uj - ujm1; 

    B = get_B(x,y);

    u = uj + 0.5*(1-nu)*B; 
    F = a * u; 
end

% New function to calculate the B term 
function B = get_B(x,y)
    if x * y < 0
        B = 0; 
    else
        if (0.5 <= x/y) && (x/y <= 2)
            B_abs = max(abs(x), abs(y)); 
        else
            B_abs = 2*min(abs(x), abs(y)); 
        end 
        B = sign(x) * B_abs; % assign sign of x 
    end 
end




            