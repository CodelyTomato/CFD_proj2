clear all; close all; clc; 

rho_max = 0.1; 
v_max   = 50; 
dx      = 100; 
cars    = 100; 

L_traffic     = cars / rho_max;       % 100 / 0.1 = 1000 ft
cells_traffic = L_traffic / dx;       % 1000 / 50 = 20 cells

L_max   = 5000; 
M       = L_max / dx;                 % 100 cells
t_final = 60; 

schemes = {'Exact Riemann'}; 
x       = linspace(0, L_max, M+1); 

rho_0          = zeros(M+1, 1); 
rho_0(1:cells_traffic) = rho_max;     % 20 jam cells

SF = 0.9; 

results = zeros(length(schemes), M+1);

for i = 1:length(schemes)
    scheme = schemes{i}; 
    rho = rho_0; 
    t = 0; 

    while t < t_final
        a_vals = abs(a(rho, rho_max, v_max));     % use |a|
        dt_max = dx / max(a_vals); 
        dt = SF * dt_max; 

        if t + dt > t_final
            dt = t_final - t;
        end

        F = zeros(M+1, 1); 

        for j = 1:M-1
            rho_L = rho(j); 
            rho_R = rho(j+1); 

            switch scheme
                case 'Exact Riemann'
                    F(j+1) = ExactRiemann(rho_L, rho_R, rho_max, v_max); 
            end 
        end 

        % boundary fluxes
        F(1)   = 0;                             % no inflow upstream
        F(M+1) = flux(rho(M), rho_max, v_max);  % outflow on the right

        rho_new      = rho; 
        rho_new(2:M) = rho(2:M) - (dt/dx) * (F(3:M+1) - F(2:M));

        rho = rho_new; 
        t   = t + dt; 
    end 

    results(i,:) = rho.';   % store as row for plotting
end 

figure; hold on;
for n = 1:length(schemes)
    plot(x, results(n,:), 'LineWidth', 2, 'DisplayName', schemes{n}); 
end
grid on; 
xlabel('x (ft)'); 
ylabel('\rho (vehicles/ft)'); 
legend show;



function F = ExactRiemann(rho_L, rho_R, rho_max, v_max)
    a_L = a(rho_L, rho_max, v_max); 
    a_R = a(rho_R, rho_max, v_max); 
    
    if (a_L <= 0) && (a_R >= 0)
        F = (rho_max * v_max) / 4; 
    elseif (a_L + a_R) >= 0
        F = flux(rho_L, rho_max, v_max); 
    else
        F = flux(rho_R, rho_max, v_max); 
    end 
end 

function a_speed = a(rho, rho_max, v_max)
    a_speed = v_max * (1 - 2 * (rho / rho_max));
end 

function q = flux(rho, rho_max, v_max)
    q = v_max * (rho .* (1 - rho / rho_max)); 
end
