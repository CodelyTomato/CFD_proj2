clear; close all; clc; 

% don't forget to change safety factor to see what it does. 

% -- Parameters -- % 
rho_max = 0.1; 
v_max = 50; 
N_cars = 100; 
t_final = 60; 
SF = 0.3; 
schemes = {'Exact Riemann', 'Lax-Wendroff', 'Roe Approximate', 'Roe Approximate w/ Entropy Fix'}; 

%{
discont_at_interface: 
    - true: initial discontinuity at interface 
    - false: initial discontinuity in interior of a cell
%}
discont_at_interface = true; 

% -- Initialize Grid -- % 
L = 10000; 
dx = 50; 
x_min = -L / 2; 
x_max = L / 2; 
x_vec = x_min:dx:x_max;
M = length(x_vec); 

% -- Setup Initial Conditions -- % 
L_traffic = N_cars / rho_max; 
N_cells_jam = round(L_traffic/dx); 
rho_0 = zeros(M, 1); 
idx_0 = round((0 - x_min)/dx) + 1; 
if discont_at_interface
    idx_jam = find(x_vec >= -L_traffic & x_vec <= 0); 
else
    shift = -dx / 2; 
    idx_jam = find(x_vec >= -L_traffic & x_vec <= shift); 
end 
rho_0(idx_jam) = rho_max;

% -- Storage -- % 
x_first_all = cell(length(schemes), 1); 
t_first_all = cell(length(schemes), 1); 

x_last_all = cell(length(schemes), 1); 
t_last_all = cell(length(schemes), 1); 

rho_20_all = cell(length(schemes), 1); 
rho_40_all = cell(length(schemes), 1); 
rho_60_all = cell(length(schemes), 1); 


for i = 1:length(schemes)
    % -- initialize first car -- %
    t_first = 0; 
    x_first = 0; 
    v_first = 0; 
    
    % -- initialize last car -- %
    t_last = 0; 
    x_last = -1000; 
    v_last = 0; 

    v_light = []; 
    t_light = []; 

    % -- initialize scheme -- % 
    scheme = schemes{i}; 
    rho = rho_0; 
    t = 0; 

    while t < t_final
        % -- compute char speeds & max dt -- % 
        a_vals = a(rho, rho_max, v_max);     
        dt_max = dx / max(abs(a_vals)); 
        dt = SF * dt_max; 
        if t + dt > t_final
            dt = t_final - t;
        end
        
        % -- compute flux at this timestep -- % 
        F = zeros(M+1, 1); 
        for j = 1:M-1
            rho_L = rho(j); 
            rho_R = rho(j+1);

            switch scheme
                case 'Exact Riemann'
                    F(j+1) = ExactRiemann(rho_L, rho_R, rho_max, v_max); 
                case 'Lax-Wendroff'
                    F(j+1) = LaxWendroff(rho_L, rho_R, rho_max, v_max, dt, dx); 
                case 'Roe Approximate'
                    F(j+1) = RoeApproximate(rho_L, rho_R, rho_max, v_max); 
                case 'Roe Approximate w/ Entropy Fix'
                    F(j+1) = RoeApproximateFix(rho_L, rho_R, rho_max, v_max); 
            end 
        end 

        % -- boundary fluxes -- % 
        F(1) = 0;                            
        F(M+1) = flux(rho(M), rho_max, v_max);  

        % -- update density -- % 
        rho_new = rho; 
        rho_new(2:M) = rho(2:M) - (dt/dx) * (F(3:M+1) - F(2:M));
        rho = rho_new; % new density 

        % -- car tracking -- % 
        idx_car = floor((x_first(end) - x_min)/dx) + 1; 
        
        rho_local = rho(idx_car); 
        v_local = v_max * (1-rho_local/rho_max); 

        x_first(end+1) = x_first(end) + v_local*dt; 
        v_first(end+1) = v_local; 
        t_first(end+1) = t_first(end) + dt; 

        % -- last car path tracking -- % 
        idx_last = floor((x_last(end) - x_min)/dx) + 1; 
        rho_local_last = rho(idx_last); 
        v_local_last = v_max * (1-rho_local_last / rho_max); 

        x_last(end+1) = x_last(end) + v_local_last*dt; 
        v_last(end+1) = v_local_last; 
        t_last(end+1) = t_last(end) + dt; 

        if t >= 20 && isempty(rho_20_all{i})
            rho_20_all{i} = rho;  
        end

        if t >= 40 && isempty(rho_40_all{i})
            rho_40_all{i} = rho; 
        end 
        
        rho_light = rho(idx_0); 
        v_light(end+1) = v_max * (1 - rho_light / rho_max); 
        t_light(end+1) = t; 


        

        t = t + dt; 
    end 
    % -- store results for plotting -- % 

    x_first_all{i} = x_first; 
    t_first_all{i} = t_first;    

    x_last_all{i} = x_last; 
    t_last_all{i} = t_last; 

    rho_60_all{i} = rho; 
end 

rho_L = rho_max;   % 0.1
rho_R = 0;

a_L = a(rho_L, rho_max, v_max);   % should be -50
a_R = a(rho_R, rho_max, v_max);   % should be +50

t_ex = t_final;   % 60 sec
xi = x_vec / t_ex;

rho_exact = zeros(size(x_vec));

% Region 1: left constant
rho_exact(xi < a_L) = rho_L;

% Region 3: right constant
rho_exact(xi > a_R) = rho_R;

% Region 2: fan
fan = (xi >= a_L) & (xi <= a_R);
rho_exact(fan) = 0.5 * rho_max * (1 - xi(fan)/v_max);



% -- FIRST/LAST CAR PATHS PLOT -- % 
figure; hold on; grid on;
% first-car path
for i = 1:length(schemes)
    plot(x_first_all{i}, t_first_all{i}, 'LineWidth', 2);
end
% last-car path
for i = 1:length(schemes)
    plot(x_last_all{i}, t_last_all{i}, 'LineWidth', 2); % dashed for clarity
end
xlabel('x (ft)');
ylabel('t (s)');
legend_strings = [ ...
    strcat("First - ", schemes), ...
    strcat("Last - ", schemes) ...
];
legend(legend_strings, 'Location', 'bestoutside');

% -- DENSITY PLOT -- %
figure;
for i = 1:length(schemes)
    subplot(2,2,i); hold on; grid on;
    plot(x_vec, rho_20_all{i}, 'LineWidth', 2);
    plot(x_vec, rho_40_all{i}, 'LineWidth', 2);
    plot(x_vec, rho_60_all{i}, 'LineWidth', 2);
    xlabel('x [ft]');
    ylabel('\rho [vehicles/ft]');
    title(schemes{i});
    legend('t = 20 s', 't = 40 s', 't = 60 s', 'Location', 'best');
end

% -- VELOCITY AT LIGHT PLOT -- % 
figure; hold on; grid on; 
plot(t_light, v_light, 'LineWidth', 2); 
xlabel('Time [s]'); 
ylabel('Velocity [ft/s'); 
title('Velocity of Vehicles Passing the Light'); 

% -- EXACT SOLUTION AT t = 60 s PLOT -- % 
figure; hold on; grid on;
plot(x_vec, rho_exact, 'LineWidth', 2)
title('Exact Riemann Solution at t = 60 s')
xlabel('x'); ylabel('\rho')



% --- HELPER FUNCTIONS --- %
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

function F = LaxWendroff(rho_L, rho_R, rho_max, v_max, dt, dx)
    flux_right = flux(rho_R, rho_max, v_max); 
    flux_left = flux(rho_L, rho_max, v_max); 
    F = flux(0.5*(rho_L + rho_R) - 0.5*(dt/dx)*(flux_right - flux_left), rho_max, v_max); 
end

function F = RoeApproximate(rho_L, rho_R, rho_max, v_max)
    a_L = a(rho_L, rho_max, v_max); 
    a_R = a(rho_R, rho_max, v_max); 

    if (a_L + a_R >= 0)
        F = flux(rho_L, rho_max, v_max); 
    else
        F = flux(rho_R, rho_max, v_max); 
    end 
end

function F = RoeApproximateFix(rho_L, rho_R, rho_max, v_max)
    a_L = a(rho_L, rho_max, v_max); 
    a_R = a(rho_R, rho_max, v_max); 
    if (a_L <= 0) && (a_R >= 0)
        F = 0.5*(flux(rho_L, rho_max, v_max) + flux(rho_R, rho_max, v_max) - ...
            0.5*(a_R - a_L)*(rho_R - rho_L)); 
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




