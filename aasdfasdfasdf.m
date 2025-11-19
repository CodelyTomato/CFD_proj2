dx = 50; % ft 
rho_max = 0.1; % vehicles / ft
v_max = 50; % ft/s
Nx = 100; % number of cells
x = 0:dx:(Nx-1)*dx; % grid
T_final = 60; % s 

% --- Initial Condition --- % 
rho_0 = zeros(1, Nx); 
% initial line is 100 vehicles, which takes up 1000 ft of road. Using dx =
% 50 ft, this corresponds to 20 cells. 
rho_0(1:20) = rho_max; % first 20 cells at max density 

safety = [0.9, 0.6, 0.3]; 
schemes = {'Exact Riemann','Lax-Wendroff','Roe','Roe-Entropy'};
results = zeros(length(schemes), length(SFs), Nx);


for i = 1:length(safety)
    sf = safety(i); 
    for j = 1:length(schemes)
        scheme = schemes{j};
        rho = rho_0; 
        t = 0; 

        while t < t_final
            % calculate dt

            % calculate interface flux
            F = zeros(1, Nx+1); 

            for k = 1:Nx-1
                rho_L = rho(k); 
                rho_R = rho(k+1); 
                switch scheme
                    case 'Exact Reimann'
                        % scheme
                    case 'Lax-Wendroff'
                        % scheme
                    case 'Roe'
                        % scheme
                    case 'Roe with Entropy Fix'
                        % scheme
                end 
            end 

            %{
F(1) = flux(0)
F(Nx+1) = flux(rho(Nx)); 
rho_new = rho; 
rho_new(2:Nx-1) = rho(2:Nx-1) - (dt/dx)*(F(3:Nx) - F(2:Nx-1)); 

rho = rho_new; 
t = t + dt; 
end 

results(k,s,:) = rho; 
            %}


        