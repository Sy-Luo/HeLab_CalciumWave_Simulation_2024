clear;

% Simulation parameters
Lx = 1000; % Length of x-axis in micrometers (ignore in lattice model)
Ly = 1000; % Length of y-axis in micrometers
Nx = 50; % Number of lattice points in x direction
Ny = 50; % Number of lattice points in y direction
delta_x = Lx / Nx; % Interval in x direction (size of lattice unit)
delta_y = Ly / Ny; % Interval in y direction (size of lattice unit)
Nt = 12000; % Number of time steps in the simulation

% Number of random repeats (Nrp) and initialization of output variables
Nrp = 200;
or_repeat = 0.04; % Scaling factor for noise
mp_repeat = zeros(Nrp, 1); % Mean peak
mh_repeat = zeros(Nrp, 1); % Mean half-width
mT_repeat = zeros(Nrp, 1); % Mean period
vp_repeat = zeros(Nrp, 1); % Variance of peak
vh_repeat = zeros(Nrp, 1); % Variance of half-width
vT_repeat = zeros(Nrp, 1); % Variance of period

% Start random repeat loop
parfor krp = 1:Nrp
    krp % Output repeat number
    tscale = 0.4; % Time scale factor
    cscale = 1000; % Concentration scale factor
    akhscale = 0.0005; % Scaling factor for akh
    akhbase = 240; % Base concentration for akh
    sigma_c = or_repeat(krp) * cscale; % Noise term for concentration fluctuations

    % Reaction rate constants (units in μM^-1 s^-1 or s^-1)
    k1 = 400; % Forward reaction rate constant for reaction 1
    k1_minus = 52; % Reverse reaction rate constant for reaction 1
    k5 = 20; % Forward reaction rate constant for reaction 5
    k5_minus = 1.64; % Reverse reaction rate constant for reaction 5

    % Other constants
    gamma = 2; % Coefficient for membrane transport
    delta = 1; % Coefficient for SERCA pump
    
    % Variables for SERCA pump and membrane transport rates
    Vs = 0.9 * tscale * cscale; % Maximum SERCA pump rate (μM/s)
    Ks = 0.1 * cscale; % Half-saturation concentration of SERCA pump (μM)
    Vp = 0.05 * tscale * cscale; % Maximum membrane transport rate (μM/s)
    Kp = 0.3 * cscale; % Half-saturation concentration for membrane transport (μM)
    
    % Variables for membrane transport
    alpha1 = 0.07 * tscale * cscale * akhscale; % Scaling factor for Ca+ transport
    alpha0 = 0.01 * tscale * cscale + akhbase * alpha1; % Base transport rate
    
    % Reaction constants related to calcium dynamics
    kf = 1.11 * tscale; % Fast reaction rate constant (s^-1)
    kleak = 0.02 * tscale; % Membrane leakage rate constant (s^-1)
    K1 = k1_minus / k1 / akhscale + akhbase; % Combined constant for reaction 1
    K5 = k5_minus / k5 * cscale; % Combined constant for reaction 5
    
    % Constant y used in current flux function
    y = 0.3;
    
    % Define flux functions for calcium ions and inhibitors
    J_irp = @(c, ce, p) kf * ((p + akhbase) .* c .* (1 - y) ./ (p + K1) ./ (c + K5)).^3 .* (ce - c);
    J_leak = @(c, ce, p) kleak * (ce - c); % Leakage flux
    J_serca = @(c, ce, p) Vs * c.^2 ./ (Ks^2 + c.^2); % SERCA pump flux
    J_mem = @(c, ce, p) delta * (alpha0 + alpha1 * p - Vp * c.^2 ./ (Kp^2 + c.^2)); % Membrane flux

    % Rate of change functions for calcium and inhibitors
    dc = @(c, ce, p) J_irp(c, ce, p) + J_leak(c, ce, p) - J_serca(c, ce, p) + J_mem(c, ce, p);
    dce = @(c, ce, p) -gamma * (J_irp(c, ce, p) + J_leak(c, ce, p) - J_serca(c, ce, p));

    % Fixed point (steady-state concentration) for calcium and inhibitors
    fixpoint = [0.235, 15.6] * cscale;
    
    % Initialize concentration arrays for calcium (var_c) and inhibitors (var_ce)
    var_c = ones(Nx, Ny, Nt) * fixpoint(1); % Calcium concentration in the lattice
    var_ce = ones(Nx, Ny, Nt) * fixpoint(2); % Inhibitor concentration in the lattice
    var_akh = zeros(Nx, Ny, Nt); % Empty array for akh (not used here)
    
    % Time step and diffusion coefficients for calcium and inhibitors
    delta_t = 0.1; % Time step (s)
    D_var_c = 20 * tscale; % Diffusion coefficient of calcium
    D_var_ce = 0 * tscale; % Diffusion coefficient of inhibitors (assumed 0 here)
    
    % Random initial values for real_a (parameter used in the system)
    real_a = ones(Nx, Ny) * 100;
    
    % Define spatial axes (not used in lattice model)
    x = Lx / Nx : Lx / Nx : Lx;
    y = Ly / Ny : Ly / Ny : Ly;
    
    % Time loop for simulation
    for ti = 2:1:Nt
        % Update calcium concentration based on diffusion and reactions
        for xi = 2:Nx-1
            for yi = 2:Ny-1
                var_c(xi, yi, ti) = var_c(xi, yi, ti-1) + delta_t * dc(var_c(xi, yi, ti-1), var_ce(xi, yi, ti-1), real_a(xi, yi)) + ...
                    delta_t * D_var_c * ((var_c(xi-1, yi, ti-1) + var_c(xi+1, yi, ti-1) - 2 * var_c(xi, yi, ti-1)) / delta_x^2 + ...
                    (var_c(xi, yi-1, ti-1) + var_c(xi, yi+1, ti-1) - 2 * var_c(xi, yi, ti-1)) / delta_y^2) + ...
                    sigma_c * sqrt(delta_t) * randn; % Add noise
            end
        end

        % Boundary condition handling for calcium at edges of the grid
        % (repeating similar updates for boundary lattice points)

        % Update inhibitor concentration similarly to calcium, without noise
        
        % Periodically plot the calcium concentration distribution
        if mod(ti, 50) == 1
            figure(1); clf;
            imagesc(x, y, var_c(:,:,ti)' / cscale); % Plot the concentration field
            daspect([1 1 1]); % Aspect ratio for the plot
            load Newcolormap.mat; % Load custom colormap
            colormap(NewColormap); colorbar; % Add color bar for scale
            clim([0, 2.5]); % Fix colorbar limits for consistency
            title(['t=', num2str((ti-1) * delta_t), 's']); % Show time on the plot
            pause(0.01); % Small pause to visualize changes
        end
    end
    
    % Calculate half-peak width, period, and other statistics
    % (This part analyzes the time dynamics of calcium concentration at each lattice point)
    
    % Results are stored in the corresponding arrays (mp_repeat, mh_repeat, etc.)
    mp_repeat(krp) = mean_peak;
    mh_repeat(krp) = mean_half;
    mT_repeat(krp) = mean_period;
    vp_repeat(krp) = sqrt(var(c_peak));
    vh_repeat(krp) = sqrt(var(c_half));
    vT_repeat(krp) = sqrt(var(c_period));

    % Save results to an Excel file
    filename = 'origin_rand_inx2i_data.xlsx'; % Output file name
    sheet = 'Sheet1'; % Excel sheet name
    
    % Code to write the calculated data to specific columns in the Excel file
end
