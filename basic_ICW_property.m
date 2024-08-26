%% IPR model amplification
clear; % Clear workspace variables
close; % Close all figure windows
clc;   % Clear command window

% Initialize variables
indu = 0:2:200; % AKH concentrations range from 0 to 200 with a step of 2
Nrp = length(indu); % Number of different AKH concentrations
max_record = zeros(1, Nrp); % Record for maximum calcium levels
min_record = zeros(1, Nrp); % Record for minimum calcium levels
max1st_record = zeros(1, Nrp); % Record for first peak calcium levels
p_record = zeros(1, Nrp); % Record for oscillation period
exci_record = zeros(1, Nrp); % Record for excitation time duration

% Loop through each concentration
for ki = 1:Nrp

    % Scaling and model parameters
    tscale = 1.5;
    cscale = 1000;
    akhscale = 0.0005;
    akhbase = 240;

    k1 = 400;          % μM^(-1) s^(-1)
    k1_minus = 52;     % s^(-1)
    k5 = 20;           % μM^(-1) s^(-1)
    k5_minus = 1.64;   % s^(-1)

    gamma = 2;
    delta = 1;

    Vs = 0.9 * tscale * cscale;  % μM s^(-1)
    Ks = 0.1 * cscale;           % μM

    Vp = 0.05 * tscale * cscale; % μM s^(-1)
    Kp = 0.3 * cscale;           % μM

    alpha1 = 0.07 * tscale * cscale * akhscale; % s^(-1)
    alpha0 = 0.01 * tscale * cscale + akhbase * alpha1; % μM s^(-1)

    kf = 1.11 * tscale;     % s^(-1)
    kleak = 0.02 * tscale;  % s^(-1)

    K1 = k1_minus / k1 / akhscale + akhbase;
    K5 = k5_minus / k5 * cscale;

    p = indu(ki); % Current AKH concentration
    y = 0.3; % A constant related to the probability of channel opening

    % Define the fluxes (J_irp: IP3 receptor flux, J_leak: Leak flux, J_serca: SERCA pump flux, J_mem: Membrane flux)
    J_irp = @(c, ce) kf * ((p + akhbase) .* c .* (1 - y) ./ (p + K1) ./ (c + K5)).^3 .* (ce - c);
    J_leak = @(c, ce) kleak * (ce - c);
    J_serca = @(c, ce) Vs * c.^2 ./ (Ks^2 + c.^2);
    J_mem = @(c, ce) delta * (alpha0 + alpha1 * p - Vp * c.^2 ./ (Kp^2 + c.^2));

    % Define the differential equations for calcium concentrations in the cytosol (dc) and ER (dce)
    dc = @(c, ce) J_irp(c, ce) + J_leak(c, ce) - J_serca(c, ce) + J_mem(c, ce);
    dce = @(c, ce) -gamma * (J_irp(c, ce) + J_leak(c, ce) - J_serca(c, ce));

    % Define the ODE system
    dxdt = @(t, x) [dc(x(1), x(2)); dce(x(1), x(2))];

    % Set the time span for simulation
    tspan = [0 2000];

    % Solve the ODE system with initial conditions (fixpoint)
    fixpoint = [0.29, 20.6] * cscale;
    [t, x] = ode23s(dxdt, tspan, fixpoint);

    % Record maximum and minimum calcium levels in the last third of the simulation
    max_record(ki) = max(x(round(length(x(:,1)) / 3 * 2):end, 1));
    mini = min(x(round(length(x(:,1)) / 3 * 2):end, 1));
    min_record(ki) = mini;

    % Record the first peak calcium level
    max1st_record(ki) = max(x(:, 1));

    % Calculate the oscillation period and excitation time
    period = 0;
    try
        peak_site = [];
        % Identify peaks in the calcium concentration time series
        for i = 2:length(t) - 1
            if x(i, 1) > x(i-1, 1) && x(i, 1) > x(i+1, 1)
                peak_site = [peak_site, i];
            end
        end
        
        N_p = length(peak_site) - 2; % Number of peaks after the initial transient
        period = (t(peak_site(end)) - t(peak_site(2))) / N_p; % Average oscillation period
        
        pe_start = zeros(N_p, 1);
        pe_end = zeros(N_p, 1);
        % Calculate the start and end times of each oscillation based on a threshold
        for i = 1:N_p
            pe_start(i) = find(x(1:peak_site(i+1), 1) < 2 * mini, 1, 'last');
            pe_end(i) = find(x(peak_site(i+1):end, 1) < 2 * mini, 1, 'first') + peak_site(i+1) - 1;
        end
        exci_record(ki) = mean(t(pe_end) - t(pe_start)); % Mean excitation duration
    end

    p_record(ki) = period; % Record oscillation period
end

% Calculate the ratio of maximum to minimum calcium levels (oscillation amplitude)
radio_record = max_record ./ min_record;

% Identify the range of AKH concentrations where oscillations occur
t_osci_start = find(radio_record > 2, 1);
t_osci_end = find(radio_record > 2, 1, 'last');

% Plot results
figure;
hold on;
plot(indu, min_record, 'g', 'LineWidth', 2);
plot(indu, max_record, 'r', 'LineWidth', 2);
plot(indu, max1st_record, 'b', 'LineWidth', 2);
legend('Resting state calcium level during oscillation', 'Excited state calcium level during oscillation', 'First peak calcium level during oscillation');
title('Amplitude changes with AKH concentration');
hold off;

figure;
plot(indu, radio_record, 'r', 'LineWidth', 2);
title('Amplification factor changes with AKH concentration');

figure;
hold on;
plot(indu(t_osci_start:t_osci_end), p_record(t_osci_start:t_osci_end), 'r', 'LineWidth', 2);
plot(indu(t_osci_start:t_osci_end), exci_record(t_osci_start:t_osci_end), 'g', 'LineWidth', 2);
plot(indu(t_osci_start:t_osci_end), p_record(t_osci_start:t_osci_end) - exci_record(t_osci_start:t_osci_end), 'b', 'LineWidth', 2);
legend('Oscillation period', 'Time in oscillatory state', 'Time in resting state');
ylim([0, 150]);
title('Oscillation period changes with AKH concentration');
hold off;
