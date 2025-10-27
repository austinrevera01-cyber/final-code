%% Midterm Project Part I.j - 12 V p-code benchmark
% This script acquires steering motor data directly from the run_Indy_car.p
% binary for a 12 V step input.  The resulting time history matches the
% measurement assumptions used in Part I.j of the midterm analysis, enabling
% an apples-to-apples comparison between the analytical gearbox model and
% the provided benchmark.  The captured motor motion is also reflected through
% the 299/14 reduction so that output shaft/tire dynamics can be compared
% directly against the Part I derivations.

clearvars; close all; clc;
clear run_Indy_car; %#ok<CLRUN> Reset persistent state inside the p-code

%% Configuration shared with the Part I.j analysis
Vstep   = 12;          % Input voltage supplied to the p-code [V]
dt      = 0.001;       % Fixed sample time enforced by run_Indy_car.p [s]
t_end   = 0.5;         % Duration of the capture window [s]
CPR     = 500 * 4;     % Encoder counts per motor revolution (quadrature)
MAXCNT  = 4096;        % Encoder rollover count reported by the p-code
gear_N  = 299/14;      % Gear ratio N = omega_motor / omega_output

time_vec = (0:dt:t_end-dt)';
steps    = numel(time_vec);

%% Ensure the p-code is reachable before attempting to execute it
pfile = which('run_Indy_car.p');
if isempty(pfile)
    error('run_Indy_car.p not found on the MATLAB path.');
end
fprintf('Part I.j: Found run_Indy_car.p at %s\n', pfile);

%% One-time initialization call (voltage, velocity, initial state, waypoint)
Vel       = 0;              % Vehicle stationary to isolate motor dynamics
X0_values = [0 0 0 0 0];    % Default initial conditions
WP_FILE   = 0;              % Waypoints disabled for this comparison run
[~, ~, ~] = run_Indy_car(0, Vel, X0_values, WP_FILE);

%% Execute the benchmark loop
acc_counts   = 0;           % Accumulated encoder counts (unwraps rollover)
last_raw     = NaN;         % Previous raw count for rollover detection
theta_counts = zeros(steps, 1);

for k = 1:steps
    [~, ~, counts] = run_Indy_car(Vstep);
    raw = double(counts);

    if isnan(last_raw)
        acc_counts = raw;
    else
        delta = raw - last_raw;
        if delta >  MAXCNT / 2, delta = delta - MAXCNT; end
        if delta < -MAXCNT / 2, delta = delta + MAXCNT; end
        acc_counts = acc_counts + delta;
    end

    last_raw = raw;
    theta_counts(k) = acc_counts;
end

%% Convert counts to motor angle (rad) and differentiate to get speed
theta_m = unwrap(theta_counts * (2 * pi / CPR));
omega_m = gradient(theta_m, dt);

%% Reflect the motion through the gearbox to obtain output variables
theta_out = theta_m / gear_N;    % Output shaft angle [rad]
omega_out = omega_m / gear_N;    % Output shaft speed [rad/s]

%% Report steady-state speed derived from the p-code alone
tail_idx      = max(1, round(0.9 * steps)):steps;
omega_ss      = mean(omega_m(tail_idx));
omega_out_ss  = mean(omega_out(tail_idx));
rpm = @(w) w * 60 / (2 * pi);
fprintf('Part I.j: Steady-state omega_m = %.3f rad/s (%.2f rpm)\n', omega_ss, rpm(omega_ss));
fprintf('          Steady-state omega_{out} = %.3f rad/s (%.2f rpm)\n', ...
        omega_out_ss, rpm(omega_out_ss));

%% Package results for downstream comparison scripts
part1j_data = struct( ...
    'time',               time_vec, ...
    'gear_ratio',         gear_N, ...
    'theta_counts',       theta_counts, ...
    'theta_motor',        theta_m, ...
    'theta_output',       theta_out, ...
    'omega_motor',        omega_m, ...
    'omega_output',       omega_out, ...
    'omega_motor_ss',     omega_ss, ...
    'omega_motor_ss_rpm', rpm(omega_ss), ...
    'omega_output_ss',    omega_out_ss, ...
    'omega_output_ss_rpm', rpm(omega_out_ss));

%% Plot the motor speed trace provided by the p-code
figure('Name', 'Part I.j run\_Indy\_car Benchmark', 'NumberTitle', 'off');
plot(time_vec, omega_out, '', 'LineWidth', 1.5, 'DisplayName', '\omega_{out}');
grid on;
xlabel('Time [s]');
ylabel('Speed [rad/s]');
title('run\_Indy\_car.p Motor vs Output Speed (12 V Step)');
legend('Location', 'best');

%% Close any file handles opened by the p-code implementation
fclose('all');
