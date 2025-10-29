clear; clc; close all;

part1_results = run_part1_motor_model();
run_part1j_benchmark(part1_results);
run_part2_pd_and_vehicle();

function part1 = run_part1_motor_model()
    %% Part I (sim_prod.m) - Motor modeling and gearbox study

    params.R      = 0.611;      % Armature resistance [Ohm]
    params.L      = 0.000119;   % Armature inductance [H]
    params.Kb     = 0.025879;   % Back EMF constant [V/(rad/s)]
    params.Ki     = 0.0259;     % Torque constant [N*m/A]
    params.Jm     = 3.35e-6;    % Rotor inertia [kg*m^2]
    params.bm     = 4.63e-6;    % Rotor viscous friction [N*m*s/rad]
    params.Tstall = 1.02;       % Stall torque [N*m]
    params.w_nl   = 922.581;    % No-load speed [rad/s]
    params.tau_e  = 0.0028327;  % Electrical time constant [s]
    params.Vs     = 12;         % Nominal comparison voltage [V]

    gear.N        = 299/14;     % Motor-to-output ratio
    gear.Jload    = 8.0e-8;     % Load inertia at output [kg*m^2]
    gear.encoder_counts = 2000; % Counts per revolution (quadrature)

    disp('--- Part I: Motor model with and without inductance ---');

    % Transfer function with inductance retained
    num_with_L = params.Ki;
    den_with_L = [params.Jm * params.L, ...
                  params.bm * params.L + params.Jm * params.R, ...
                  params.bm * params.R + params.Ki * params.Kb];
    G_omega_with_L = tf(num_with_L, den_with_L, 'Variable', 's');

    % Transfer function neglecting inductance
    num_no_L = params.Ki;
    den_no_L = [params.Jm * params.R, params.bm * params.R + params.Ki * params.Kb];
    G_omega_no_L = tf(num_no_L, den_no_L, 'Variable', 's');

    J_eff = params.Jm;
    b_eff = params.bm + (params.Ki * params.Kb) / params.R;
    fprintf('J_eff = %.4e kg*m^2, b_eff = %.4e N*m*s/rad\n', J_eff, b_eff);

    eig_with_L = pole(G_omega_with_L);
    eig_no_L   = pole(G_omega_no_L);

    tau_with_L = -1 ./ real(eig_with_L);
    tau_no_L   = -1 ./ real(eig_no_L);

    fprintf('Eigenvalues with inductance: %s\n', mat2str(eig_with_L, 4));
    fprintf('Time constants with inductance: %s s\n', mat2str(tau_with_L, 4));
    fprintf('Eigenvalues without inductance: %s\n', mat2str(eig_no_L, 4));
    fprintf('Time constants without inductance: %s s\n', mat2str(tau_no_L, 4));

    % Step responses for the two plant models
    time_span = 0:1e-4:0.1;
    [y_with_L, t_with_L] = step(params.Vs * G_omega_with_L, time_span);
    [y_no_L,   t_no_L]   = step(params.Vs * G_omega_no_L, time_span);

    figure('Name','Part I.f Step Response','NumberTitle','off');
    plot(t_with_L, y_with_L, 'LineWidth', 1.5); hold on;
    plot(t_no_L, y_no_L, '--', 'LineWidth', 1.5);
    grid on;
    title('Motor Speed Response to 12 V Step');
    xlabel('Time [s]');
    ylabel('\omega_m [rad/s]');
    legend('With inductance','Without inductance','Location','southeast');
    xlim([0,0.02]);

    stepinfo_with_L = stepinfo(params.Vs * G_omega_with_L);
    stepinfo_no_L   = stepinfo(params.Vs * G_omega_no_L);
    fprintf('With inductance: steady-state %.2f rad/s, settling %.4f s\n', ...
        dcgain(params.Vs * G_omega_with_L), stepinfo_with_L.SettlingTime);
    fprintf('Without inductance: steady-state %.2f rad/s, settling %.4f s\n', ...
        dcgain(params.Vs * G_omega_no_L), stepinfo_no_L.SettlingTime);

    tau_electrical = params.L / params.R;
    tau_mech_model = J_eff / b_eff;
    ratio_tau = tau_electrical / tau_mech_model;

    % Gearbox extension
    % The gearbox specification provides the load inertia at the OUTPUT
    % shaft.  To express the plant dynamics in terms of the motor-side speed
    % (which is what the electrical equations couple to), the reflected
    % inertia must be divided by N^2.  The analytical model previously
    % multiplied by N^2, which effectively treated the 8e-8 kg*m^2 load as if
    % it were mounted directly on the motor shaft scaled up by the gear
    % ratio.  That inflated the effective inertia by ~four orders of
    % magnitude and produced a much slower step response than the physical
    % benchmark.  Using the correct reflected inertia keeps the analytical
    % model consistent with the measured run_Indy_car behaviour.
    J_total_motor = params.Jm + gear.Jload / gear.N^2;

    den_with_L_gear = [J_total_motor * params.L, ...
                       params.bm * params.L + J_total_motor * params.R, ...
                       params.bm * params.R + params.Ki * params.Kb];
    G_omega_with_L_gear = tf(params.Ki, den_with_L_gear, 'Variable', 's');
    G_omega_out_with_L = (1/gear.N) * G_omega_with_L_gear;

    fprintf('Effective inertia with gearbox = %.4e kg*m^2\n', J_total_motor);
    time_span_gear = 0:1e-3:2;
    [omega_out_step, t_gear] = step(params.Vs * G_omega_out_with_L, time_span_gear);

    figure('Name','Part I.i Gearbox Motor Model','NumberTitle','off');
    plot(t_gear, omega_out_step, 'LineWidth', 1.5);
    grid on;
    title('Motor Model with Gearbox (12 V Step)');
    xlabel('Time [s]');
    ylabel('\omega_{out} [rad/s]');
    xlim([0, 0.5]);

    part1 = struct();
    part1.params                 = params;
    part1.gear                   = gear;
    part1.G_omega_with_L         = G_omega_with_L;
    part1.G_omega_no_L           = G_omega_no_L;
    part1.G_omega_with_L_gear    = G_omega_with_L_gear;
    part1.G_omega_out_with_L     = G_omega_out_with_L;
    part1.stepinfo_with_L        = stepinfo_with_L;
    part1.stepinfo_no_L          = stepinfo_no_L;
    part1.tau_mech_model         = tau_mech_model;
    part1.tau_electrical         = tau_electrical;
    part1.tau_ratio              = ratio_tau;
    part1.gear_step_time         = t_gear(:);
    part1.gear_step_response     = omega_out_step(:);
end

function part1j = run_part1j_benchmark(part1_model)
    %% Part I.j (pcode_test.m) - Pull benchmark data from run_Indy_car.p
    disp('--- Part I.j: Capturing run_Indy_car benchmark ---');

    if nargin < 1
        part1_model = struct();
    end

    Vstep   = 12;          % Input voltage [V]
    dt      = 0.001;       % Sample time [s]
    t_end   = 0.5;         % Capture window [s]
    CPR     = 500 * 4;     % Encoder counts per revolution (quadrature)
    MAXCNT  = 4096;        % Encoder rollover count
    gear_N  = 299/14;      % Gear ratio

    time_vec = (0:dt:t_end-dt)';
    steps    = numel(time_vec);

    pfile = which('run_Indy_car.p');
    if isempty(pfile)
        error('run_Indy_car.p not found on the MATLAB path.');
    end
    fprintf('Found run_Indy_car.p at %s\n', pfile);

    clear run_Indy_car;

    Vel       = 0;
    X0_values = [0 0 0 0 0];
    WP_FILE   = 0;
    run_Indy_car(0, Vel, X0_values, WP_FILE);

    acc_counts   = 0;
    last_raw     = NaN;
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

    counts2rad = 2 * pi / CPR;

    % The encoder in run_Indy_car is mounted after the gearbox, so the
    % accumulated counts already reflect the output shaft position.  Convert
    % that measurement to radians first, then infer the motor-side quantities
    % by multiplying with the gear ratio instead of dividing.  This keeps the
    % benchmark output trajectory consistent with the Part I.i analytical
    % model, which predicts the output speed.
    theta_out = unwrap(theta_counts * counts2rad);
    omega_out = gradient(theta_out, dt);

    theta_m = theta_out * gear_N;
    omega_m = omega_out * gear_N;

    tail_idx      = max(1, round(0.9 * steps)):steps;
    omega_ss      = mean(omega_m(tail_idx));
    omega_out_ss  = mean(omega_out(tail_idx));
    rpm = @(w) w * 60 / (2 * pi);

    fprintf('Steady-state omega_m = %.3f rad/s (%.2f rpm)\n', omega_ss, rpm(omega_ss));
    fprintf('Steady-state omega_out = %.3f rad/s (%.2f rpm)\n', omega_out_ss, rpm(omega_out_ss));

    figure('Name', 'Part I.j run_Indy_car Benchmark', 'NumberTitle', 'off');
    plot(time_vec, omega_out, 'LineWidth', 1.5, 'DisplayName', '\omega_{out}');
    grid on;
    xlabel('Time [s]');
    ylabel('Speed [rad/s]');
    title('run\_Indy\_car.p Motor vs Output Speed (12 V Step)');
    legend('Location', 'best');

    if isfield(part1_model, 'gear_step_time') && isfield(part1_model, 'gear_step_response')
        figure('Name', 'Part I.i vs Part I.j Output Speed', 'NumberTitle', 'off');
        plot(part1_model.gear_step_time, part1_model.gear_step_response, 'LineWidth', 1.5, ...
            'DisplayName', 'Part I.i Model');
        hold on;
        plot(time_vec, omega_out, '--', 'LineWidth', 1.5, 'DisplayName', 'Part I.j Benchmark');
        grid on;
        xlabel('Time [s]');
        ylabel('Speed [rad/s]');
        title('Part I.i Model vs. Part I.j Benchmark Output Speed');
        legend('Location', 'best');
        xlim([0, min(max(part1_model.gear_step_time), max(time_vec))]);
    end

    fclose('all');

    part1j = struct( ...
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
end

function part2 = run_part2_pd_and_vehicle()
    %% Part II (simprodpt2.m) - PD design, bench test, and vehicle sim
    disp('--- Part II: PD design and vehicle simulation ---');

    % Motor + gear parameters
    R  = 0.611;
    L  = 0.000119;
    Ki = 0.0259;
    Kb = 0.025879;
    Jm = 3.35e-6;

    N = 299/14;
    Jg_ref = (8.0e-8) * N^2;
    Jload = Jm + Jg_ref;
    J_out = Jload;

    CPR    = 500 * 4;
    MAXCNT = 4096;

    % Transfer function model
    s = tf('s');
    P_motor = Ki / ((Jload * s) * (L * s + R) + Ki * Kb);
    P_theta = P_motor / s;

    % Design specs
    zeta = 0.69;
    Ts   = 0.025;
    wn   = 266.6666666667;
    fprintf('Target wn = %.1f rad/s\n', wn);

    wbw = wn * sqrt(1 - 2 * zeta^2 + sqrt(2 + 4 * zeta^4));
    fprintf('Closed-loop bandwidth ≈ %.1f rad/s (%.1f Hz)\n', wbw, wbw / (2 * pi));

    % Controller design
    Kp = (wn^2 * J_out) / Ki;
    Kd = (2 * zeta * wn * J_out - Kb * Ki / R) / Ki;
    fprintf('Kp = %.2f, Kd = %.5f\n', Kp, Kd);

    C = Kp + Kd * s;
    T_closed = feedback(C * P_theta, 1);
    info = stepinfo(T_closed);
    fprintf('Closed-loop settling time ≈ %.4f s\n', info.SettlingTime);

    % Sinusoidal tracking test
    f_ref = 10;
    t = 0:1e-5:0.2;
    ref = 0.1 * sin(2 * pi * f_ref * t);
    y = lsim(T_closed, ref, t);

    figure('Name', 'Part II.b Responses', 'NumberTitle', 'off');
    subplot(1, 2, 1);
    step(T_closed);
    title('Step Response (Constant Ref)');
    ylabel('Position [rad]');
    grid on;

    subplot(1, 2, 2);
    plot(t, y, 'b', 'LineWidth', 1.3); hold on;
    plot(t, ref, 'r--', 'LineWidth', 1.2);
    title(sprintf('Sinusoidal Tracking (%.1f Hz)', f_ref));
    xlabel('Time [s]');
    ylabel('Position [rad]');
    legend('Output', 'Reference');
    grid on;

    figure('Name', 'Part II.b Closed-loop Bode', 'NumberTitle', 'off');
    bode(T_closed);
    grid on;
    title('Bode Plot of Closed-loop Position System');

    clear run_Indy_car;

    dt    = 0.001;
    t_end = 10;
    steps = round(t_end / dt);
    t_vec = (0:steps - 1).' * dt;

    Vmax       = 24;
    counts2rad = 2 * pi / CPR;

    theta_ref_deg = 30;
    theta_ref     = deg2rad(theta_ref_deg);

    Vel       = 0;
    X0_values = [0 0 0 0 0];
    WP_FILE   = 0;
    run_Indy_car(0, Vel, X0_values, WP_FILE);

    acc_counts     = 0;
    last_raw       = NaN;
    motor_counts   = zeros(steps, 1);
    theta_motor    = zeros(steps, 1);
    theta_output   = zeros(steps, 1);
    control_volts  = zeros(steps, 1);
    kp_component   = zeros(steps, 1);
    kd_component   = zeros(steps, 1);

    err_prev       = theta_ref;
    theta_prev_out = 0;

    fprintf('Running PD bench test with Kp = %.2f, Kd = %.5f\n', Kp, Kd);

    for k = 1:steps
        error_k  = theta_ref - theta_prev_out;
        deriv_k  = (error_k - err_prev) / dt;
        u_p      = Kp * error_k;
        u_d      = Kd * deriv_k;
        u_unsat  = u_p + u_d;
        u        = max(min(u_unsat, Vmax), -Vmax);

        control_volts(k) = u;
        kp_component(k)  = u_p;
        kd_component(k)  = u_d;

        [~, ~, counts] = run_Indy_car(u);
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

        motor_counts(k) = acc_counts;
        theta_motor(k)  = acc_counts * counts2rad;
        theta_output(k) = theta_motor(k) / N;

        err_prev       = error_k;
        theta_prev_out = theta_output(k);
    end

    theta_output_deg = rad2deg(theta_output);
    tol              = 0.02 * abs(theta_ref);
    settle_idx       = find(abs(theta_output - theta_ref) > tol, 1, 'last');
    if isempty(settle_idx)
        Ts_measured = 0;
    else
        Ts_measured = t_vec(min(settle_idx + 1, steps));
    end

    figure('Name', 'Part II.c Steering Response', 'NumberTitle', 'off');
    plot(t_vec, theta_output_deg, 'b', 'LineWidth', 1.5); hold on;
    plot(t_vec, theta_ref_deg * ones(size(t_vec)), 'r--', 'LineWidth', 1.3);
    grid on;
    legend('run\_Indy\_car output', 'Reference', 'Location', 'southeast');
    xlabel('Time [s]');
    ylabel('Steering Angle [deg]');
    title('Part II.c:run\_Indy\_car.p with controller');
    xlim([0, 1]);

    fclose('all');
    clear run_Indy_car;

    % Vehicle state simulation
    clear run_Indy_car; 

    V_cmd       = 12;
    vehicle_vel = 10;
    dt_vehicle  = 0.001;
    t_final     = 5;

    encoder_cpr = CPR;
    steer_ratio = 15;
    delta_max   = deg2rad(20);

    num_steps_vehicle = round(t_final / dt_vehicle);
    time_vec_vehicle  = (0:num_steps_vehicle - 1)' * dt_vehicle;

    X0      = [0 0 0 0 0];
    WP_FILE = 0;
    run_Indy_car(0, vehicle_vel, X0, WP_FILE);

    motor_counts_vehicle = zeros(num_steps_vehicle, 1);
    yaw_rate             = zeros(num_steps_vehicle, 1);
    heading              = zeros(num_steps_vehicle, 1);

    for k = 1:num_steps_vehicle
        [gps, yaw_k, counts_k] = run_Indy_car(V_cmd);
        motor_counts_vehicle(k) = double(counts_k);
        yaw_rate(k)             = yaw_k;
        heading(k)              = gps(3);
    end

    motor_angle_vehicle = unwrap(motor_counts_vehicle * (2 * pi / encoder_cpr));
    delta_tire = motor_angle_vehicle / steer_ratio;
    delta_tire = min(max(delta_tire, -delta_max), delta_max);

    ma_window = ones(5, 1) / 5;
    yaw_rate  = filter(ma_window, 1, yaw_rate);

    figure('Name', 'Part II.d Vehicle Output States', 'NumberTitle', 'off');
    subplot(3, 1, 1);
    plot(time_vec_vehicle, rad2deg(delta_tire), 'LineWidth', 1.5);
    ylabel('\delta_{tire} [deg]');
    grid on;
    title('Steering Response for 12 V Command at 10 m/s');
    xlim([0, 0.07]);

    subplot(3, 1, 2);
    plot(time_vec_vehicle, yaw_rate, 'LineWidth', 1.5);
    ylabel('Yaw rate [rad/s]');
    grid on;
    xlim([0, 5]);

    subplot(3, 1, 3);
    plot(time_vec_vehicle, rad2deg(heading), 'LineWidth', 1.5);
    xlabel('Time [s]');
    ylabel('\psi [deg]');
    grid on;
    xlim([0, 5]);

    steady_idx_vehicle = round(0.9 * num_steps_vehicle):num_steps_vehicle;
    delta_ss           = mean(delta_tire(steady_idx_vehicle));
    yaw_ss             = mean(yaw_rate(steady_idx_vehicle));
    psi_ss             = mean(heading(steady_idx_vehicle));

    fclose('all');
    clear run_Indy_car;

    part2 = struct();
    part2.Kp               = Kp;
    part2.Kd               = Kd;
    part2.closed_loop_tf   = T_closed;
    part2.stepinfo         = info;
    part2.t_tracking       = t(:);
    part2.ref_tracking     = ref(:);
    part2.output_tracking  = y(:);
    part2.t_pd             = t_vec;
    part2.theta_output_deg = theta_output_deg;
    part2.control_volts    = control_volts;
    part2.kp_component     = kp_component;
    part2.kd_component     = kd_component;
    part2.Ts_measured      = Ts_measured;
    part2.vehicle_time     = time_vec_vehicle;
    part2.vehicle_delta    = delta_tire;
    part2.vehicle_yaw_rate = yaw_rate;
    part2.vehicle_heading  = heading;
    part2.delta_ss         = delta_ss;
    part2.yaw_ss           = yaw_ss;
    part2.psi_ss           = psi_ss;
end