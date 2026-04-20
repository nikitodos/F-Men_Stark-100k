function Engine = Firing(Geometry, a, n, C_star, rho, dt, toll_d, m_tot, Gas, lambda, T_nom, rb_nom, plot_flag)
% FIRING  Quasi-steady 0-D transient simulation of a solid rocket motor.
%
%   Engine = Firing(Geometry, a, n, C_star, rho, dt, toll_d, m_tot, ...
%                   Gas, lambda, T_nom, rb_nom, plot_flag)
%
%   Simulates the burning of a cylindrical-perforated grain (or multiple
%   grains) using the Vieille burning rate law r_b = a * p^n.  The chamber
%   pressure is computed from the mass balance p = (a·ρ·c*·Ab/At)^(1/(1-n)).
%   Integration stops when all grains are consumed (web thickness reaches
%   zero or length vanishes).  Performance metrics (Isp, total impulse,
%   delivered c*) are computed by trapezoidal integration.
%
%   -----------------------------------------------------------------------
%   Inputs
%   -----------------------------------------------------------------------
%   Geometry    struct with fields:
%                 .L_g        grain length                       [m]
%                 .N          number of grains                   [-]
%                 .D_perf     initial port diameter              [m]
%                 .D_out      grain outer diameter               [m]
%                 .L_cc       chamber length (including gaps)    [m]
%                 .D_throut   throat diameter (scalar or vector) [m]
%                 .D_exit     nozzle exit diameter (scalar/vec)  [m]
%   a           Vieille law coefficient                          [m/s / Pa^n]
%   n           Vieille law exponent                             [-]
%   C_star      characteristic velocity (ideal, from CEA)        [m/s]
%   rho         propellant density                               [kg/m³]
%   dt          time integration step                            [s]
%   toll_d      fraction of initial web / length for burn‑out    [-]
%   m_tot       total propellant mass                            [kg]
%   Gas         struct with fields:
%                 .gamma   specific heat ratio                   [-]
%                 .R       gas constant                          [J/(kg·K)]
%   lambda      nozzle divergence loss factor (= (1+cos(α))/2)   [-]
%   T_nom       nominal thrust (for error evaluation)            [N]
%   rb_nom      nominal burning rate (for error evaluation)      [m/s]
%   plot_flag   (optional) true to generate figures, false to    [logical]
%               suppress (default = true)
%
%   -----------------------------------------------------------------------
%   Output (struct Engine with length = number of throats)
%   -----------------------------------------------------------------------
%   Engine.t_burn      burn time                                 [s]
%   Engine.Time        time vector                               [s]
%   Engine.P_trace     chamber pressure history                  [Pa]
%   Engine.p_mean      mean chamber pressure                     [bar]
%   Engine.rb_mean     mean burning rate                         [mm/s]
%   Engine.T_mean      mean thrust                               [kN]
%   Engine.T_max       maximum thrust                            [kN]
%   Engine.I_tot       total impulse                             [N·s]
%   Engine.Isp         delivered specific impulse                [s]
%   Engine.c_star      delivered characteristic velocity         [m/s]
%   Engine.eta_cstar   c* efficiency (ideal / delivered)         [-]
%   Engine.C_T         thrust coefficient (constant)             [-]
%
%   -----------------------------------------------------------------------
%   Notes
%   -----------------------------------------------------------------------
%   - The simulation is quasi-steady: chamber pressure responds
%     instantaneously to changes in burning area.
%   - The nozzle is assumed choked throughout the burn.
%   - Erosive burning and two-phase flow effects are neglected.
%
%   See also: SolidUnitDesign, mach_from_area.

if nargin < 13
    plot_flag = true;
end

% Extract geometry
L      = Geometry.L_g;
N      = Geometry.N;
D_perf = Geometry.D_perf;
D_out  = Geometry.D_out;
L_cc   = Geometry.L_cc;

g0 = 9.81;

for ii = 1:length(C_star)

    % Nozzle geometry
    d_t    = Geometry.D_throut(ii);
    c_star = C_star(ii);
    A_t    = pi * (d_t/2)^2;
    A_e    = pi * (Geometry.D_exit(ii)/2)^2;
    eps    = A_e / A_t;

    % Gas properties
    gamma  = Gas(ii).gamma;
    R_gas  = Gas(ii).R;

    % Exit Mach number (isentropic relation)
    obj_M  = @(M) (2/(gamma+1) * (1+(gamma-1)/2*M^2))^((gamma+1)/(2*(gamma-1))) / M - eps;
    Me     = fzero(obj_M, [1.01 50]);

    % Pressure ratio and thrust coefficient
    pr  = (1 + (gamma-1)/2 * Me^2)^(-gamma/(gamma-1));
    C_T = sqrt( 2*gamma^2/(gamma-1) * (2/(gamma+1))^((gamma+1)/(gamma-1)) * ...
        (1 - pr^((gamma-1)/gamma)) ) + pr * eps;

    % Initial grain state
    d_in   = D_perf * ones(1, N);
    l      = L      * ones(1, N);
    active = true(1, N);

    % Initialize histories
    P = 0; Rb = 0; T = 0; Time = 0; t = 0;

    % Burning surface function
    Ab_fun = @(di, li, act) sum( act .* (pi*di.*li + pi/2*(D_out^2 - di.^2)) );
    A_b = Ab_fun(d_in, l, active);

    % Time integration loop
    while any(active)
        t  = t + dt;

        % Chamber pressure and regression rate
        p  = (a * rho * c_star * A_b / A_t)^(1/(1-n));
        rb = a * p^n;

        % Store history
        Time = [Time; t];
        P    = [P;    p];
        Rb   = [Rb;   rb];
        T    = [T;    C_T * A_t * p * lambda];

        % Grain regression
        for g = 1:N
            if ~active(g), continue; end

            l(g)    = l(g)    - 2*dt*rb;
            d_in(g) = d_in(g) + 2*dt*rb;

            % Deactivate grain if burned out
            if (D_out - d_in(g))/D_out < toll_d || l(g) < toll_d
                active(g) = false;
                l(g) = 0;
                d_in(g) = D_out;
            end
        end

        % Update burning area
        A_b = Ab_fun(d_in, l, active);
    end

    % Add final zero point
    Time = [Time; Time(end)+dt];
    P    = [P; 0];
    Rb   = [Rb; 0];
    T    = [T; 0];

    % Performance metrics
    c_star_th = trapz(Time, P) * A_t / m_tot;
    Isp       = trapz(Time, T) / (m_tot * g0);
    I_tot     = trapz(Time, T);

    % Relative errors
    RE_rb = abs(Rb(2:end-1) - rb_nom) / rb_nom * 100;
    RE_T  = abs(T(2:end-1) - T_nom) / T_nom * 100;

    %% -------- Improved plots (only if plot_flag is true) -------- <--- MODIFICA
    if plot_flag
        % Pressure
        figure('Name','Chamber Pressure','NumberTitle','off');
        plot(Time, P*1e-5, 'LineWidth', 1.5);
        xlabel('Time [s]');
        ylabel('Pressure [bar]');
        title('Combustion Chamber Pressure');
        grid on;
        xlim([0, Time(end)])

        % Burning rate
        figure('Name','Burning Rate','NumberTitle','off');
        subplot(2,1,1)
        plot(Time, Rb*1e3, 'LineWidth', 1.5); hold on;
        yline(rb_nom*1e3, '--', 'LineWidth', 1.5);
        xlabel('Time [s]');
        ylabel('r_b [mm/s]');
        title('Burning Rate');
        legend('r_b','r_{b,nom}','Location','best');
        grid on;
        xlim([0, Time(end)])

        subplot(2,1,2)
        plot(Time(2:end-1), RE_rb, 'LineWidth', 1.5);
        xlabel('Time [s]');
        ylabel('Error [%]');
        title('Burning Rate Relative Error');
        grid on;
        xlim([Time(1), Time(end-1)])

        % Thrust
        figure('Name','Thrust','NumberTitle','off');
        subplot(2,1,1)
        plot(Time, T*1e-3, 'LineWidth', 1.5); hold on;
        yline(T_nom*1e-3, '--', 'LineWidth', 1.5);
        xlabel('Time [s]');
        ylabel('Thrust [kN]');
        title('Thrust Profile');
        legend('T','T_{nom}','Location','best');
        grid on;
        xlim([0, Time(end)])

        subplot(2,1,2)
        plot(Time(2:end-1), RE_T, 'LineWidth', 1.5);
        xlabel('Time [s]');
        ylabel('Error [%]');
        title('Thrust Relative Error');
        grid on;
        xlim([Time(1), Time(end-1)])
    end

    %% Results
    Engine(ii).t_burn    = Time(end);
    Engine(ii).Time      = Time;
    Engine(ii).P_trace   = P;
    Engine(ii).p_mean    = mean(P)*1e-5;
    Engine(ii).rb_mean   = mean(Rb*1e3);
    Engine(ii).T_mean    = mean(T)*1e-3;
    Engine(ii).T_max     = max(T)*1e-3;
    Engine(ii).I_tot     = I_tot;
    Engine(ii).Isp       = Isp;
    Engine(ii).c_star    = c_star_th;
    Engine(ii).eta_cstar = c_star / c_star_th;
    Engine(ii).C_T       = C_T;

end
end