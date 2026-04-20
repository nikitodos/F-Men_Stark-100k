function MotorEngine = SolidUnitDesign(Design, Case, Case_Metal, TP_Case, th_TPC, N_grains, T_amb, CF_rb)
% SOLIDUNITDESIGN  Parametric design loop for solid rocket motor geometry.
%
%   MotorEngine = SolidUnitDesign(Design, Case, Case_Metal, TP_Case, ...
%                                 th_TPC, N_grains, safehook, T_amb, CF_rb)
%
%   For each propellant formulation defined in the Case array, the function:
%     - Computes ideal nozzle performance (c*, C_T, I_sp).
%     - Determines the required propellant mass and volume.
%     - Solves the grain geometry (port radius, chamber radius, length)
%       using solve_geometry_fzero2, assuming a cylindrical BATES grain
%       burning from both ends and internal surface.
%     - Performs a 1‑D transient thermal analysis of the combustion chamber
%       wall (insulation + case) to verify that temperature limits are met.
%     - Stores all viable configurations in the output struct array.
%
%   -----------------------------------------------------------------------
%   Inputs
%   -----------------------------------------------------------------------
%   Design      struct with mission requirements:
%                 .Thrust   nominal thrust                         [N]
%                 .P_cc     chamber pressure                       [Pa]
%                 .I_tot    total impulse                          [N·s]
%                 .alpha    nozzle divergence half‑angle           [deg]
%                 .beta     nozzle convergence half‑angle          [deg]
%                 .P_e      nozzle exit pressure (ideal expanded)  [Pa]
%                 .a_bar    Vieille coefficient (mm/s/bar^n)       [-]
%                 .n_bar    Vieille exponent                       [-]
%                 .a_SI     Vieille coefficient (m/s/Pa^n)         [-]
%                 .n_SI     Vieille exponent (same as n_bar)       [-]
%   Case        array of structs (one per propellant) with fields:
%                 .Gas      gas properties (gamma, R, rho, mu, k, Pr)
%                 .T_cc     chamber temperature                    [K]
%                 .rho      propellant density                     [kg/m³]
%                 .k        propellant thermal conductivity        [W/(m·K)]
%                 .cp       propellant specific heat               [J/(kg·K)]
%                 .AP       AP mass fraction (for printout)        [%]
%                 .HTPB     HTPB mass fraction                     [%]
%   Case_Metal  struct with case material properties:
%                 .k        thermal conductivity                   [W/(m·K)]
%                 .rho      density                                [kg/m³]
%                 .cp       specific heat                          [J/(kg·K)]
%                 .yield    yield strength                         [Pa]
%                 .T_melt   melting temperature                    [K]
%   TP_Case     struct with thermal protection properties (same fields)
%   th_TPC      insulation thickness                               [m]
%   N_grains    number of grains                                   [-]
%   safehook    structural safety factor (yield / allowable)       [-]
%   T_amb       ambient temperature                                [K]
%   CF_rb       correction factor for regression burning rate      [-]
%
%   -----------------------------------------------------------------------
%   Output
%   -----------------------------------------------------------------------
%   MotorEngine   array of structs, each containing:
%                 - Case, AP, HTPB, Tcc, rb, tb, I_sp, c_star, c_T, ...
%                 - Geom (all geometric parameters)
%                 - Thermal analysis results (printed, not stored)
%
%   See also: solve_geometry_fzero2, Firing.

% Number of cases
NN = length(Case);

% Design inputs
Thrust_des = Design.Thrust;
P_cc = Design.P_cc;
I_tot = Design.I_tot;
alpha = Design.alpha;
beta = Design.beta;
P_e = Design.P_e;
a_bar = Design.a_bar;
n_bar = Design.n_bar;
a_SI = Design.a_SI;
n_SI = Design.n_SI;

lambda = (1 - cosd(alpha))/2;
g0 = 9.81;

jj = 1;

for ii = 1:NN

    % Gas properties
    gamma_g = Case(ii).Gas.gamma;
    R_g     = Case(ii).Gas.R;
    rho_g   = Case(ii).Gas.rho;
    mu_g    = Case(ii).Gas.mu;
    k_g     = Case(ii).Gas.k;
    Pr_g    = Case(ii).Gas.Pr;
    T_cc_i  = Case(ii).T_cc;

    Gas = Case(ii).Gas;

    % Case material
    k_case  = Case_Metal.k;
    rho_case = Case_Metal.rho;
    cp_case = Case_Metal.cp;

    % Grain properties
    k_grain = Case(ii).k;
    rho_grain = Case(ii).rho;
    cp_grain = Case(ii).cp;

    % Insulation properties
    k_iso   = TP_Case.k;
    rho_iso = TP_Case.rho;
    cp_iso  = TP_Case.cp;

    % Nozzle and performance
    c_star = sqrt(R_g*T_cc_i/gamma_g) * ((gamma_g+1)/2)^((gamma_g+1)/(2*(gamma_g-1)));
    c_T    = sqrt(2*gamma_g^2/(gamma_g-1) * (2/(gamma_g+1))^((gamma_g+1)/(gamma_g-1)) * ...
        (1-(P_e/P_cc)^((gamma_g-1)/gamma_g)));

    u_e    = sqrt(2*gamma_g/(gamma_g-1)*R_g*T_cc_i*(1-(P_e/P_cc)^((gamma_g-1)/gamma_g)));
    T_e    = T_cc_i*(P_e/P_cc)^((gamma_g-1)/gamma_g);
    rho_e  = P_e/(R_g*T_e);
    Me     = u_e/sqrt(gamma_g*R_g*T_e);

    I_sp   = c_T*c_star/g0;
    m_sp   = I_tot/(I_sp*g0);
    V_sp   = m_sp/rho_grain;
    t_b    = I_tot / Thrust_des;

    mdot = I_tot / (I_sp * g0 * t_b);

    A_e = mdot/(rho_e*u_e);
    eps = 1/(((gamma_g+1)/2)^(1/(gamma_g-1))*(P_e/P_cc)^(1/gamma_g)* ...
        sqrt((gamma_g+1)/(gamma_g-1)*(1-(P_e/P_cc)^((gamma_g-1)/gamma_g))));

    A_t = A_e/eps;
    r_e = sqrt(A_e/pi);
    r_t = sqrt(A_t/pi);
    L_div = (r_e-r_t)/tand(alpha);

    % Burning rate
    rb = (a_SI * P_cc ^ n_SI)*CF_rb;

    % Solve grain geometry
    Solution = solve_geometry_fzero2(V_sp, rb, t_b, N_grains);

    if Solution.sol ~= 0

        r_p  = Solution.r_p;
        r_cc = Solution.r_cc;
        L_g  = Solution.L_g;
        web  = Solution.web;

        % Chamber geometry
        L_conv  = (r_cc-r_t)/tand(beta);
        L_cc    = L_g;

        % Flow in chamber
        u_cc = mdot / ( pi * r_cc^2 * rho_g );
        M_cc = u_cc / sqrt( gamma_g * R_g * T_cc_i);

        % Case thickness
        th_case = 2 * P_cc * ( r_cc + th_TPC ) / (Case_Metal.yield );

        % Burning area
        Ab = N_grains * ( 2 * pi * r_p * L_g + 2 * pi * ( r_cc^2 - r_p^2 ) );

        % Mass flow consistency check
        mdot_check = rb * rho_grain * Ab;
        mdot_rd = abs( mdot_check - mdot ) / mdot * 100;

        if mdot_rd >= 5
            warning('Mass flow mismatch exceeds 5%%');
        end

        % Geometric ratios
        V_cc    = pi*r_cc^2*L_cc;
        LD      = L_cc/(2*r_cc);
        b_f     = web/r_cc;
        V_f     = V_sp/V_cc;
        eps_cc  = pi*r_cc^2/A_t;

        %% -------- THERMAL ANALYSIS (WITH & WITHOUT TPS) --------

        % Setup parameters
        Nz = 1e3;
        z = linspace(0, L_cc, Nz);
        dz = L_cc / Nz;
        dt = 1e-2;
        time_vec = 0:dt:t_b;
        Nt = length(time_vec);

        % --- GEOMETRIC DEFINITIONS ---

        % 1. CASE WITH TPS
        r_inner_tps    = r_cc;
        r_inter_tps    = r_cc + th_TPC;        % Interface TBC-Metal
        r_outer_tps    = r_cc + th_TPC + th_case; % Outer metal radius

        % 2. CASE NO TPS (Direct Contact)
        r_inner_no     = r_cc;
        r_outer_no     = r_cc + th_case;       % Outer metal radius (no TBC)

        % --- CALCULATE SEGMENT MASSES ---
        % With TPS: Metal layer is between r_inter_tps and r_outer_tps
        M_iso_seg      = rho_iso  * (pi * (r_inter_tps^2 - r_inner_tps^2) * dz);
        M_case_seg     = rho_case * (pi * (r_outer_tps^2 - r_inter_tps^2) * dz);

        % No TPS: Metal layer is between r_inner_no and r_outer_no
        M_case_seg_no_tps = rho_case * (pi * (r_outer_no^2 - r_inner_no^2) * dz);

        % Initialize temperature arrays
        T_grain_surf = ones(Nt, Nz) * T_amb;
        T_iso_int    = ones(Nt, Nz) * T_amb;
        T_case_int   = ones(Nt, Nz) * T_amb;
        T_case_ext   = ones(Nt, Nz) * T_amb;

        T_case_int_no_tps = ones(Nt, Nz) * T_amb;
        T_case_ext_no_tps = ones(Nt, Nz) * T_amb;

        % Main time loop
        for t_idx = 1:Nt-1
            t_curr = time_vec(t_idx);
            front_head = rb * t_curr;
            front_tail = L_cc - rb * t_curr;

            for i = 1:Nz
                r_grain_l = r_p + rb * t_curr;

                % --- 1. CASE: WITH TPS ---
                is_exposed_tps = (z(i) < front_head) || (z(i) > front_tail) || (r_grain_l >= r_inner_tps);

                if is_exposed_tps
                    r_flow_tps = r_inner_tps;
                    R_grain_local_tps = 0;
                else
                    r_flow_tps = r_grain_l;
                    R_grain_local_tps = log(r_inner_tps / r_grain_l) / (2 * pi * k_grain * dz);
                end

                % Convection heat transfer (Bartz/Standard)
                u_local_tps = mdot / (pi * r_flow_tps^2 * rho_g);
                Re_local_tps = (rho_g * 2 * r_flow_tps * u_local_tps) / mu_g;
                f_tps = (1.82 * log10(Re_local_tps) - 1.64)^ (-2);
                Nu_tps = (f_tps/8 * Re_local_tps * Pr_g)/ (1.07 + 12.7 * sqrt(f_tps/8) * (Pr_g^(2/3) - 1));
                h_int_tps = Nu_tps * k_g / (2 * r_flow_tps);

                % Resistances
                R_conv_i_tps = 1 / (h_int_tps * 2 * pi * r_flow_tps * dz);
                R_iso        = log(r_inter_tps / r_inner_tps) / (2 * pi * k_iso * dz);
                R_case_tps   = log(r_outer_tps / r_inter_tps) / (2 * pi * k_case * dz);
                R_conv_e_tps = 1 / (10 * 2 * pi * r_outer_tps * dz);

                % Heat Fluxes
                Q_in_to_iso   = (T_cc_i - T_iso_int(t_idx, i)) / (R_conv_i_tps + R_grain_local_tps);
                Q_iso_to_case = (T_iso_int(t_idx, i) - T_case_int(t_idx, i)) / R_iso;
                Q_case_to_ext = (T_case_int(t_idx, i) - T_case_ext(t_idx, i)) / R_case_tps;
                Q_env         = (T_case_ext(t_idx, i) - T_amb) / R_conv_e_tps;

                % Update Temperatures
                T_grain_surf(t_idx, i) = T_cc_i - Q_in_to_iso * R_conv_i_tps;
                T_iso_int(t_idx+1, i) = T_iso_int(t_idx, i) + ((Q_in_to_iso - Q_iso_to_case) / (0.5 * M_iso_seg * cp_iso)) * dt;
                T_case_int(t_idx+1, i) = T_case_int(t_idx, i) + ((Q_iso_to_case - Q_case_to_ext) / (0.5 * M_case_seg * cp_case)) * dt;
                T_case_ext(t_idx+1, i) = T_case_ext(t_idx, i) + ((Q_case_to_ext - Q_env) / (0.5 * M_case_seg * cp_case)) * dt;

                % --- 2. CASE: NO TPS (Direct Contact) ---
                is_exposed_no_tps = (z(i) < front_head) || (z(i) > front_tail) || (r_grain_l >= r_inner_no);

                if is_exposed_no_tps
                    r_flow_no_tps = r_inner_no;
                    R_grain_local_no_tps = 0;
                else
                    r_flow_no_tps = r_grain_l;
                    R_grain_local_no_tps = log(r_inner_no / r_grain_l) / (2 * pi * k_grain * dz);
                end

                % Convection (using r_inner_no)
                u_local_no = mdot / (pi * r_flow_no_tps^2 * rho_g);
                Re_local_no = (rho_g * 2 * r_flow_no_tps * u_local_no) / mu_g;
                f_no = (1.82 * log10(Re_local_no) - 1.64)^ (-2);
                Nu_no = (f_no/8 * Re_local_no * Pr_g)/ (1.07 + 12.7 * sqrt(f_no/8) * (Pr_g^(2/3) - 1));
                h_int_no = Nu_no * k_g / (2 * r_flow_no_tps);

                % Resistance and Fluxes (Using No-TPS Radii)
                R_conv_i_no = 1 / (h_int_no * 2 * pi * r_flow_no_tps * dz);
                R_case_no   = log(r_outer_no / r_inner_no) / (2 * pi * k_case * dz);
                R_conv_e_no = 1 / (10 * 2 * pi * r_outer_no * dz);

                Q_in_no_tps      = (T_cc_i - T_case_int_no_tps(t_idx, i)) / (R_conv_i_no + R_grain_local_no_tps);
                Q_case_to_ext_no = (T_case_int_no_tps(t_idx, i) - T_case_ext_no_tps(t_idx, i)) / R_case_no;
                Q_env_no         = (T_case_ext_no_tps(t_idx, i) - T_amb) / R_conv_e_no;

                % Update temperatures (No TPS)
                T_case_int_no_tps(t_idx+1, i) = T_case_int_no_tps(t_idx, i) + ...
                    ((Q_in_no_tps - Q_case_to_ext_no) / (0.5 * M_case_seg_no_tps * cp_case)) * dt;
                T_case_ext_no_tps(t_idx+1, i) = T_case_ext_no_tps(t_idx, i) + ...
                    ((Q_case_to_ext_no - Q_env_no) / (0.5 * M_case_seg_no_tps * cp_case)) * dt;
            end
        end

        %% -------- PLOTS --------

        figure('Color', 'w', 'Name', 'Thermal Analysis');

        subplot(1,3,1);
        imagesc(time_vec, z, T_iso_int');
        set(gca, 'YDir', 'normal'); colormap('hot'); colorbar;
        title('Grain/Insulation Interface [K]');
        xlabel('Time [s]'); ylabel('Axial Position z [m]');

        subplot(1,3,2);
        imagesc(time_vec, z, T_case_int');
        set(gca, 'YDir', 'normal'); colormap('hot'); colorbar;
        title('Insulation/Case Interface [K]');
        xlabel('Time [s]'); ylabel('Axial Position z [m]');

        subplot(1,3,3);
        imagesc(time_vec, z, T_case_ext');
        set(gca, 'YDir', 'normal'); colormap('hot'); colorbar;
        title('External Case Surface [K]');
        xlabel('Time [s]'); ylabel('Axial Position z [m]');

        % Max temperatures
        T_max_metal = max(max(T_case_int));

        fprintf('\n--- Thermal Analysis Case%d without Isolation ---\n', ii);
        fprintf('Max case inside temperature = %.1f K  ', max(max(T_case_int_no_tps)));
        fprintf('Max case outside temperature = %.1f K  ', max(max(T_case_ext_no_tps)));
        fprintf('Critical temperature = %.1f K\n', Case_Metal.Tmax);

        fprintf('\n--- Thermal Analysis Case%d with Isolation ---\n', ii);
        fprintf('Max insulation temperature = %.1f K  ', max(max(T_iso_int)));
        fprintf('Max case temperature = %.1f K  ', T_max_metal);
        fprintf('Critical temperature = %.1f K\n\n', Case_Metal.Tmax);

        if T_max_metal >= Case_Metal.Tmax
            warning('Case temperature exceeds allowable limit');
        end

        if M_cc >= 0.07
            warning('Chamber Mach number exceeds 0.07');
        end

        %% -------- STORE RESULTS --------

        if r_p > 0 && L_g > 0 && r_cc > r_p

            E_tmp = struct();

            E_tmp.Case           = Case(ii);
            E_tmp.AP             = Case(ii).AP;
            E_tmp.HTPB           = Case(ii).HTPB;
            E_tmp.Tcc            = T_cc_i;
            E_tmp.rb             = rb;
            E_tmp.tb             = t_b;
            E_tmp.I_sp           = I_sp;
            E_tmp.c_star         = c_star;
            E_tmp.c_T            = c_T;
            E_tmp.rho_p          = rho_grain;
            E_tmp.m_p            = m_sp;
            E_tmp.mdot           = mdot;
            E_tmp.M_cc           = M_cc;
            E_tmp.Me             = Me;

            % Geometry
            E_tmp.Geom.r_p       = r_p;
            E_tmp.Geom.r_cc      = r_cc;
            E_tmp.Geom.web       = web;
            E_tmp.Geom.r_t       = r_t;
            E_tmp.Geom.r_e       = r_e;
            E_tmp.Geom.epsilon_e = eps;
            E_tmp.Geom.epsilon_cc= eps_cc;
            E_tmp.Geom.L_g       = L_g;
            E_tmp.Geom.N         = N_grains;
            E_tmp.Geom.L_cc      = L_cc;
            E_tmp.Geom.t_cc      = th_case;
            E_tmp.Geom.t_iso     = th_TPC;
            E_tmp.Geom.L_div     = L_div;
            E_tmp.Geom.L_conv    = L_conv;
            E_tmp.Geom.L_tot     = L_cc + L_conv + L_div;
            E_tmp.Geom.b_f       = b_f;
            E_tmp.Geom.LD        = LD;
            E_tmp.Geom.V_f       = V_f;

            MotorEngine(jj) = E_tmp;
            jj = jj + 1;
        end
    end
end
end