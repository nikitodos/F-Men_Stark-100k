function refined = refine_cooling_support(T_bulk_p2, dP_channel, mdot_w, t_op, rho_w, g0, P_a, ...
                                           m_water_total, V_tank_water, P_tank_req, m_tank_shell)
% REFINE_COOLING_SUPPORT   Refined pressure cascade, tank sizing, mass roll-up
%   based on correct ONB criterion (bulk fluid temperature, not wall temperature).
%
%   Inputs (from main workspace after thermal analysis):
%       T_bulk_p2      - coolant bulk temperature profile [K], vector
%       dP_channel     - jacket friction pressure drop [Pa]
%       mdot_w         - coolant mass flow rate [kg/s]
%       t_op           - total operation time (t_burn + 4 s) [s]
%       rho_w          - coolant density at mean bulk T [kg/m³]
%       g0             - standard gravity [m/s²]
%       P_a            - ambient pressure [Pa]
%       m_water_total  - Rev7 water mass [kg] (for comparison)
%       V_tank_water   - Rev7 tank volume [m³] (for comparison)
%       P_tank_req     - Rev7 tank pressure [Pa] (for comparison)
%       m_tank_shell   - Rev7 shell mass [kg] (for comparison)
%
%   Output (struct 'refined'):
%       .P_tank_req        - refined water tank delivery pressure [Pa]
%       .m_water_load      - consumable water mass (with SF) [kg]
%       .V_water_tank      - internal volume [m³]
%       .m_shell_water     - water tank shell mass [kg]
%       .t_wall_water      - wall thickness [m]
%       .m_PG              - N₂ pressurant mass [kg]
%       .V_PG_tank         - N₂ tank internal volume [m³]
%       .m_shell_PG        - N₂ tank shell mass [kg]
%       .t_wall_PG         - N₂ tank wall thickness [m]
%       .m_dry_total       - total dry mass (shells only) [kg]
%       .m_wet_total       - total wet mass (water + N₂ + shells) [kg]
%       .P_out_jacket      - jacket exit pressure [Pa]
%       .P_in_jacket       - jacket inlet pressure [Pa]
%       .MAWP_water        - MAWP of water tank [Pa]
%       .P_design_water    - design pressure (MAWP + static head) [Pa]
%       .expansion_ratio   - N₂ blowdown expansion ratio [-]
%       .T_PG_f            - final N₂ temperature [K]
%       .cond_margin_PG    - condensation margin (T_f - T_b,N2) [K]

    %% 17. Pressure cascade (correct ONB criterion)
    dP_feed_total  = 2.5e5;          % Pa – lumped feed‑line allowance
    dT_boil_margin = 10.0;           % K   – subcooling margin (Sutton & Biblarz)

    T_bulk_max    = max(T_bulk_p2);                       % K
    T_target_boil = T_bulk_max + dT_boil_margin;          % K
    P_min_boil    = P_sat_water(T_target_boil);           % Pa

    refined.P_out_jacket = max(P_min_boil, P_a + 0.1e5);   % Pa
    refined.P_in_jacket  = refined.P_out_jacket + dP_channel;
    refined.P_tank_req   = refined.P_in_jacket + dP_feed_total;

    %% 18. Water tank (SS316L, Moss & Basic [10])
    SF_m_water = 1.10;
    SF_V_water = 1.10;
    m_water_req = SF_m_water * mdot_w * t_op;             % kg
    V_water_req = SF_V_water * m_water_req / rho_w;       % m³
    refined.V_water_tank = V_water_req;
    refined.m_water_load = rho_w * refined.V_water_tank;  % consumable water mass

    LoD_water = 2.0;
    R_water_tank = ( refined.V_water_tank / (2*pi*LoD_water) )^(1/3);
    D_water_tank = 2 * R_water_tank;
    L_water_tank = LoD_water * D_water_tank;

    S_allow_water  = 115e6;        % Pa (SA-240-316L, Moss Fig.11-1)
    eta_weld_water = 0.85;
    rho_shell_water = 8000;        % kg/m³
    CA_water       = 0.25e-3;      % m
    f_MAWP         = 1.25;

    refined.MAWP_water     = f_MAWP * refined.P_tank_req;
    P_hydro_water          = rho_w * g0 * L_water_tank;
    refined.P_design_water = refined.MAWP_water + P_hydro_water;

    t_c_water  = refined.P_design_water * R_water_tank / ...
                 (S_allow_water*eta_weld_water - 0.6*refined.P_design_water);
    t_l_water  = refined.P_design_water * R_water_tank / ...
                 (2*S_allow_water*eta_weld_water + 0.4*refined.P_design_water);
    t_wall_water = max([t_c_water, t_l_water]) + CA_water;
    refined.t_wall_water = max(t_wall_water, 1e-3);

    V_shell_water = 2*pi*R_water_tank*refined.t_wall_water*L_water_tank + ...
                    4*pi*R_water_tank^2*refined.t_wall_water;
    refined.m_shell_water = rho_shell_water * V_shell_water;

    %% 19. Pressurant N₂ (isentropic blowdown) and its tank
    MM_PG    = 28.014e-3;          % kg/mol
    gamma_PG = 1.40;
    T_b_PG   = 77.4;               % K (NIST)
    R_u_gas  = 8.3144626;
    R_s_PG   = R_u_gas / MM_PG;

    P_PG_i   = 300e5;              % Pa – 300 bar initial
    T_PG_i   = 293.15;             % K

    P_PG_f   = refined.MAWP_water;
    T_PG_f   = T_PG_i * (P_PG_f / P_PG_i) ^ ((gamma_PG - 1) / gamma_PG);
    refined.T_PG_f = T_PG_f;

    expansion_ratio = (P_PG_f / P_PG_i) * (T_PG_i / T_PG_f);
    refined.expansion_ratio = expansion_ratio;

    V_PG_final = refined.V_water_tank / (1 - expansion_ratio);
    refined.V_PG_tank = V_PG_final - refined.V_water_tank;
    refined.m_PG = P_PG_f * V_PG_final / (R_s_PG * T_PG_f);

    refined.cond_margin_PG = T_PG_f - T_b_PG;

    % PG tank (SA-516-70 carbon steel)
    S_allow_PG   = 138e6;          % Pa (Moss Fig.11-1 curve1)
    eta_weld_PG  = 0.85;
    rho_PG_shell = 7850;           % kg/m³
    LoD_PG       = 2.0;
    P_design_PG  = f_MAWP * P_PG_i;

    R_PG = ( refined.V_PG_tank / (2*pi*LoD_PG) )^(1/3);
    D_PG = 2 * R_PG;
    L_PG = LoD_PG * D_PG;
    t_c_PG   = P_design_PG * R_PG / (S_allow_PG*eta_weld_PG - 0.6*P_design_PG);
    t_l_PG   = P_design_PG * R_PG / (2*S_allow_PG*eta_weld_PG + 0.4*P_design_PG);
    t_wall_PG = max([t_c_PG, t_l_PG, 1e-3]);
    refined.t_wall_PG = t_wall_PG;

    V_shell_PG = 2*pi*R_PG*t_wall_PG*L_PG + 4*pi*R_PG^2*t_wall_PG;
    refined.m_shell_PG = rho_PG_shell * V_shell_PG;

    %% 20. Mass roll‑up
    refined.m_dry_total = refined.m_shell_water + refined.m_shell_PG;
    refined.m_wet_total = refined.m_dry_total + refined.m_water_load + refined.m_PG;

    % Store Rev7 values for comparison (optional, but useful for print)
    refined.m_water_total_rev7 = m_water_total;
    refined.V_tank_water_rev7  = V_tank_water;
    refined.P_tank_req_rev7    = P_tank_req;
    refined.m_tank_shell_rev7  = m_tank_shell;
end