function fl = feed_line_design(mdot_w, rho_w, mu_w, P_tank, P_in_jacket, ...
                                T_fluid, g0, roughness_pipe, material)
% FEED_LINE_DESIGN   Physical sizing of the cooling-water feed line.
%
%   fl = feed_line_design(mdot_w, rho_w, mu_w, P_tank, P_in_jacket, ...
%                         T_fluid, g0, roughness_pipe, material)
%
%   Performs three analyses on the straight run of pipe that connects the
%   water storage tank to the cooling-jacket inlet manifold:
%
%   (A) PIPE SIZING  — iterative diameter selection so that the distributed
%       Darcy-Weisbach pressure drop along the total feed length does not
%       exceed the available budget (P_tank - P_in_jacket - dP_fittings).
%       Wall thickness is sized to ASME B31.3 Process Piping, using
%       SS316L allowable stress.
%
%   (B) DISTRIBUTED PRESSURE DROP — Darcy-Weisbach with Colebrook-White
%       friction factor, broken down into:
%         * straight pipe friction
%         * minor losses (fittings: 2 elbows, 1 gate valve, 1 check valve)
%       The balance confirms that the chosen diameter leaves adequate margin
%       in the pressure cascade.
%
%   (C) WATER HAMMER  — Joukowsky instantaneous pressure surge and the
%       acoustic time scale of the line (round-trip wave travel time).
%       The pipe bulk modulus model (Korteweg equation) is used to account
%       for pipe-wall elasticity.  A sudden valve closure (worst case) is
%       assumed.  Per Wylie & Streeter [14] the design must satisfy:
%         dP_WH < SF_WH * P_design_pipe
%       where SF_WH = 0.20 (surge not to exceed 20 % of design pressure).
%
%   -----------------------------------------------------------------------
%   Inputs
%   -----------------------------------------------------------------------
%   mdot_w        mass flow rate of coolant         [kg/s]
%   rho_w         coolant density                   [kg/m³]
%   mu_w          coolant dynamic viscosity         [Pa·s]
%   P_tank        water-tank delivery pressure      [Pa]
%   P_in_jacket   jacket inlet pressure             [Pa]
%   T_fluid       coolant bulk temperature          [K]   (for printout)
%   g0            standard gravity                  [m/s²]
%   roughness_pipe pipe internal roughness          [m]   (default: 1.5e-6 m = drawn SS)
%   material      struct with fields:
%                   .name     string label
%                   .S_allow  allowable stress  [Pa]  (ASME B31.3, Table A-1)
%                   .E_mod    Young's modulus   [Pa]  (for Korteweg)
%                   .rho      density           [kg/m³]
%
%   -----------------------------------------------------------------------
%   Outputs  (struct fl)
%   -----------------------------------------------------------------------
%   fl.D_pipe          selected inner diameter             [m]
%   fl.t_pipe          wall thickness                      [m]
%   fl.L_pipe          assumed total pipe length           [m]
%   fl.u_pipe          mean flow velocity                  [m/s]
%   fl.Re_pipe         Reynolds number                     [-]
%   fl.f_pipe          Darcy friction factor               [-]
%   fl.dP_friction     straight-pipe friction drop         [Pa]
%   fl.dP_fittings     minor-loss drop                     [Pa]
%   fl.dP_total        total feed-line pressure drop       [Pa]
%   fl.dP_budget       available pressure budget           [Pa]
%   fl.margin_Pa       remaining margin  (budget - total)  [Pa]
%   fl.margin_pct      margin as % of budget               [%]
%   fl.c_wave          pressure-wave speed (Korteweg)      [m/s]
%   fl.t_round         acoustic round-trip time            [s]
%   fl.dP_WH           Joukowsky peak surge pressure       [Pa]
%   fl.P_design_pipe   pipe design pressure (ASME B31.3)   [Pa]
%   fl.WH_ratio        dP_WH / P_design_pipe               [-]
%   fl.WH_flag         string: 'OK' / 'MARGINAL' / 'CRITICAL'
%   fl.m_pipe          pipe mass (shell only)              [kg]
%
%   -----------------------------------------------------------------------
%   Physical model notes
%   -----------------------------------------------------------------------
%   Pipe length:  assumed 1.5 m straight run (tank → manifold).
%                 Conservative for a test-bench layout; adjust L_pipe_assumed
%                 below for the actual installation.
%
%   Minor losses (sum of K-factors, Idelchik [15]):
%     2 × 90° standard elbow  K = 0.75 each  → 1.50
%     1 × gate valve (open)   K = 0.13
%     1 × swing check valve   K = 2.50
%     1 × sharp-edge entrance K = 0.50
%     1 × exit into manifold  K = 1.00
%     Total K_fittings = 6.63
%
%   Korteweg wave speed (Halliwell 1963, see also Wylie & Streeter [14]):
%     c = sqrt( K_f / rho_w ) / sqrt( 1 + (K_f/E_mod)*(D/t) )
%   where K_f = bulk modulus of water (2.2 GPa at ~20 °C).
%
%   Joukowsky surge (instantaneous valve closure, worst case):
%     dP_WH = rho_w * c * u_pipe
%
%   ASME B31.3 wall thickness (straight pipe under internal pressure,
%   Eq. 304.1.2):
%     t_min = P_design * D / (2 * (S_allow * E_qual + P_design * Y))
%   where E_qual = 1.0 (seamless), Y = 0.4 (austenitic SS, T < 482°C).
%   A manufacturing/corrosion allowance of 0.5 mm is added.
%
%   See also: refine_cooling_support, darcy_friction

% =========================================================================
%  0.  Defaults and geometry assumptions
% =========================================================================
L_pipe_assumed = 1.5;          % m   — conservative test-bench pipe length
K_fittings     = 6.63;         % [-] — lumped minor-loss coefficient (see header)
K_f_water      = 2.2e9;        % Pa  — bulk modulus of liquid water at ~20°C
SF_design      = 1.5;          % ASME B31.3 design-pressure safety factor
Y_coeff        = 0.4;          % ASME B31.3 Eq.304.1.2 Y-coefficient (austenic SS)
E_qual         = 1.0;          % weld-quality factor (seamless)
CA_pipe        = 0.5e-3;       % m   — corrosion/manufacturing allowance
SF_WH          = 0.20;         % water-hammer surge limit (fraction of P_design)

% =========================================================================
%  A.  PIPE DIAMETER SELECTION
%      Iterate over a standard DN schedule; pick the smallest DN that
%      keeps the total pressure drop within the available budget.
% =========================================================================

% Available budget: tank pressure minus jacket-inlet pressure
dP_budget = P_tank - P_in_jacket;   % Pa

% Standard DN inner diameters (SS316L schedule 40S, approx.)
% DN10 DN15 DN20 DN25 DN32 DN40 DN50 DN65 DN80  [mm]
D_candidates = [10.3, 15.8, 20.9, 26.6, 35.1, 40.9, 52.5, 62.7, 77.9] * 1e-3; % m

D_pipe = NaN;
for k = 1:length(D_candidates)
    D_try  = D_candidates(k);
    A_try  = pi/4 * D_try^2;
    u_try  = mdot_w / (rho_w * A_try);
    Re_try = rho_w * u_try * D_try / mu_w;
    f_try  = darcy_friction(Re_try, roughness_pipe, D_try);
    dP_fr  = f_try * (L_pipe_assumed / D_try) * 0.5 * rho_w * u_try^2;
    dP_mi  = K_fittings * 0.5 * rho_w * u_try^2;
    dP_tot = dP_fr + dP_mi;
    if dP_tot < dP_budget * 0.80   % keep at least 20% margin
        D_pipe = D_try;
        break;
    end
end

% Fallback: if no standard DN fits, use the largest candidate
if isnan(D_pipe)
    D_pipe = D_candidates(end);
    warning(['FEED_LINE_DESIGN: no standard DN satisfies the 20%% margin ' ...
             'criterion.  Using DN80 (D_i=%.1f mm).  Consider increasing ' ...
             'P_tank or shortening the feed line.'], D_pipe*1e3);
end

% =========================================================================
%  B.  DISTRIBUTED PRESSURE DROP  (with selected diameter)
% =========================================================================
A_pipe  = pi/4 * D_pipe^2;
u_pipe  = mdot_w / (rho_w * A_pipe);
Re_pipe = rho_w * u_pipe * D_pipe / mu_w;
f_pipe  = darcy_friction(Re_pipe, roughness_pipe, D_pipe);

dP_friction = f_pipe * (L_pipe_assumed / D_pipe) * 0.5 * rho_w * u_pipe^2;
dP_fittings = K_fittings * 0.5 * rho_w * u_pipe^2;
dP_total    = dP_friction + dP_fittings;

margin_Pa  = dP_budget - dP_total;
margin_pct = 100 * margin_Pa / dP_budget;

% =========================================================================
%  C.  WALL THICKNESS  (ASME B31.3, Eq. 304.1.2)
% =========================================================================
P_design_pipe = SF_design * P_tank;   % Pa — design pressure for the pipe
t_min = P_design_pipe * (D_pipe/2) / ...
        (material.S_allow * E_qual + P_design_pipe * Y_coeff);
t_pipe = t_min + CA_pipe;
t_pipe = max(t_pipe, 1.0e-3);        % 1 mm manufacturing floor

D_ext   = D_pipe + 2*t_pipe;
L_pipe  = L_pipe_assumed;
m_pipe  = material.rho * pi/4 * (D_ext^2 - D_pipe^2) * L_pipe;

% =========================================================================
%  D.  WATER HAMMER  (Joukowsky + Korteweg)
% =========================================================================
% Korteweg wave speed
c_wave = sqrt(K_f_water / rho_w) / ...
         sqrt(1 + (K_f_water / material.E_mod) * (D_pipe / t_pipe));

% Acoustic round-trip time (2 × pipe length)
t_round = 2 * L_pipe / c_wave;   % s

% Joukowsky surge (instantaneous closure)
dP_WH = rho_w * c_wave * u_pipe;  % Pa

% Check against allowable surge
WH_ratio = dP_WH / P_design_pipe;
if WH_ratio <= SF_WH
    WH_flag = 'OK';
elseif WH_ratio <= 0.50
    WH_flag = 'MARGINAL';
else
    WH_flag = 'CRITICAL';
end

% =========================================================================
%  E.  PACK OUTPUT STRUCT
% =========================================================================
fl.D_pipe          = D_pipe;
fl.D_ext           = D_ext;
fl.t_pipe          = t_pipe;
fl.L_pipe          = L_pipe;
fl.u_pipe          = u_pipe;
fl.Re_pipe         = Re_pipe;
fl.f_pipe          = f_pipe;
fl.dP_friction     = dP_friction;
fl.dP_fittings     = dP_fittings;
fl.dP_total        = dP_total;
fl.dP_budget       = dP_budget;
fl.margin_Pa       = margin_Pa;
fl.margin_pct      = margin_pct;
fl.P_design_pipe   = P_design_pipe;
fl.t_min_ASME      = t_min;
fl.m_pipe          = m_pipe;
fl.c_wave          = c_wave;
fl.t_round         = t_round;
fl.dP_WH           = dP_WH;
fl.WH_ratio        = WH_ratio;
fl.WH_flag         = WH_flag;
fl.K_fittings      = K_fittings;
fl.L_pipe_assumed  = L_pipe_assumed;
fl.material        = material;

end