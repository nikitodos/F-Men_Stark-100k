%% STARK-100k - INTEGRATED MOTOR, THERMAL & FEED-SYSTEM ANALYSIS
%  Politecnico di Milano - Space Propulsion - A.Y. 2025-26
%
%  AUTHORS:
%    THE "F-MEN: DAYS OF FUTURE BLAST" TEAM
%      D'Aloisio Giovanni Nicola
%      Gattone Marco
%      Pacheco Nikolaas Valentin
%      Rossi Simone
%      Tempesta Tommaso Elia
%
%  INPUTS:
%    None (engine properties defined inside the script)
%
%  OUTPUTS:
%    Fig 01  Burning Rate Fit (Saint-Robert / Vieille law)
%    Fig 02  SRM Normal Thickness — preliminary render
%    Fig 03  Combustion Chamber Pressure (Firing)
%    Fig 04  Burning Rate & relative error (Firing)
%    Fig 05  Thrust Profile & relative error (Firing)
%    Fig 06  Nozzle Isentropic Flow (Mach, T, P)
%    Fig 07  Nozzle Transient Heatmap (Gas/Coating, Coating/Metal, External)
%    Fig 08  Bartz HTC & Heat Flux along jacket + Water bulk temperature
%    Fig 09  Monte Carlo Convergence
%    Fig 10  Uncertainty Distributions (Isp, Itot, pc, t_burn)
%    Fig 11  SRM Normal Thickness — final render (with structural t_nozzle)
%    Fig 12  Feed-line Pressure Budget
%
%  HELPER FUNCTIONS REQUIRED IN THE SAME FOLDER:
%   bartz_curv.m            — Bartz (1957) HTC with curvature correction
%   conductivity_water.m    — liquid water thermal conductivity fit
%   cp_water.m              — liquid water specific heat fit
%   darcy_friction.m        — Darcy friction (laminar / Colebrook-White)
%   density_water.m         — liquid water density fit
%   Firing.m                — quasi-steady SRM transient simulation
%   mach_from_area.m        — area-Mach bisection solver
%   P_sat_water.m           — Antoine equation for water saturation pressure
%   refine_cooling_support.m — refined pressure cascade, tank sizing, mass
%                              roll-up (Sections 17-20, Nikolas)
%   render_engine_section_2D.m — 2-D axisymmetric section plot
%   sigma_bartz.m           — Bartz sigma correction factor
%   solve_geometry_fzero2.m — grain geometry solver (fzero bisection)
%   SolidUnitDesign.m       — parametric SRM design loop
%   thermal_chain.m         — iterative thermal chain at each nozzle station
%   viscosity_water.m       — liquid water dynamic viscosity fit
%
%  References:
%    [01] NASA SP-8039 (1971) — gas emissivity correlations
%    [02] Modest (2013), Radiative Heat Transfer, 3rd ed., Academic Press
%    [03] Rohsenow, Hartnett & Cho (1998), Handbook of Heat Transfer, 3rd ed.
%    [04] Churchill & Chu (1975), Int. J. Heat Mass Transfer 18:1323-1329
%    [05] Incropera et al. (2011), Fundamentals of HMT, 7th ed., Wiley
%    [06] Barabash et al. (2010), ITER Materials Properties Handbook,
%         CuCrZr database
%    [07] Marple et al. (2007), J. Therm. Spray Technol. 16:791
%    [08] Bartz (1957), Jet Propulsion 27(7):558-566
%    [09] Gnielinski (1976), Int. Chem. Eng. 16:359
%    [10] Moss & Basic (2013), Pressure Vessel Design Manual, 4th ed.,
%         Elsevier — used for ASME BPVC VIII-1 wall-thickness rules and
%         material allowables (SA-240-316L, SA-516-70)
%    [11] Sutton G.P. & Biblarz O. (2016), Rocket Propulsion Elements,
%         9th ed., Wiley — Sec. 8.5, ONB criterion for liquid coolants
%    [12] NIST Chemistry WebBook SRD 69 (Linstrom & Mallard, eds.) —
%         N2 normal boiling point, https://webbook.nist.gov/chemistry/
%    [13] S. Morelli et al., Surface & Coatings Technology 513 (2025) 132498
%         — YSZ APS TBC thermal conductivity data

clear; clc; close all;

tic;   % Start the timer

fprintf('\n');
fprintf('  ███████╗      ███╗   ███╗███████╗███╗   ██╗\n');
fprintf('  ██╔════╝      ████╗ ████║██╔════╝████╗  ██║\n');
fprintf('  █████╗  █████╗██╔████╔██║█████╗  ██╔██╗ ██║\n');
fprintf('  ██╔══╝  ╚════╝██║╚██╔╝██║██╔══╝  ██║╚██╗██║\n');
fprintf('  ██║           ██║ ╚═╝ ██║███████╗██║ ╚████║\n');
fprintf('  ╚═╝           ╚═╝     ╚═╝╚══════╝╚═╝  ╚═══╝\n');
fprintf('\n');
fprintf('  ██████╗  █████╗ ██╗   ██╗███████╗     ██████╗ ███████╗\n');
fprintf('  ██╔══██╗██╔══██╗╚██╗ ██╔╝██╔════╝    ██╔═══██╗██╔════╝\n');
fprintf('  ██║  ██║███████║ ╚████╔╝ ███████╗    ██║   ██║█████╗  \n');
fprintf('  ██║  ██║██╔══██║  ╚██╔╝  ╚════██║    ██║   ██║██╔══╝  \n');
fprintf('  ██████╔╝██║  ██║   ██║   ███████║    ╚██████╔╝██║     \n');
fprintf('  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝     ╚═════╝ ╚═╝    \n');
fprintf('\n');
fprintf('  ███████╗██╗   ██╗████████╗██╗   ██╗██████╗ ███████╗    ██████╗ ██╗      █████╗ ███████╗████████╗\n');
fprintf('  ██╔════╝██║   ██║╚══██╔══╝██║   ██║██╔══██╗██╔════╝    ██╔══██╗██║     ██╔══██╗██╔════╝╚══██╔══╝\n');
fprintf('  █████╗  ██║   ██║   ██║   ██║   ██║██████╔╝█████╗      ██████╔╝██║     ███████║███████╗   ██║   \n');
fprintf('  ██╔══╝  ██║   ██║   ██║   ██║   ██║██╔══██╗██╔══╝      ██╔══██╗██║     ██╔══██║╚════██║   ██║   \n');
fprintf('  ██║     ╚██████╔╝   ██║   ╚██████╔╝██║  ██║███████╗    ██████╔╝███████╗██║  ██║███████║   ██║   \n');
fprintf('  ╚═╝      ╚═════╝    ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚══════╝    ╚═════╝ ╚══════╝╚═╝  ╚═╝╚══════╝   ╚═╝  \n');
fprintf('\n');
fprintf('  STARK-100k  |  Politecnico di Milano  |  Space Propulsion A.Y. 2025-26 \n');
fprintf('          Team: D''Aloisio · Gattone · Pacheco · Rossi · Tempesta        \n');
fprintf('  -----------------------------------------------------------------------\n');
fprintf('  Initialising... \n\n');

% Uncomment the following lines to set a custom plotting style (e.g., for the report).
% Use only if Utopia font is installed on your system; otherwise, MATLAB will fall back to the default font.
% Install Utopia font from https://github.com/greyscalepress/font-specimens/blob/master/fonts/fonts-perso/Laurent/Utopia/utopia-regular.ttf

% % Put the LaTeX interpreter for the mathematics and the text
% set(0, 'defaultTextInterpreter', 'tex');
% set(0, 'defaultAxesTickLabelInterpreter', 'tex');
% set(0, 'defaultLegendInterpreter', 'tex');

% % Set Utopia as the font for the axes
% set(0, 'DefaultAxesFontName', 'Utopia');
% set(0, 'DefaultTextFontName', 'Utopia');

% % Set a dimensions coherent with the report
% set(0, 'DefaultAxesFontSize', 16);
% set(0, 'DefaultTextFontSize', 16);

% fprintf('=================================================================\n');
% fprintf('          STARK-100k  Integrated Motor, Thermal Analysis         \n');
% fprintf('      MDN250 liner + EPDM/CCF/KP TBC + active water cooling      \n');
% fprintf('=================================================================\n\n');

%% -----------------------------------------------------------------------
%% 1. SAINT-ROBERT BALLISTIC FIT
%% -----------------------------------------------------------------------
data.p_bar   = [10.1, 9.7, 10.3, 30.2, 31.0, 29.8, 50.2, 51.0, 50.3, ...
    70.2, 69.0, 69.8, 91.2, 89.1, 89.0];
data.rb_mm_s = [4.0,  3.8,  4.1,  5.6,  6.0,  5.7,  7.0,  7.2,  7.1, ...
    8.4,  8.3,  8.6,  8.8,  9.0,  9.2];

data.p_Pa  = data.p_bar  * 1e5;
data.rb_ms = data.rb_mm_s * 1e-3;

coeff_bar = polyfit(log(data.p_bar),  log(data.rb_mm_s), 1);
n_bar     = coeff_bar(1);
a_bar     = exp(coeff_bar(2));

coeff_SI  = polyfit(log(data.p_Pa),   log(data.rb_ms),   1);
n_SI      = coeff_SI(1);
a_SI      = exp(coeff_SI(2));

SS_res = sum((log(data.rb_mm_s) - polyval(coeff_bar, log(data.p_bar))).^2);
SS_tot = sum((log(data.rb_mm_s) - mean(log(data.rb_mm_s))).^2);
R2     = 1 - SS_res / SS_tot;

[a_SI_dossi, Inc_a_SI_dossi, n_SI_dossi, Inc_n_SI_dossi, R2_SI_dossi] = Uncertainty(data.p_Pa, data.rb_ms);
[a_bar_dossi, Inc_a_bar_dossi, n_bar_dossi, Inc_n_bar_dossi, R2_bar_dossi] = Uncertainty(data.p_bar, data.rb_mm_s);

Inc_a_SI_dossi_perc = Inc_a_SI_dossi / a_SI_dossi * 100;
Inc_n_SI_dossi_perc = Inc_n_SI_dossi / n_SI_dossi * 100;

Inc_a_bar_dossi_perc = Inc_a_bar_dossi / a_bar_dossi * 100;
Inc_n_bar_dossi_perc = Inc_n_bar_dossi / n_bar_dossi * 100;

fprintf('--- Saint-Robert Fit ---\n');
fprintf('a = %.4f mm/s/bar^n   n = %.4f   R2 = %.4f\n\n', a_bar, n_bar, R2);

figure('Name','Burning Rate Fit','NumberTitle','off');
subplot(2,1,1);
scatter(data.p_bar, data.rb_mm_s, 'filled'); hold on;
p_fit_vec = linspace(min(data.p_bar), max(data.p_bar), 100);
plot(p_fit_vec, a_bar * p_fit_vec.^n_bar, 'r-', 'LineWidth',2);
xlabel('Pressure [bar]'); ylabel('Burning rate [mm/s]');
legend('Experimental data','Vieille fit'); grid on;
title('Saint-Robert law fit');
subplot(2,1,2);
rel_err = (a_bar*data.p_bar.^n_bar - data.rb_mm_s) ./ data.rb_mm_s * 100;
scatter(data.p_bar, rel_err, 'filled');
xlabel('Pressure [bar]'); ylabel('Relative error [%]');
yline(0,'--r'); grid on; title('Relative percentage error');

%% -----------------------------------------------------------------------
%% 2. CEA GAS PROPERTIES — AP/HTPB 75/25 to 80/20
%% -----------------------------------------------------------------------
mf_AP   = [75,  76,  77,  78,  79,  80];
mf_HTPB = [25,  24,  23,  22,  21,  20];

rho_AP   = 1900;
rho_HTPB =  950;
k_AP = 0.405;
k_HTPB = 0.276;
cp_AP = 1255;
cp_HTPB = 1900;

Mm_gas  = [19.933, 20.329, 20.743, 21.174, 21.622, 22.087];
cp_gas  = [1.8997, 1.9005, 1.9006, 1.8998, 1.8983, 1.8961] * 1e3;
rho_gas = [9.1663, 8.8595, 8.5865, 8.3437, 8.1279, 7.9368];
mu_gas  = [0.63800, 0.66604, 0.69441, 0.72305, 0.75196, 0.78107] * 1e-4;
k_gas   = [2.2457,  2.3034,  2.3574,  2.4079,  2.4550,  2.4989] * 1e-1;
Pr_gas  = [0.5397,  0.5495,  0.5598,  0.5705,  0.5815,  0.5927];
T_cc    = [1830.79, 1931.87, 2033.85, 2136.52, 2239.64, 2342.92];

for ii = 1:6
    Case(ii).AP        = mf_AP(ii);
    Case(ii).HTPB      = mf_HTPB(ii);
    Case(ii).rho       = 1/(mf_AP(ii)/100/rho_AP + mf_HTPB(ii)/100/rho_HTPB);
    Case(ii).k         = 1/(mf_AP(ii)/100/k_AP + mf_HTPB(ii)/100/k_HTPB);
    Case(ii).cp        = 1/(mf_AP(ii)/100/cp_AP + mf_HTPB(ii)/100/cp_HTPB);
    Case(ii).T_cc      = T_cc(ii);
    Case(ii).Gas.Mm    = Mm_gas(ii);
    Case(ii).Gas.rho   = rho_gas(ii);
    Case(ii).Gas.cp    = cp_gas(ii);
    Case(ii).Gas.R     = 8314/Mm_gas(ii);
    Case(ii).Gas.cv    = cp_gas(ii) - 8314/Mm_gas(ii);
    Case(ii).Gas.gamma = cp_gas(ii)/(cp_gas(ii) - 8314/Mm_gas(ii));
    Case(ii).Gas.mu    = mu_gas(ii);
    Case(ii).Gas.k     = k_gas(ii);
    Case(ii).Gas.Pr    = Pr_gas(ii);
end

%% -----------------------------------------------------------------------
%% 3. MISSION REQUIREMENTS & MATERIAL PROPERTIES
%% -----------------------------------------------------------------------
Thrust_real = 100e3;
I_tot       = 2.5e6;
P_cc        = 70e5;
T_amb       = 293.15;
P_a         = 101325;
g0          = 9.80665;
alpha       = 10;
beta        = 35;
N_grains    = 1;

lambda     = (1 + cosd(alpha)) / 2;
Thrust_des = Thrust_real / lambda;

% --- MDN250 case material ---
MDN250.T_melt = 1430 + 273.15;
MDN250.Tmax_C = 460;
MDN250.Tmax   = MDN250.Tmax_C + 273.15;
MDN250.yield  = 1640*1e6;
MDN250.rho    = 7880;
MDN250.k      = 25.1;
MDN250.cp     = 481;

% Thermal Protection
th_TPC = 3e-3;

% EPDM + 50% CCF
EPDMCF.rr = 0.01;
EPDMCF.k = 0.453;
EPDMCF.cp = 1549;
EPDMCF.rho = 1293;

% EPDM + 30% KP
EPDMKP.rr = 0.015;
EPDMKP.k = 0.171;
EPDMKP.cp = 1778;
EPDMKP.rho = 1185;

% EPDM + 25% KP + 25% CCF
EPDMCFKP.rr = 0.005;
EPDMCFKP.k = 0.198;
EPDMCFKP.cp = 1973;
EPDMCFKP.rho = 1256;

% --- Nozzle liner material: MDN250 ---
Liner.k     = MDN250.k;
Liner.T_max = MDN250.Tmax;
Liner.yield = MDN250.yield;
Liner.rho   = MDN250.rho;
Liner.cp    = MDN250.cp;

% --- TBC: EPDM + 25%-CCF + 25%-KP ---
t_TBC   = th_TPC;
k_TBC   = EPDMCFKP.k;
rho_TBC = EPDMCFKP.rho;
cp_TBC  = EPDMCFKP.cp;
% R_TBC is defined after Sec.12 once t_liner (structural) is known.

% --- Radiation properties ---
eps_gas = 0;
sigma   = 5.670374419e-8;
eps_ext = 0;

% --- Coolant parameters ---
% These are kept for the energy-balance pass (Sec.9b) and feed-system sizing.
T_cool_in  = 291.15;   % K (18°C)
dT_target  = 20;       % K
T_out_max  = T_cool_in + dT_target;

% Channel geometry (used in Sec.9b energy balance and Sec.17-21)
N_ch     = 120;
b_ch     = 3e-3;
h_ch     = 2e-3;
A_ch     = b_ch * h_ch;
D_h      = 2*b_ch*h_ch / (b_ch + h_ch);
roughness = 3.2e-6;

%% -----------------------------------------------------------------------
%% 4. NOZZLE SIZING + GRAIN GEOMETRY PER CEA CASE
%% -----------------------------------------------------------------------
Design.Thrust = Thrust_des;
Design.P_cc = P_cc;
Design.I_tot = I_tot;
Design.alpha = alpha;
Design.beta = beta;
Design.P_e = P_a;
Design.a_bar = a_bar;
Design.n_bar = n_bar;
Design.a_SI = a_SI;
Design.n_SI = n_SI;

CF_rb = 1.02;

MotorEngine = SolidUnitDesign(Design, Case, MDN250, EPDMCFKP, th_TPC, N_grains, T_amb, CF_rb);

if ~exist('MotorEngine','var') || isempty(MotorEngine)
    error('No valid grain geometry found.');
end

%% -----------------------------------------------------------------------
%% 5. SCORING — compact design weights
%% -----------------------------------------------------------------------

w_cstar = 0.15; w_cT = 0.15; w_m = 0.20; w_L = 0.20; w_r = 0.30;

for ii = 1:length(MotorEngine)
    R_arr(ii)     = MotorEngine(ii).Geom.r_cc;
    L_arr(ii)     = MotorEngine(ii).Geom.L_tot;
    Cstar_arr(ii) = MotorEngine(ii).c_star;
    CT_arr(ii)    = MotorEngine(ii).c_T;
    M_arr(ii)     = MotorEngine(ii).m_p;
end

norm_fn = @(x) (x - min(x)) / max(max(x)-min(x), eps);

Score = w_cstar*(1-norm_fn(Cstar_arr)) + w_cT*(1-norm_fn(CT_arr)) + ...
    w_m*norm_fn(M_arr) + w_L*norm_fn(L_arr) + w_r*norm_fn(R_arr);

[Score_sorted, idx_s] = sort(Score);
Engine_sorted = MotorEngine(idx_s);
for ii = 1:length(Engine_sorted)
    Engine_sorted(ii).Score = Score_sorted(ii);
end

E1 = Engine_sorted(1);
fprintf('\n\n--- Best Configuration: AP/HTPB %d/%d ---\n', E1.AP, E1.HTPB);
fprintf('Isp=%.1f s   c*=%.1f m/s   cT=%.4f\n', E1.I_sp, E1.c_star, E1.c_T);
fprintf('r_cc=%.1f mm  r_t=%.1f mm  L_g=%.1f mm  m_p=%.2f kg\n\n', ...
    E1.Geom.r_cc*1e3, E1.Geom.r_t*1e3, E1.Geom.L_g*1e3, E1.m_p);

%% -----------------------------------------------------------------------
%% 6. 2D MOTOR RENDERING — PRELIMINARY
%% -----------------------------------------------------------------------

E1.Geom.t_nozzle = 0.02 * E1.Geom.r_cc;
[X_case_profile, R_case_profile] = render_engine_section_2D(E1, 0);

%% -----------------------------------------------------------------------
%% 7. TRANSIENT FIRING — Firing.m
%% -----------------------------------------------------------------------

Geometry.L_g      = E1.Geom.L_g;
Geometry.D_perf   = E1.Geom.r_p  * 2;
Geometry.D_out    = E1.Geom.r_cc * 2;
Geometry.D_exit   = E1.Geom.r_e  * 2;
Geometry.D_throut = E1.Geom.r_t  * 2;
Geometry.L_cc     = E1.Geom.L_cc;
Geometry.N        = N_grains;

Fire = Firing(Geometry, a_SI, n_SI, E1.c_star, E1.Case.rho, ...
    1e-3, 1e-3, E1.m_p, E1.Case.Gas, lambda, Thrust_real, E1.rb, true);

fprintf('--- Firing Results ---\n');
fprintf('t_burn  = %.1f s\n',    Fire.t_burn);
fprintf('T_mean  = %.1f kN\n',   Fire.T_mean);
fprintf('T_max   = %.1f kN\n',   Fire.T_max);
fprintf('Isp     = %.1f s\n',    Fire.Isp);
fprintf('Itot    = %.1f MN*s\n', Fire.I_tot*1e-6);
fprintf('c*      = %.1f m/s\n',  Fire.c_star);
fprintf('eta_c*  = %.3f\n',      Fire.eta_cstar);
fprintf('C_T     = %.4f\n\n',    Fire.C_T);

%% -----------------------------------------------------------------------
%% 8. NOZZLE ISENTROPIC FLOW RECONSTRUCTION (from render profile)
%% -----------------------------------------------------------------------

T0    = E1.Case.T_cc;
P0    = P_cc;
gamma = E1.Case.Gas.gamma;
R_g   = E1.Case.Gas.R;
N_pts = 300;

L_conv_val = E1.Geom.L_conv;
L_div_val  = E1.Geom.L_div;
L_noz      = L_conv_val + L_div_val;

x_end        = X_case_profile(end);
x_conv_start = x_end - L_noz;
x_th_abs     = x_conv_start + L_conv_val;

mask      = (X_case_profile >= x_conv_start) & (X_case_profile <= x_end);
X_noz_raw = X_case_profile(mask);
R_noz_raw = R_case_profile(mask);

[X_noz_raw, ia] = unique(X_noz_raw);
R_noz_raw       = R_noz_raw(ia);

R_th = interp1(X_noz_raw, R_noz_raw, x_th_abs, 'pchip');

x_noz = linspace(0, L_noz, N_pts);
r_noz = interp1(X_noz_raw - x_conv_start, R_noz_raw, x_noz, 'pchip');

A_noz  = pi * r_noz.^2;
A_th   = pi * R_th^2;
AR_noz = A_noz / A_th;

[~, i_th] = min(abs(x_noz - L_conv_val));

M_noz = zeros(1, N_pts);
for i = 1:N_pts
    M_noz(i) = mach_from_area(AR_noz(i), gamma, i > i_th);
end
M_noz(i_th) = 1.0;

cp_g  = gamma * R_g / (gamma-1);
T_noz = T0 ./ (1 + (gamma-1)/2 * M_noz.^2);
P_noz = P0 .* (T_noz/T0).^(gamma/(gamma-1));
V_noz = M_noz .* sqrt(gamma*R_g*T_noz);

fprintf('--- Nozzle Flow Table ---\n');
fprintf('%6s %8s %8s %10s %10s %10s\n','x[m]','A[m2]','M','T[K]','P[kPa]','V[m/s]');
idx_tbl = unique([1, round(linspace(1,N_pts,12)), N_pts]);
for i = idx_tbl
    fprintf('%6.3f %8.5f %8.4f %10.2f %10.3f %10.2f\n', ...
        x_noz(i), A_noz(i), M_noz(i), T_noz(i), P_noz(i)/1e3, V_noz(i));
end
fprintf('\n');

%% -----------------------------------------------------------------------
%% 9. NTTA - Nozzle Transient Thermal Analysis
%% -----------------------------------------------------------------------
% This is the thermal analysis of the nozzle wall.
%
% After the transient run, a supplementary quasi-steady Bartz pass is
% performed along the jacket span (Sec.9b) to obtain:
%   - hg_jacket   : Bartz HTC at end-of-burn wall temperature [W/(m²K)]
%   - q_jacket    : heat flux at end-of-burn [W/m²]
%   - T_bulk_water: coolant bulk temperature profile [K]
% These are used only for plotting and for the feed-system sizing.

Time        = Fire.Time;
Pcc_trace   = Fire.P_trace;
T_cc_val    = E1.Tcc;

Geo.X = x_noz;
Geo.R = r_noz;

Mat.t_coating   = th_TPC;
Mat.k_coating   = EPDMCFKP.k;
Mat.rho_coating = EPDMCFKP.rho;
Mat.cp_coating  = EPDMCFKP.cp;
Mat.t_metal     = E1.Geom.t_cc;
Mat.k_metal     = MDN250.k;
Mat.rho_metal   = MDN250.rho;
Mat.cp_metal    = MDN250.cp;

Gas        = E1.Case.Gas;

Air.rho    = 1.1614;
Air.mu     = 1.846e-5;
Air.k      = 0.02551;
Air.cp     = 1005;
Air.T      = T_amb;
Air.alpha  = Air.k / (Air.rho * Air.cp);
Air.nu     = Air.mu  / Air.rho;
Air.Pr     = Air.nu / Air.alpha;

T_cr_metal = MDN250.Tmax;

[DataOut] = Nozzle_Transient_Heatmap(Time, Pcc_trace, T_cc_val, Geo, Mat, Gas, Air, T_cr_metal);

% --- Extract peak temperatures from NTTA for summaries ---
T_TBC_max_NTTA    = max(DataOut.T_inner(:));   % K — max TBC hot-face temp
T_metal_max_NTTA  = max(DataOut.T_inter(:));   % K — max metal hot-face temp

%% -----------------------------------------------------------------------
%% 12. NOZZLE STRUCTURAL THICKNESS + FINAL RENDER
%% -----------------------------------------------------------------------
th_nozzle         = 2 * P_noz .* ( r_noz + th_TPC ) / (MDN250.yield);
E1.Geom.t_nozzle  = max(th_nozzle);
E1.Geom.t_coat    = th_TPC;
[X_case_profile, R_case_profile] = render_engine_section_2D(E1, 1);

% --- Define t_liner and R_liner here, after structural sizing ---
% They are also needed by the feed-system sections (Sec.17-21).
t_liner = E1.Geom.t_nozzle;
R_liner = t_liner / Liner.k;
R_TBC   = t_TBC   / k_TBC;          % [m²K/W] — also defined here for completeness

%% -----------------------------------------------------------------------
%% 9b. JACKET SPAN + QUASI-STEADY BARTZ PASS (for plots & feed sizing)
%% -----------------------------------------------------------------------
% Jacket span: from A/At = 2 subsonic to A/At = 2 supersonic.
% A quasi-steady Bartz evaluation is performed at the end-of-burn wall
% temperature (DataOut.T_inner at the last time step) to obtain physically
% consistent hg_jacket and q_jacket without assuming steady-state through
% the full thermal chain. The energy balance then gives T_bulk_water(x).

AR_jacket = 2.0;

i_sub = find(AR_noz(1:i_th) >= AR_jacket, 1, 'last');
if isempty(i_sub), i_sub = 1; end

i_sup_local = find(AR_noz(i_th:end) >= AR_jacket, 1, 'first');
if isempty(i_sup_local)
    i_sup = N_pts;
else
    i_sup = i_th - 1 + i_sup_local;
end

x_sub  = x_noz(i_sub);
x_sup  = x_noz(i_sup);
L_cool = x_sup - x_sub;

T_inner = DataOut.T;

idx_cool = i_sub:i_sup;
x_cool   = x_noz(idx_cool);
r_cool   = r_noz(idx_cool);
AR_cool  = AR_noz(idx_cool);
M_cool   = M_noz(idx_cool);
% T_cool_s = T_inner(idx_cool)';   % static gas temperature along jacket

fprintf('--- Jacket span ---\n');
fprintf('A/At=%.1f subsonic:   x = %.1f mm\n', AR_jacket, x_sub*1e3);
fprintf('A/At=%.1f supersonic: x = %.1f mm\n', AR_jacket, x_sup*1e3);
fprintf('Jacket length = %.1f mm  (%d stations)\n\n', L_cool*1e3, length(idx_cool));

% --- Bartz HTC at end-of-burn wall temperature ---
% Extract the TBC hot-face temperature at the last time step along the
% jacket stations. DataOut.T_inner is [Nx x Nt]; rows = axial stations.
T_wall_g = DataOut.T_inner(idx_cool, end)';   % [1 x N_jacket], K

Dt_noz  = 2 * R_th;
At_noz  = A_th;
mdot_g_eob = P_cc * At_noz / E1.c_star;   % mass flux at nominal MEOP
G_t_eob    = mdot_g_eob / At_noz;
rc_th      = 1.5 * R_th;

mu0_g  = E1.Case.Gas.mu;
k0_g   = E1.Case.Gas.k;
Pr_g   = E1.Case.Gas.Pr;
omega_B = 0.6;

dx_cool = x_noz(2) - x_noz(1);


mdot_w_g = 0.05; % Kg/s
T_w_in = 18 + 273.15;
mdot_w = mdot_w_g;

flag = 0;

while  flag == 0

T_w_bulk(1) = T_w_in;

for ii = 1 : length(T_wall_g)

    % Bartz sigma with end-of-burn wall temperature
    sigma_B_eob = sigma_bartz(T_wall_g(ii), T0, gamma, M_cool(ii), omega_B);

    % Bartz HTC
    hg_jacket(ii) = bartz_curv(Dt_noz, mu0_g, cp_g, Pr_g, G_t_eob, rc_th, AR_cool(ii), sigma_B_eob);

    % --- Water properties at mean bulk temperature ---
    rho_w = density_water(T_w_bulk(ii));
    mu_w  = viscosity_water(T_w_bulk(ii));
    k_w   = conductivity_water(T_w_bulk(ii));
    cp_w  = cp_water(T_w_bulk(ii));
    Pr_w  = cp_w * mu_w / k_w;

    % Channel velocity and Reynolds number (for feed-system sizing)
    u_w  = mdot_w / (N_ch * rho_w * A_ch);
    Re_w = rho_w * u_w * D_h / mu_w;

    if Re_w < 2300
    % Laminar regime — Shah & London (1978)
        alpha_star = min(b_ch, h_ch) / max(b_ch, h_ch);
        Nu_w   = 8.235 * (1 - 2.0421*alpha_star + 3.0853*alpha_star^2 ...
                        - 2.4765*alpha_star^3 + 1.0578*alpha_star^4 ...
                        - 0.1861*alpha_star^5);
        f_w    = 64 / Re_w;   % Laminar Darcy 
    else
        % Transitional — Gnielinski
        f_w   = darcy_friction(Re_w, roughness, D_h);
        Nu_w  = (f_w/8)*(Re_w-1000)*Pr_w / (1 + 12.7*sqrt(f_w/8)*(Pr_w^(2/3)-1));
    end
    h_c_w(ii) = Nu_w * k_w / D_h;

    R_conv_i = 1 ./ (hg_jacket(ii) .* (2 * pi * r_cool(ii) .* dx_cool));
    R_cond_c = log((r_cool(ii) + th_TPC) ./ r_cool(ii)) ./ (2 * pi * EPDMCFKP.k .* dx_cool);
    R_cond_m = log((r_cool(ii) + th_TPC + t_liner) ./ (r_cool(ii) + th_TPC)) ./ (2 * pi * Liner.k .* dx_cool);
    R_conv_e = 1 ./ (h_c_w(ii) .* (2 * pi * (r_cool(ii) + th_TPC + t_liner) .* dx_cool));

    % Heat flux: convection only (eps_gas = 0 for non-aluminised AP/HTPB)
    q_jacket(ii) = (T_wall_g(ii) - T_w_bulk(ii)) / (R_conv_i+R_cond_c+R_cond_m+R_conv_e);

    q_dot(ii) = q_jacket(ii) / r_cool(ii);
    dT_w = q_jacket(ii) / (mdot_w * cp_w);

    T_w_bulk(ii+1) = T_w_bulk(ii) + dT_w;

    dp_channel(ii) = f_w * (dx_cool/D_h) * 0.5 * rho_w * u_w^2;

end

if abs(T_w_bulk(end)-T_out_max) < 0.01
    flag = 1;
else
    
    mdot_w = mdot_w * T_w_bulk(end)/T_out_max;

end

end

dP_channel = sum(dp_channel);

% Worst-case station (maximum heat flux)
[q_max, i_wc] = max(q_dot);
x_wc   = x_cool(i_wc);
AR_wc  = AR_cool(i_wc);
hg_wc  = hg_jacket(i_wc);

fprintf('--- Quasi-steady Bartz pass (end-of-burn wall T from NTTA) ---\n');
fprintf('eps_gas = %.2f (non-aluminised AP/HTPB)\n', eps_gas);
fprintf('TBC: EPDM/CCF/KP  t=%.0f mm  k=%.3f W/mK\n', t_TBC*1e3, k_TBC);
fprintf('Worst-case: x=%.1f mm  A/At=%.3f  hg=%.0f W/(m^2K)  q=%.4f MW/m^2\n', ...
    x_wc*1e3, AR_wc, hg_wc, q_max/1e6);
fprintf('\n--- Coolant system (energy balance) ---\n');
fprintf('mdot_w            = %.4f kg/s\n', mdot_w);
fprintf('T_in              = %.1f C\n', T_w_bulk(1) - 273.15);
fprintf('T_out (actual)    = %.1f C  [target: %.0f C]\n', ...
    T_w_bulk(end) - 273.15, T_out_max - 273.15);
fprintf('N_channels        = %d,  size: %.0fx%.0f mm,  Dh=%.2f mm\n', ...
    N_ch, b_ch*1e3, h_ch*1e3, D_h*1e3);
fprintf('Pressure drop     = %.4f bar\n', dP_channel/1e5);

T_cool_out_actual = T_w_bulk(end);

%% -----------------------------------------------------------------------
%% 11. SATURATION PRESSURE CHECK & EXTERNAL HEAT LOSS
%% -----------------------------------------------------------------------
% Uses T_bulk_water from Sec.9b.

P_sat_out = P_sat_water(T_cool_out_actual);
fprintf('\nP_sat at T_out = %.1f C  ->  %.3f bar\n', ...
    T_cool_out_actual - 273.15, P_sat_out/1e5);
fprintf('  Coolant pressure must exceed this to avoid boiling.\n');

P_vent_margin = 0.1e5;
P_jacket_exit = max(P_sat_out * 1.5, P_a + P_vent_margin);
P_jacket_in   = P_jacket_exit + dP_channel;
feed_line_dP  = 0.5e5;
P_tank_req    = P_jacket_in + feed_line_dP;

% External natural convection + radiation
% We use the NTTA outer surface temperature at end-of-burn
T_surf_ext = DataOut.T_outer(idx_cool, end)';   % [1 x N_jacket], K
T_film_ext = (T_surf_ext + T_amb) / 2;
rho_air    = 1.1614 * (293.15 ./ T_film_ext);
mu_air     = 1.846e-5 * (T_film_ext/293.15).^0.7;
k_air      = 0.02551 * (T_film_ext/293.15).^0.82;
cp_air     = 1005;
alpha_air  = k_air ./ (rho_air * cp_air);
nu_air     = mu_air ./ rho_air;
Pr_air     = nu_air ./ alpha_air;
beta_air   = 1 ./ T_film_ext;
DT_ext     = T_surf_ext - T_amb;

Ra_nc  = 9.80665 * beta_air .* DT_ext * L_noz^3 ./ (nu_air .* alpha_air);
psi_nc = (1 + (0.492./Pr_air).^(9/16)).^(-16/9);
Nu_nc  = (0.825 + 0.387*(Ra_nc.*psi_nc).^(1/6)).^2;
h_nc   = Nu_nc .* k_air / L_noz;
q_rad_ext = eps_ext * sigma * (T_surf_ext.^4 - T_amb^4);
h_ext  = h_nc + q_rad_ext ./ max(DT_ext, 1);

r_noz_out  = r_cool + t_liner;
dA_ext_vec = 2*pi*r_noz_out .* sqrt(1 + gradient(r_noz_out,x_cool).^2) .* gradient(x_cool);
Q_ext = sum((h_ext .* DT_ext + q_rad_ext) .* dA_ext_vec);

fprintf('\n--- External heat loss (nat. conv. + radiation, Churchill-Chu) ---\n');
fprintf('Mean h_nc = %.1f W/(m^2K),  mean h_rad_ext = %.1f W/(m^2K)\n', ...
    mean(h_nc), mean(q_rad_ext./max(DT_ext,1)));

%% -----------------------------------------------------------------------
%% 13. COOLANT TANK SIZING — APPROXIMATE BASELINE
%% -----------------------------------------------------------------------
t_pre  = 2;  t_post = 2;
t_op   = Fire.t_burn + t_pre + t_post;
V_water = mdot_w * t_op / rho_w;

P_jacket_exit = P_sat_out * 1.5;
P_jacket_in   = P_jacket_exit + dP_channel;
feed_line_dP  = 0.5e5;

gamma_N2 = 1.4;
P_N2_i   = 30e5;
if P_tank_req >= P_N2_i
    warning('Required tank pressure exceeds N2 initial pressure. Increase P_N2_i.');
    P_N2_i = P_tank_req * 1.5;
end
V_tank_water = V_water / (1 - (P_tank_req/P_N2_i)^(1/gamma_N2));

rho_tank = 8000;
S_tank   = 137e6;
eta_weld = 0.85;
L_over_D = 2.0;
D_tank = (4*V_tank_water/(pi*L_over_D))^(1/3);
L_tank = L_over_D * D_tank;
P_design_tank = 1.25 * P_tank_req;
t_tank = P_design_tank * (D_tank/2) / (S_tank * eta_weld - 0.6*P_design_tank);
t_tank = max(t_tank, 1e-3);
m_tank_shell = rho_tank * (2*pi*(D_tank/2)*t_tank*L_tank + 4*pi*(D_tank/2)^2*t_tank);
m_water_total = rho_w * V_water;

fprintf('\n--- Coolant tank sizing (baseline — see Sec.17 for refined) ---\n');
fprintf('Operation time   = %.1f s  (t_burn=%.1f + 2+2 s buffer)\n', t_op, Fire.t_burn);
fprintf('Water volume     = %.1f L  (%.1f kg)\n', V_water*1e3, m_water_total);
fprintf('Tank volume      = %.1f L  (incl. pressurisation headroom)\n', V_tank_water*1e3);
fprintf('P_tank_req       = %.2f bar\n', P_tank_req/1e5);
fprintf('P_N2_initial     = %.0f bar\n', P_N2_i/1e5);
fprintf('D_tank           = %.0f mm   L_tank = %.0f mm\n', D_tank*1e3, L_tank*1e3);
fprintf('Wall thickness   = %.1f mm   (SS316L ASME UG-27)\n', t_tank*1e3);
fprintf('Tank shell mass  = %.1f kg\n', m_tank_shell);
fprintf('Total system mass (water+tank) = %.1f kg\n', m_water_total + m_tank_shell);

%% -----------------------------------------------------------------------
%% 14. PLOTS
%% -----------------------------------------------------------------------

% --- Nozzle isentropic flow ---
r_norm_noz = r_noz / max(r_noz);
figure('Name','Nozzle Isentropic Flow','NumberTitle','off','Position',[100 80 1200 500]);
vars_noz  = {M_noz, P_noz/1e3};
ylbls_noz = {'M [-]','P [kPa]'};
ttls_noz  = {'Mach at MEOP','Static Pressure at MEOP'};
for k = 1:2
    subplot(1,2,k);
    plot(x_noz*1e3, vars_noz{k},'LineWidth',1.8,'Color',[0.17 0.45 0.72]);
    yl = ylim; hold on;
    sil_top = yl(1) + r_norm_noz*(yl(2)-yl(1))*0.18;
    fill([x_noz*1e3,fliplr(x_noz*1e3)],[sil_top,fliplr(yl(1)*ones(1,N_pts))], ...
        [0.75 0.75 0.75],'FaceAlpha',0.4,'EdgeColor',[0.5 0.5 0.5]);
    plot(x_noz*1e3, vars_noz{k},'LineWidth',1.8,'Color',[0.17 0.45 0.72]);
    xline(L_conv_val*1e3,'--','Color',[0.4 0.4 0.4],'LineWidth',1, ...
        'Label','throat','LabelVerticalAlignment','bottom');
    hold off; xlabel('x [mm]'); ylabel(ylbls_noz{k}); title(ttls_noz{k}); grid on;
end

% --- Bartz HTC, heat flux, and water bulk temperature (jacket zone) ---
% Three subplots: HTC, heat flux, T_bulk_water
figure('Name','Bartz HTC, Heat Flux & Water Temperature (jacket zone)', ...
    'NumberTitle','off','Position',[100 80 800 700]);

subplot(4,1,1);
plot(x_cool*1e3, hg_jacket/1e3, 'r-', 'LineWidth', 2);
xline(x_noz(i_th)*1e3,'--','Color',[0.4 0.4 0.4],'LineWidth',1,'Label','throat');
xline(x_wc*1e3,'-.','Color',[0.8 0.2 0.2],'LineWidth',1.5,'Label','worst-case');
xlabel('x [mm]'); ylabel('h_g [kW/(m^2K)]');
title('Bartz HTC along nozzle (jacket zone)'); grid on;

subplot(4,1,2);
plot(x_cool*1e3, h_c_w/1e3, 'r-', 'LineWidth', 2);
xline(x_noz(i_th)*1e3,'--','Color',[0.4 0.4 0.4],'LineWidth',1,'Label','throat');
xline(x_wc*1e3,'-.','Color',[0.8 0.2 0.2],'LineWidth',1.5,'Label','worst-case');
xlabel('x [mm]'); ylabel('h_c_w [kW/(m^2K)]');
title('Water convenction coefficient along the nozzle'); grid on;

subplot(4,1,3);
plot(x_cool*1e3, q_dot/1e6, 'b-', 'LineWidth', 2);
xline(x_noz(i_th)*1e3,'--','Color',[0.4 0.4 0.4],'LineWidth',1,'Label','throat');
xline(x_wc*1e3,'-.','Color',[0.8 0.2 0.2],'LineWidth',1.5,'Label','worst-case');
xlabel('x [mm]'); ylabel('q [MW/m^2]');
title('Heat flux along nozzle (jacket zone)'); grid on;

subplot(4,1,4);
plot(x_cool*1e3, T_w_bulk(2:end) - 273.15, 'c-', 'LineWidth', 2);
xline(x_noz(i_th)*1e3,'--','Color',[0.4 0.4 0.4],'LineWidth',1,'Label','throat');
xline(x_wc*1e3,'-.','Color',[0.8 0.2 0.2],'LineWidth',1.5,'Label','worst-case');
yline(T_w_bulk(1) - 273.15, '--k', 'LineWidth', 1, ...
    'Label', sprintf('T_{in} = %.0f C', T_w_bulk(1)-273.15));
yline(T_w_bulk(end) - 273.15, '--r', 'LineWidth', 1, ...
    'Label', sprintf('T_{out} = %.1f C', T_w_bulk(end)-273.15));
xlabel('x [mm]'); ylabel('T_{bulk} [^{\circ}C]');
title('Water bulk temperature along jacket'); grid on;

%% -----------------------------------------------------------------------
%% 15. UNCERTAINTY ANALYSIS — Monte Carlo using Firing.m
%% -----------------------------------------------------------------------

N_MC = 2000;
rng(42);

sigma_a_rel = Inc_a_SI_dossi_perc / 100;
sigma_n     = Inc_n_SI_dossi;

p_max_MC   = zeros(N_MC, 1);
t_burn_MC  = zeros(N_MC, 1);
Isp_MC     = zeros(N_MC, 1);
Itot_MC    = zeros(N_MC, 1);

Geometry   = struct('L_g', E1.Geom.L_g, 'N', N_grains, ...
    'D_perf', E1.Geom.r_p*2, 'D_out', E1.Geom.r_cc*2, ...
    'D_exit', E1.Geom.r_e*2, 'D_throut', E1.Geom.r_t*2, ...
    'L_cc', E1.Geom.L_cc);
rho_prop   = E1.Case.rho;
m_prop     = E1.m_p;
Gas_mc     = E1.Case.Gas;
lambda_div = lambda;
T_nom      = Thrust_real;
rb_nom     = E1.rb;
c_star_nom = E1.c_star;

fprintf('Running Monte Carlo (N = %d) with Firing.m ...\n', N_MC);
for i = 1:N_MC
    a_i = a_SI * exp(sigma_a_rel * randn());
    n_i = n_SI + sigma_n * randn();
    try
        Eng = Firing(Geometry, a_i, n_i, c_star_nom, rho_prop, ...
            1e-3, 1e-3, m_prop, Gas_mc, lambda_div, T_nom, rb_nom, false);
        p_max_MC(i) = max(Eng.P_trace) / 1e5;
        t_burn_MC(i) = Eng.t_burn;
        Isp_MC(i)    = Eng.Isp;
        Itot_MC(i)   = Eng.I_tot;
    catch
        p_max_MC(i) = NaN;
        t_burn_MC(i) = NaN;
        Isp_MC(i)    = NaN;
        Itot_MC(i)   = NaN;
    end
end

valid = ~isnan(p_max_MC);
p_max_MC  = p_max_MC(valid);
t_burn_MC = t_burn_MC(valid);
Isp_MC    = Isp_MC(valid);
Itot_MC   = Itot_MC(valid);
N_valid = sum(valid);
fprintf('Valid samples: %d / %d\n', N_valid, N_MC);

MEOP        = mean(p_max_MC);
MEOP_std    = std(p_max_MC);
tburn       = mean(t_burn_MC);
tburn_std   = std(t_burn_MC);
Isp         = mean(Isp_MC);
Isp_std     = std(Isp_MC);
Itot        = mean(Itot_MC) / 1e6;
Itot_std    = std(Itot_MC) / 1e6;

fprintf('\n--- Uncertainty Analysis Results (N_valid = %d) ---\n', N_valid);
fprintf('Sources: a (log-normal, sigma=%.2f%%), n (normal, sigma=%.4f)\n', ...
    sigma_a_rel*100, sigma_n);
fprintf('c* ideal and constant = %.1f m/s\n', c_star_nom);
fprintf('MEOP (mean +/- std)     = %.2f +/- %.2f bar  (design P_cc = %.0f bar)\n', ...
    MEOP, MEOP_std, P_cc/1e5);
fprintf('t_b     (mean +/- std)  = %.2f +/- %.2f s    (nominal = %.2f s)\n', ...
    tburn, tburn_std, E1.tb);
fprintf('Isp     (mean +/- std)  = %.1f +/- %.1f s    (nominal = %.1f s)\n', ...
    Isp, Isp_std, Fire.Isp);
fprintf('Itot    (mean +/- std)  = %.2f +/- %.2f MN*s (nominal = %.2f MN*s)\n\n', ...
    Itot, Itot_std, Fire.I_tot/1e6);

% ----- Convergence monitoring -----
N_check  = round(linspace(50, N_valid, 100));
run_mean = zeros(size(N_check));
run_std  = zeros(size(N_check));
for k = 1:length(N_check)
    sub         = t_burn_MC(1:N_check(k));
    run_mean(k) = mean(sub);
    run_std(k)  = std(sub);
end

figure('Name','Monte Carlo Convergence','NumberTitle','off');
subplot(2,1,1);
plot(N_check, run_mean, 'b-', 'LineWidth', 1.5);
yline(mean(t_burn_MC), 'r--', 'Label', 'converged mean');
xlabel('N samples'); ylabel('Mean t_{burn} [s]');
title('MC convergence — mean burn time (Firing model)'); grid on;
subplot(2,1,2);
plot(N_check, run_std, 'r-', 'LineWidth', 1.5);
yline(std(t_burn_MC), 'b--', 'Label', 'converged std');
xlabel('N samples'); ylabel('Std t_{burn} [s]');
title('MC convergence — std'); grid on;

% ----- Uncertainty distributions -----
colorBlue   = [0 0.4470 0.7410];
colorOrange = [0.8500 0.3250 0.0980];
colorGreen  = [0.4660 0.6740 0.1880];
colorPurple = [0.4940 0.1840 0.5560];
colorMean = 'k-';
colorStd  = 'k--';
colorNom  = 'r-';

figure('Name','Uncertainty Distributions: Performance','NumberTitle','off');
subplot(1,2,1);
histogram(Isp_MC, 40, 'Normalization','pdf', 'FaceColor', colorBlue, 'EdgeColor', 'w');
xline(Isp, colorMean, 'LineWidth', 1.8, 'Label', sprintf('mean = %.1f s', Isp));
xline(Isp+Isp_std, colorStd, 'LineWidth', 1.2, 'Label', sprintf('+1sigma = %.1f', Isp+Isp_std));
xline(Isp-Isp_std, colorStd, 'LineWidth', 1.2, 'Label', sprintf('-1sigma = %.1f', Isp-Isp_std));
xline(Fire.Isp, colorNom, 'LineWidth', 1.5, 'Label', sprintf('nominal = %.1f s', Fire.Isp));
xlabel('Specific Impulse [s]'); ylabel('PDF');
title('I_{sp} distribution'); grid on;
subplot(1,2,2);
histogram(Itot_MC/1e6, 40, 'Normalization','pdf', 'FaceColor', colorOrange, 'EdgeColor', 'w');
xline(Itot, colorMean, 'LineWidth', 1.8, 'Label', sprintf('mean = %.2f MN*s', Itot));
xline(Itot+Itot_std, colorStd, 'LineWidth', 1.2, 'Label', sprintf('+1sigma = %.2f', Itot+Itot_std));
xline(Itot-Itot_std, colorStd, 'LineWidth', 1.2, 'Label', sprintf('-1sigma = %.2f', Itot-Itot_std));
xline(Fire.I_tot/1e6, colorNom, 'LineWidth', 1.5, 'Label', sprintf('nominal = %.2f MN*s', Fire.I_tot/1e6));
xlabel('Total Impulse [MN*s]'); ylabel('PDF');
title('I_{tot} distribution'); grid on;

figure('Name','Uncertainty Distributions: Chamber Conditions','NumberTitle','off');
subplot(1,2,1);
histogram(p_max_MC, 40, 'Normalization','pdf', 'FaceColor', colorGreen, 'EdgeColor', 'w');
xline(MEOP, colorMean, 'LineWidth', 1.8, 'Label', sprintf('mean = %.1f bar', MEOP));
xline(MEOP+MEOP_std, colorStd, 'LineWidth', 1.2, 'Label', sprintf('+1sigma = %.1f', MEOP+MEOP_std));
xline(MEOP-MEOP_std, colorStd, 'LineWidth', 1.2, 'Label', sprintf('-1sigma = %.1f', MEOP-MEOP_std));
xline(P_cc/1e5, colorNom, 'LineWidth', 1.5, 'Label', sprintf('design = %.0f bar', P_cc/1e5));
xlabel('Maximum Chamber Pressure [bar]'); ylabel('PDF');
title('P_{c,max} distribution'); grid on;
subplot(1,2,2);
histogram(t_burn_MC, 40, 'Normalization','pdf', 'FaceColor', colorPurple, 'EdgeColor', 'w');
xline(tburn, colorMean, 'LineWidth', 1.8, 'Label', sprintf('mean = %.1f s', tburn));
xline(tburn+tburn_std, colorStd, 'LineWidth', 1.2, 'Label', sprintf('+1sigma = %.1f', tburn+tburn_std));
xline(tburn-tburn_std, colorStd, 'LineWidth', 1.2, 'Label', sprintf('-1sigma = %.1f', tburn-tburn_std));
xline(E1.tb, colorNom, 'LineWidth', 1.5, 'Label', sprintf('nominal = %.1f s', E1.tb));
xlabel('Burn Time [s]'); ylabel('PDF');
title('t_b distribution'); grid on;

%% -----------------------------------------------------------------------
%% 16. FINAL DESIGN SUMMARY
%% -----------------------------------------------------------------------
% Temperatures sourced from NTTA (DataOut)
fprintf('=================================================================\n');
fprintf('                       FINAL DESIGN SUMMARY                      \n');
fprintf('             (MDN250 + EPDM/CCF/KP TBC + active water cooling)   \n');
fprintf('=================================================================\n');
fprintf('Propellant           AP/HTPB %d/%d\n',        E1.AP, E1.HTPB);
fprintf('Thrust               %.0f kN\n',              Thrust_real/1e3);
fprintf('Chamber pressure     %.0f bar\n',             P_cc/1e5);
fprintf('Throat diameter      %.1f mm\n',              2*R_th*1e3);
fprintf('Chamber radius       %.1f mm\n',              E1.Geom.r_cc*1e3);
fprintf('Port radius          %.1f mm\n',              E1.Geom.r_p*1e3);
fprintf('Grain length         %.1f mm\n',              E1.Geom.L_g*1e3);
fprintf('L_cc (incl. L_ad)    %.1f mm\n',              E1.Geom.L_cc*1e3);
fprintf('Isp (ideal)          %.1f s\n',               E1.I_sp);
fprintf('t_burn               %.1f s\n',               Fire.t_burn);
fprintf('T_mean               %.1f kN\n',              Fire.T_mean);
fprintf('Isp (Firing)         %.1f s\n',               Fire.Isp);
fprintf('eta_c*               %.3f\n',                 Fire.eta_cstar);
fprintf('MEOP (mean +/- std)  %.2f +/- %.2f bar\n',    MEOP, MEOP_std);
fprintf('t_burn (mean+/-std)  %.2f +/- %.2f s\n',      tburn, tburn_std);
fprintf('Isp (mean +/- std)   %.1f +/- %.1f s\n',      Isp, Isp_std);
fprintf('Itot (mean +/- std)  %.2f +/- %.2f MN*s\n',   Itot, Itot_std);
fprintf('-----------------------------------------------------------------\n');
fprintf('THERMAL — NTTA (transient, MDN250 + EPDM/CCF/KP TBC)\n');
fprintf('  [temperatures = max over full time history and axial space]\n');
fprintf('Jacket span          A/At=%.1f sub -> A/At=%.1f sup\n', AR_jacket, AR_jacket);
fprintf('Jacket length        %.1f mm\n',              L_cool*1e3);
fprintf('Worst-case x (q_max) x=%.1f mm  A/At=%.3f\n', x_wc*1e3, AR_wc);
fprintf('q_max (end-of-burn)  %.4f MW/m^2\n',          q_max/1e6);
fprintf('-----------------------------------------------------------------\n');
fprintf('COOLING JACKET (energy balance from NTTA q-profile)\n');
fprintf('mdot_w               %.4f kg/s\n',            mdot_w);
fprintf('T_out (actual)       %.1f C  [target %.0f C]\n', ...
    T_cool_out_actual - 273.15, T_out_max - 273.15);
fprintf('N_channels           %d x %.0fx%.0f mm\n',    N_ch, b_ch*1e3, h_ch*1e3);
fprintf('Pressure drop        %.4f bar\n',             dP_channel/1e5);
fprintf('-----------------------------------------------------------------\n');
fprintf('SUPPORT APPARATUS  (baseline — see Sec.17-20 for refined)\n');
fprintf('Water volume         %.1f L\n',               V_water*1e3);
fprintf('Tank volume          %.1f L  (SS316L)\n',     V_tank_water*1e3);
fprintf('Tank pressure        %.2f bar\n',             P_tank_req/1e5);
fprintf('Tank mass (shell)    %.1f kg\n',              m_tank_shell);
fprintf('Total mass (w+tank)  %.1f kg\n',              m_water_total+m_tank_shell);
fprintf('=================================================================\n');

%% -----------------------------------------------------------------------
%%  17-20. REFINED PRESSURE CASCADE, TANK SIZING & MASS ROLL-UP
%% -----------------------------------------------------------------------
fprintf('\n=================================================================\n');
fprintf('   Sections 17-20 — Refined pressure cascade, tanks, mass\n');
fprintf('=================================================================\n');

T_bulk_water = T_w_bulk;

refined = refine_cooling_support( ...
    T_bulk_water, dP_channel, mdot_w, t_op, rho_w, g0, P_a, ...
    m_water_total, V_tank_water, P_tank_req, m_tank_shell);

fprintf('\n--- Lumped feed-line pressure cascade ---\n');
fprintf('  dP_feed_total (lumped)  = %.3f bar\n', 2.5e5/1e5);
fprintf('  dP_channel              = %.4f bar\n', dP_channel/1e5);
fprintf('  T_bulk_max (coolant outlet) = %.1f C\n', ...
    max(T_bulk_water)-273.15);
fprintf('  P_sat(T_bulk_max + 10 K) = %.4f bar  [11]\n', ...
    P_sat_water(max(T_bulk_water)+10)/1e5);
fprintf('P_out_jacket    = %.3f bar   P_in_jacket = %.3f bar\n', ...
    refined.P_out_jacket/1e5, refined.P_in_jacket/1e5);
fprintf('P_tank_required = %.3f bar   [Approx. P_tank_req = %.3f bar]\n', ...
    refined.P_tank_req/1e5, P_tank_req/1e5);
fprintf('  dP_feed / P_tank_required = %.1f%%\n', ...
    100*2.5e5/refined.P_tank_req);

fprintf('\n--- Water tank mechanical integrity (Moss & Basic [10]) ---\n');
fprintf('Material: SS316L  S=115 MPa [10]  E=0.85 [10]  CA=0.25 mm\n');
fprintf('Water mass (SF=1.10)  = %.1f kg    [Approx. = %.1f kg]\n', ...
    refined.m_water_load, m_water_total);
fprintf('Tank volume (SF=1.10) = %.1f L     [Approx. = %.1f L]\n', ...
    refined.V_water_tank*1e3, V_tank_water*1e3);
fprintf('MAWP = 1.25*P_tank = %.3f bar\n', refined.MAWP_water/1e5);
fprintf('  t_wall = %.2f mm   m_shell_water = %.1f kg   [Approx. = %.1f kg]\n', ...
    refined.t_wall_water*1e3, refined.m_shell_water, m_tank_shell);

fprintf('\n--- Pressurant gas system (N2 isentropic blowdown) ---\n');
fprintf('Gas: N2   MM=28.014 g/mol   gamma=1.40   T_b(1 atm)=77.4 K  [12]\n');
fprintf('P_i = 300 bar  T_i = 293.15 K  -->  P_f = %.3f bar  T_f = %.1f K\n', ...
    refined.MAWP_water/1e5, refined.T_PG_f);
fprintf('Expansion ratio = %.4f  (feasibility: < 1)\n', refined.expansion_ratio);
fprintf('V_PG_tank = %.3f L   m_PG = %.4f kg\n', ...
    refined.V_PG_tank*1e3, refined.m_PG);
fprintf('Condensation margin = T_f - T_b = %.1f K  [12]', refined.cond_margin_PG);
if refined.cond_margin_PG > 50
    fprintf('  --> OK\n');
elseif refined.cond_margin_PG > 0
    fprintf('  --> MARGINAL\n');
else
    fprintf('  *** CONDENSATION RISK ***\n');
end
fprintf('PG tank (SA-516-70): t=%.2f mm  m_shell=%.2f kg\n', ...
    refined.t_wall_PG*1e3, refined.m_shell_PG);
fprintf('  P_design = 1.25*300 bar = 375 bar\n');

fprintf('\n--- System mass and volume roll-up ---\n');
fprintf('                             mass [kg]    volume [L]\n');
fprintf('  water (consumable)          %7.2f       %7.3f\n', ...
    refined.m_water_load, refined.V_water_tank*1e3);
fprintf('  water-tank shell SS316L     %7.2f\n', refined.m_shell_water);
fprintf('  N2 pressurant               %7.4f       %7.3f (at final state)\n', ...
    refined.m_PG, refined.V_PG_tank*1e3);
fprintf('  PG-tank shell SA-516-70     %7.2f\n', refined.m_shell_PG);
fprintf('  -------------------------------------------------\n');
fprintf('  TOTAL WET (refined)         %7.2f\n', refined.m_wet_total);
fprintf('  TOTAL DRY (refined)         %7.2f\n', refined.m_dry_total);
fprintf('\n  Approx. baseline (water+shell only) : %.1f kg\n', ...
    m_water_total + m_tank_shell);
fprintf('  Refined water+shell              : %.1f kg  (delta = %+.2f kg)\n', ...
    refined.m_water_load + refined.m_shell_water, ...
    (refined.m_water_load + refined.m_shell_water) - (m_water_total + m_tank_shell));
fprintf('=================================================================\n');

%% -----------------------------------------------------------------------
%% 21. FEED-LINE PHYSICAL DESIGN
%% -----------------------------------------------------------------------

fprintf('\n=================================================================\n');
fprintf('   Section 21 — Feed-line physical design\n');
fprintf('=================================================================\n');

pipe_mat.name    = 'SS316L';
pipe_mat.S_allow = 115e6;
pipe_mat.E_mod   = 193e9;
pipe_mat.rho     = 8000;

roughness_pipe = 1.5e-6; % Assumed

fl = feed_line_design(mdot_w, rho_w, mu_w, ...
    refined.P_tank_req, refined.P_in_jacket, ...
    T_cool_in, g0, roughness_pipe, pipe_mat);

fprintf('\n--- (A) Pipe sizing ---\n');
fprintf('Material          : %s  (S=%.0f MPa, E=%.0f GPa)  [16]\n', ...
    pipe_mat.name, pipe_mat.S_allow/1e6, pipe_mat.E_mod/1e9);
fprintf('Assumed length    : %.2f m\n', fl.L_pipe_assumed);
fprintf('Selected DN       : D_i = %.1f mm   D_e = %.1f mm   t = %.2f mm\n', ...
    fl.D_pipe*1e3, fl.D_ext*1e3, fl.t_pipe*1e3);
fprintf('  t_min (ASME B31.3 Eq.304.1.2) = %.2f mm   + CA 0.5 mm\n', ...
    fl.t_min_ASME*1e3);
fprintf('Pipe shell mass   : %.3f kg\n', fl.m_pipe);

fprintf('\n--- (B) Distributed pressure drop ---\n');
fprintf('Flow velocity     : %.3f m/s\n', fl.u_pipe);
fprintf('Reynolds number   : %.0f  (%s)\n', fl.Re_pipe, ...
    ternary_str(fl.Re_pipe > 4000, 'turbulent', 'laminar/transitional'));
fprintf('Darcy f           : %.5f\n', fl.f_pipe);
fprintf('dP straight pipe  : %.4f bar\n', fl.dP_friction/1e5);
fprintf('dP fittings       : %.4f bar  (K_tot=%.2f)\n', ...
    fl.dP_fittings/1e5, fl.K_fittings);
fprintf('dP TOTAL          : %.4f bar\n', fl.dP_total/1e5);
fprintf('Budget available  : %.4f bar\n', fl.dP_budget/1e5);
fprintf('Margin            : %.4f bar  (%.1f%% of budget)  ', ...
    fl.margin_Pa/1e5, fl.margin_pct);
if fl.margin_pct >= 20
    fprintf('--> OK\n');
elseif fl.margin_pct >= 0
    fprintf('--> MARGINAL (< 20%%)\n');
else
    fprintf('--> INSUFFICIENT\n');
end

figure('Name','Feed-line Pressure Budget','NumberTitle','off', ...
    'Color','w','Position',[150 150 1000 420]);
bar_labels = {'dP friction','dP fittings','dP TOTAL','Budget','Margin'};
bar_vals   = [fl.dP_friction, fl.dP_fittings, fl.dP_total, fl.dP_budget, fl.margin_Pa]/1e5;
bar_colors = [0.2 0.5 0.8; 0.4 0.7 0.4; 0.9 0.6 0.1; 0.6 0.6 0.6; 0.2 0.7 0.3];
b = bar(bar_vals, 'FaceColor','flat');
for kk = 1:5; b.CData(kk,:) = bar_colors(kk,:); end
set(gca,'XTickLabel',bar_labels,'XTickLabelRotation',20);
ylabel('Pressure [bar]');
title('Feed-line pressure budget');
grid on; box on;

%% -----------------------------------------------------------------------
%% EXTENDED FINAL SUMMARY
%% -----------------------------------------------------------------------
fprintf('\n=================================================================\n');
fprintf('               EXTENDED FINAL SUMMARY                            \n');
fprintf('          MDN250 + EPDM/CCF/KP TBC + water cooling (refined)     \n');
fprintf('=================================================================\n');
fprintf('--- Ballistics & Performance ---\n');
fprintf('Propellant           AP/HTPB %d/%d\n',       E1.AP, E1.HTPB);
fprintf('Thrust               %.0f kN\n',              Thrust_real/1e3);
fprintf('Chamber pressure     %.0f bar\n',             P_cc/1e5);
fprintf('Isp (ideal)          %.1f s\n',               E1.I_sp);
fprintf('Isp (Firing)         %.1f s\n',               Fire.Isp);
fprintf('eta_c*               %.3f\n',                 Fire.eta_cstar);
fprintf('t_burn               %.1f s\n',               Fire.t_burn);
fprintf('MEOP (mean +/- std)  %.2f +/- %.2f bar\n',    MEOP, MEOP_std);
fprintf('t_burn (mean+/-std)  %.2f +/- %.2f s\n',      tburn, tburn_std);
fprintf('Isp (mean +/- std)   %.1f +/- %.1f s\n',      Isp, Isp_std);
fprintf('Itot (mean +/- std)  %.2f +/- %.2f MN*s\n',   Itot, Itot_std);
fprintf('--- Geometry ---\n');
fprintf('Throat diameter      %.1f mm\n',              2*R_th*1e3);
fprintf('Chamber radius       %.1f mm\n',              E1.Geom.r_cc*1e3);
fprintf('Port radius          %.1f mm\n',              E1.Geom.r_p*1e3);
fprintf('Grain length         %.1f mm\n',              E1.Geom.L_g*1e3);
fprintf('--- Thermal (NTTA — max over time & axial space) ---\n');
fprintf('q_max (end-of-burn)  %.4f MW/m^2\n',          q_max/1e6);
fprintf('T_TBC hot face       %.0f K  (%.0f C)\n', ...
    T_TBC_max_NTTA, T_TBC_max_NTTA - 273.15);
fprintf('T_metal hot face     %.0f K  (%.0f C)  [MDN250 limit %.0f C]\n', ...
    T_metal_max_NTTA, T_metal_max_NTTA - 273.15, MDN250.Tmax_C);
fprintf('--- Refined Feed System ---\n');
fprintf('P_out_jacket         %.3f bar\n',             refined.P_out_jacket/1e5);
fprintf('P_in_jacket          %.3f bar\n',             refined.P_in_jacket/1e5);
fprintf('P_tank_required      %.3f bar\n',             refined.P_tank_req/1e5);
fprintf('Water mass (consm.)  %.1f kg\n',              refined.m_water_load);
fprintf('Water tank shell     %.1f kg\n',              refined.m_shell_water);
fprintf('N2 pressurant        %.4f kg\n',              refined.m_PG);
fprintf('PG tank shell        %.2f kg\n',              refined.m_shell_PG);
fprintf('System mass WET      %.2f kg\n',              refined.m_wet_total);
fprintf('System mass DRY      %.2f kg\n',              refined.m_dry_total);
fprintf('=================================================================\n');
fprintf('--- Feed Line (Sec.21, ASME B31.3) ---\n');
fprintf('Pipe (SS316L)        D_i=%.1f mm  t=%.2f mm  L=%.2f m\n', ...
    fl.D_pipe*1e3, fl.t_pipe*1e3, fl.L_pipe);
fprintf('Flow velocity        %.3f m/s   Re=%.0f\n',   fl.u_pipe, fl.Re_pipe);
fprintf('dP distributed       %.4f bar  (friction + fittings)\n', fl.dP_total/1e5);
fprintf('Pressure margin      %.1f%%  of budget\n',    fl.margin_pct);
fprintf('Pipe mass            %.3f kg\n',              fl.m_pipe);
fprintf('=================================================================\n');

elapsed = toc;
fprintf('\n  -----------------------------------------------------------------------\n');
fprintf('  Total execution time: %d min %05.2f s\n', floor(elapsed/60), mod(elapsed,60));
fprintf('  -----------------------------------------------------------------------\n\n');