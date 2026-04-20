function [DataOut] = Nozzle_Transient_Heatmap(time, Pc_trace, Tc_const, Geo, Mat, GasProp, Air, T_cr_metal)
% NOZZLE_TRANSIENT_HEATMAP  Quasi-steady 1D transient thermal analysis of a nozzle.
%
%   DataOut = Nozzle_Transient_Heatmap(time, Pc_trace, Tc_const, Geo, ...
%                                      Mat, GasProp, Air, T_cr_metal)
%
%   Performs a transient thermal analysis of an axisymmetric nozzle
%   protected by an internal coating (TBC) and a metallic liner.  The
%   model is quasi‑steady on the gas side: at each time step the Bartz
%   heat transfer coefficient is recomputed from the instantaneous chamber
%   pressure.  The wall temperature evolution is obtained by explicit
%   time integration of the lumped thermal capacitances.
%
%   -----------------------------------------------------------------------
%   Inputs
%   -----------------------------------------------------------------------
%   time        [1 × Nt] time vector                            [s]
%   Pc_trace    [1 × Nt] chamber pressure history               [Pa]
%   Tc_const    chamber stagnation temperature (constant)       [K]
%   Geo         struct with fields:
%                 .X   [Nx × 1] axial coordinates               [m]
%                 .R   [Nx × 1] wall radius at each station     [m]
%   Mat         struct with fields:
%                 .t_coating     coating thickness              [m]
%                 .k_coating     coating thermal conductivity   [W/(m·K)]
%                 .rho_coating   coating density                [kg/m³]
%                 .cp_coating    coating specific heat          [J/(kg·K)]
%                 .t_metal       metal liner thickness          [m]
%                 .k_metal       metal thermal conductivity     [W/(m·K)]
%                 .rho_metal     metal density                  [kg/m³]
%                 .cp_metal      metal specific heat            [J/(kg·K)]
%   GasProp     struct with fields:
%                 .gamma   specific heat ratio                  [-]
%                 .R       gas constant                         [J/(kg·K)]
%                 .Pr      gas Prandtl number                   [-]
%                 .mu      gas dynamic viscosity                [Pa·s]
%                 .cp      gas specific heat                    [J/(kg·K)]
%   Air         struct with fields:
%                 .Pr   Prandtl number                          [-]
%                 .nu   kinematic viscosity                     [m²/s]
%                 .k    thermal conductivity                    [W/(m·K)]
%                 .T    ambient temperature                     [K]
%   T_cr_metal  critical temperature limit for the metal liner  [K]
%
%   -----------------------------------------------------------------------
%   Output
%   -----------------------------------------------------------------------
%   DataOut     struct with fields:
%                 .T_inner  [Nx × Nt] coating inner surface temp [K]
%                 .T_inter  [Nx × Nt] coating/metal interface    [K]
%                 .T_outer  [Nx × Nt] metal outer surface        [K]
%
%   -----------------------------------------------------------------------
%   Notes
%   -----------------------------------------------------------------------
%   - External heat transfer is modelled as natural convection over a
%     horizontal cylinder (Churchill-Chu correlation) and radiation to
%     ambient (assumed black surroundings).
%   - Gas radiation is NOT included in this transient model.
%
%   See also: bartz_curv, mach_from_area, sigma_bartz.

% --- 1. INITIALIZATION & DATA EXTRACTION ---
Nt = length(time);
Nx = length(Geo.X);
X = Geo.X(:);
R = Geo.R(:);
Pc_trace = Pc_trace(:)';
dt = time(2) - time(1);

% Gas properties
gamma = GasProp.gamma;
R_gas = GasProp.R;
Pr_gas = GasProp.Pr;
mu_gas = GasProp.mu;
cp_gas = GasProp.cp;

% Material properties
t_c = Mat.t_coating;
k_c = Mat.k_coating;
rho_c = Mat.rho_coating;
cp_c = Mat.cp_coating;

t_m = Mat.t_metal;
k_m = Mat.k_metal;
rho_m = Mat.rho_metal;
cp_m = Mat.cp_metal;

T_a = Air.T;


% --- 2. GEOMETRY & THROAT ANALYSIS ---
[R_th, i_th] = min(R);
A_th = pi * R_th^2;
Dt = 2 * R_th;

dX = diff(X);
dR = diff(R);
d2R = diff(dR./dX);

rc = 1 / abs(d2R(max(1,i_th-1)) / dX(max(1,i_th-1)));
if isinf(rc) || rc <= 0, rc = 1.5 * R_th; end


% --- 3. ISENTROPIC FLOW ---
AR = (pi .* R.^2) ./ A_th;
M = zeros(Nx, 1);

for i = 1:Nx
    if i < i_th,      M(i) = mach_from_area(AR(i), gamma, false);
    elseif i == i_th, M(i) = 1.0;
    else,             M(i) = mach_from_area(AR(i), gamma, true);
    end
end

T_aw = (Tc_const ./ (1 + (gamma-1)/2 .* M.^2)) .* (1 + Pr_gas^(1/3) * (gamma-1)/2 .* M.^2);


% --- 4. THERMAL CAPACITANCE & RESISTANCE ---
dz_seg = (max(X) - min(X)) / (Nx - 1);

% Radii definition for cylindrical calculation
r_gas_wall = R;
r_tbc_int  = r_gas_wall;
r_tbc_ext  = r_gas_wall + t_c;
r_met_ext  = r_tbc_ext + t_m;

r_no_int   = r_gas_wall;
r_no_ext   = r_gas_wall + t_m;

% Volumetric Capacitance (J/K) = rho * cp * Volume
Cap_c = rho_c * cp_c * (pi * (r_tbc_ext.^2 - r_tbc_int.^2) * dz_seg);
Cap_m = rho_m * cp_m * (pi * (r_met_ext.^2 - r_tbc_ext.^2) * dz_seg);
Cap_m_no = rho_m * cp_m * (pi * (r_no_ext.^2 - r_no_int.^2) * dz_seg);

% Temp matrices
T_inner = ones(Nx, Nt) * T_a;
T_inter = ones(Nx, Nt) * T_a;
T_outer = ones(Nx, Nt) * T_a;

T_inner_no = ones(Nx, Nt) * T_a;
T_outer_no = ones(Nx, Nt) * T_a;


% --- 5. TRANSIENT LOOP ---
fprintf('Running Transient Comparison Analysis...\n');

for t = 1:Nt-1

    P0 = max(Pc_trace(t), 101325);
    mdot = (A_th * P0 / sqrt(Tc_const)) * sqrt(gamma/R_gas * (2/(gamma+1))^((gamma+1)/(gamma-1)));
    G_t = mdot / A_th;

    % --- WITH COATING ---
    Tw_est = T_inner(:, t);
    sigma_B = sigma_bartz(Tw_est, Tc_const, gamma, M, 0.6);
    hg = bartz_curv(Dt, mu_gas, cp_gas, Pr_gas, G_t, rc, AR, sigma_B);

    % External Natural Convection
    dT_ext = max(T_outer(:, t) - T_a, 0.1);
    Ra = (9.81 * (1/T_a) .* dT_ext .* (2 .* R).^3) ./ (Air.nu^2) .* Air.Pr;
    h_ext = (0.53 .* Ra.^0.25) .* Air.k ./ (2 .* R);

    % Resistance cylindrical calculation
    R_conv_i = 1 ./ (hg .* (2 * pi * r_tbc_int .* dz_seg));
    R_cond_c = log(r_tbc_ext ./ r_tbc_int) ./ (2 * pi * k_c .* dz_seg);
    R_cond_m = log(r_met_ext ./ r_tbc_ext) ./ (2 * pi * k_m .* dz_seg);
    R_conv_e = 1 ./ (h_ext .* (2 * pi * r_met_ext .* dz_seg));

    % Heat Fluxes [Watts]
    Q_in = (T_aw - T_inner(:, t)) ./ R_conv_i;
    Q_cond_c = (T_inner(:, t) - T_inter(:, t)) ./ R_cond_c;
    Q_cond_m = (T_inter(:, t) - T_outer(:, t)) ./ R_cond_m;
    Q_out = (T_outer(:, t) - T_a) ./ R_conv_e;

    T_inner(:, t+1) = T_inner(:, t) + ((Q_in - Q_cond_c) ./ (0.5 * Cap_c)) * dt;
    T_inter(:, t+1) = T_inter(:, t) + ((Q_cond_c - Q_cond_m) ./ (0.5 * Cap_c + 0.5 * Cap_m)) * dt;
    T_outer(:, t+1) = T_outer(:, t) + ((Q_cond_m - Q_out) ./ (0.5 * Cap_m)) * dt;


    % --- NO COATING (DIRECT) ---
    Tw_est_no = T_inner_no(:, t);
    sigma_B_no = sigma_bartz(Tw_est_no, Tc_const, gamma, M, 0.6);
    hg_no = bartz_curv(Dt, mu_gas, cp_gas, Pr_gas, G_t, rc, AR, sigma_B_no);

    % External Convection (shared)
    dT_ext_no = max(T_outer_no(:, t) - T_a, 0.1);
    Ra_no = (9.81 * (1/T_a) .* dT_ext_no .* (2 .* R).^3) ./ (Air.nu^2) .* Air.Pr;
    h_ext_no = (0.53 .* Ra_no.^0.25) .* Air.k ./ (2 .* R);

    % Resistance cylindrical calculation
    R_conv_i_no = 1 ./ (hg_no .* (2 * pi * r_no_int .* dz_seg));
    R_cond_m_no = log(r_no_ext ./ r_no_int) ./ (2 * pi * k_m .* dz_seg);
    R_conv_e_no = 1 ./ (h_ext_no .* (2 * pi * r_no_ext .* dz_seg));

    % Heat Fluxes [Watts]
    Q_in_no = (T_aw - T_inner_no(:, t)) ./ R_conv_i_no;
    Q_cond_no = (T_inner_no(:, t) - T_outer_no(:, t)) ./ R_cond_m_no;
    Q_out_no = (T_outer_no(:, t) - T_a) ./ R_conv_e_no;

    T_inner_no(:, t+1) = T_inner_no(:, t) + ((Q_in_no - Q_cond_no) ./ (0.5 * Cap_m_no)) * dt;
    T_outer_no(:, t+1) = T_outer_no(:, t) + ((Q_cond_no - Q_out_no) ./ (0.5 * Cap_m_no)) * dt;

end

% --- 6. RESULTS & CRITICAL TEMP CHECK ---
max_T_coat  = max(T_inner(:,end));
max_T_metal = max(T_inter(:,end));

max_T_metal_no = max(T_inner_no(:,end));

% Calculate safety margins
margin_metal = T_cr_metal - max_T_metal;
margin_metal_no = T_cr_metal - max_T_metal_no;

fprintf('\n--- Nozzle Thermal Analysis Results without Coating ---\n');
fprintf('Max Metal Temp:   %.1f K | Limit: %.1f K | Margin: %.1f K\n', ...
    max_T_metal_no, T_cr_metal, margin_metal_no);

fprintf('\n--- Nozzle Thermal Analysis Results with Coating ---\n');
fprintf('Max Coating Temp: %.1f K \n', ...
    max_T_coat);
fprintf('Max Metal Temp:   %.1f K | Limit: %.1f K | Margin: %.1f K\n', ...
    max_T_metal, T_cr_metal, margin_metal);

% Status Message logic
if max_T_metal < T_cr_metal
    fprintf('STATUS: [SUCCESS] All temperatures are within safety limits.\n');
else
    fprintf('STATUS: [DANGER] Critical temperature exceeded!\n');
    if max_T_metal >= T_cr_metal
        warning('Metal failure: T_metal > %.1f K', T_cr_metal);
    end
end
fprintf('---------------------------------------\n');

% --- 7. VISUALIZATION (Independent Colorbars & Throat Line) ---
figure('Color','w','Name','Nozzle Transient Heatmap','Position', [100 100 1500 450]);
mats = {T_inner, T_inter, T_outer};
titles = {'Gas/Coating Interface [K]', 'Coating/Metal Interface [K]', 'External Surface [K]'};

% Find throat axial position
[~, idx_throat] = min(R);
x_throat = X(idx_throat);

for j = 1:3
    subplot(1,3,j);
    imagesc(time, X, mats{j});
    set(gca, 'YDir', 'normal');
    hold on;

    % Add dashed line for throat position
    line([time(1), time(end)], [x_throat, x_throat], ...
        'Color', [0.8 0.8 0.8], 'LineStyle', '--', 'LineWidth', 1.5);

    % Individual colorbar and styling
    colorbar;
    colormap(hot);

    local_min = min(mats{j}(:));
    local_max = max(mats{j}(:));
    if local_max > local_min
        clim([local_min, local_max]);
    end

    title(titles{j});
    xlabel('Time [s]');
    ylabel('X [m]');

    % Optional: label the throat line in the first plot
    if j == 1
        text(time(round(end/10)), x_throat, ' Throat', ...
            'VerticalAlignment', 'bottom', 'Color', 'w', 'FontSize', 8);
    end
    hold off;
end

% Final Data Export
DataOut.T_inner = T_inner;
DataOut.T_inter = T_inter;
DataOut.T_outer = T_outer;
DataOut.T = T_inner(:,end);
end