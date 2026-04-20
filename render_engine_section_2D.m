function [X_tot, R_tot] = render_engine_section_2D(E, plot_flag)
% RENDER_ENGINE_SECTION_2D  2‑D axisymmetric cross‑section of the SRM.
%
%   [X_tot, R_tot] = render_engine_section_2D(E, plot_flag)
%
%   Builds the internal and external profiles of the solid rocket motor,
%   including the combustion chamber (with insulation) and the convergent‑
%   divergent nozzle (with internal coating and structural wall).  When
%   plot_flag = 1, a figure is generated showing the propellant grains,
%   insulation, coating, and metal structure.
%
%   -----------------------------------------------------------------------
%   Inputs
%   -----------------------------------------------------------------------
%   E           struct containing geometry fields (output of SolidUnitDesign)
%               required fields:
%                 .Geom.r_cc      chamber inner radius            [m]
%                 .Geom.L_g       grain length                    [m]
%                 .Geom.N         number of grains                [-]
%                 .Geom.r_t       throat radius                   [m]
%                 .Geom.r_e       exit radius                     [m]
%                 .Geom.L_conv    convergent section length       [m]
%                 .Geom.L_div     divergent section length        [m]
%                 .Geom.L_cc      chamber length                  [m]
%               optional fields (default values used if missing):
%                 .Geom.t_cc      chamber wall thickness          [m]
%                 .Geom.t_nozzle  nozzle wall thickness           [m]
%                 .Geom.t_iso     insulation thickness (CC)       [m]
%                 .Geom.t_coat    internal coating thickness (NZ) [m]
%                 .Geom.r_p       initial port radius (for grain) [m]
%   plot_flag   logical: 1 = create figure, 0 = silent            [-]
%
%   -----------------------------------------------------------------------
%   Outputs
%   -----------------------------------------------------------------------
%   X_tot       [1 × N] axial coordinates of the internal profile  [m]
%   R_tot       [1 × N] radial coordinates of the internal profile [m]
%
%   -----------------------------------------------------------------------
%   Notes
%   -----------------------------------------------------------------------
%   - The nozzle profile uses circular arcs at the throat (R_u upstream,
%     R_d downstream) and straight conical walls.
%   - The coating and structural layers are obtained by a normal offset
%     of the internal contour.

% 2D axisymmetric SRM section
% Layers: Propellant -> Insulation (CC) -> Coating (Nozzle) -> Structure

%% --- Geometry ---
R_cc   = E.Geom.r_cc;
L_g    = E.Geom.L_g;
N      = E.Geom.N;
R_t    = E.Geom.r_t;
R_e    = E.Geom.r_e;
L_conv = E.Geom.L_conv;
L_div  = E.Geom.L_div;
L_cc   = E.Geom.L_cc;

% Thicknesses (default values if missing)
t_cc   = 0.02 * R_cc; if isfield(E.Geom, 't_cc'),     t_cc = E.Geom.t_cc; end
t_nz   = 0.02 * R_cc; if isfield(E.Geom, 't_nozzle'), t_nz = E.Geom.t_nozzle; end
t_iso  = 0;           if isfield(E.Geom, 't_iso'),    t_iso = E.Geom.t_iso; end
t_coat = 0;           if isfield(E.Geom, 't_coat'),   t_coat = E.Geom.t_coat; end

% Angles and radii
beta_calc  = atand((R_cc - R_t) / L_conv);
alpha_calc = atand((R_e - R_t) / L_div);
R_curv_cc  = 0.3 * (R_cc - R_t);
R_u        = 1.5 * R_t;
R_d        = 0.5 * R_t;

individual_gap = (L_cc - N * L_g) / (N + 1);

%% --- Internal profile ---
% Combustion chamber
X_CC_int = [0, L_cc];
R_CC_int = [R_cc, R_cc];

% Nozzle profile
theta_cc   = linspace(0, beta_calc, 30);
x_round_cc = L_cc + R_curv_cc * sind(theta_cc);
r_round_cc = (R_cc - R_curv_cc) + R_curv_cc * cosd(theta_cc);

x_t_line   = L_cc + L_conv;
x_tang_u   = x_t_line - R_u * sind(beta_calc);
r_tang_u   = R_t + R_u * (1 - cosd(beta_calc));

theta_u    = linspace(beta_calc, 0, 30);
x_round_u  = x_t_line - R_u * sind(theta_u);
r_round_u  = R_t + R_u * (1 - cosd(theta_u));

theta_d    = linspace(0, alpha_calc, 30);
x_round_d  = x_t_line + R_d * sind(theta_d);
r_round_d  = R_t + R_d * (1 - cosd(theta_d));

x_end_nozzle = L_cc + L_conv + L_div;

X_NZ_int = [x_round_cc, x_tang_u, x_round_u, x_round_d, x_end_nozzle];
R_NZ_int = [r_round_cc, r_tang_u, r_round_u, r_round_d, R_e];

%% --- Combustion chamber layers ---
X_CC_ins = X_CC_int;
R_CC_ins = R_CC_int + t_iso;

X_CC_ext = X_CC_int;
R_CC_ext = R_CC_ins + t_cc;

%% --- Nozzle layers (normal offset) ---
dx = gradient(X_NZ_int);
dr = gradient(R_NZ_int);
mag = sqrt(dx.^2 + dr.^2);

nx = -dr ./ mag;
nr =  dx ./ mag;

% Coating (inner layer)
X_NZ_coat = X_NZ_int + t_coat * nx;
R_NZ_coat = R_NZ_int + t_coat * nr;

% Structure (outer layer)
X_NZ_strut = X_NZ_coat + t_nz * nx;
R_NZ_strut = R_NZ_coat + t_nz * nr;

%% --- Output profile ---
X_tot = [X_CC_int, X_NZ_int];
R_tot = [R_CC_int, R_NZ_int];

%% --- Plot ---
if plot_flag == 1
    figure('Color', 'w', 'Name', 'SRM 2D Section'); hold on;

    % --- Propellant grains ---
    if isfield(E.Geom, 'r_p')
        R_p = E.Geom.r_p;
        for k = 1:N
            x_start = k * individual_gap + (k-1) * L_g;
            x_end   = x_start + L_g;

            patch([x_start x_end x_end x_start], ...
                  [R_p R_p R_cc R_cc], ...
                  [0.9 0.9 0.9], 'EdgeColor', [0.6 0.6 0.6]);
        end
    end

    % --- Combustion chamber ---
    if t_iso > 0
        patch([X_CC_int, fliplr(X_CC_ins)], ...
              [R_CC_int, fliplr(R_CC_ins)], ...
              [0.95 0.9 0.6], 'EdgeColor', 'none');
    end

    patch([X_CC_ins, fliplr(X_CC_ext)], ...
          [R_CC_ins, fliplr(R_CC_ext)], ...
          [0.3 0.3 0.3], 'EdgeColor', 'k');

    % --- Nozzle coating (inner) ---
    if t_coat > 0
        patch([X_NZ_int, fliplr(X_NZ_coat)], ...
              [R_NZ_int, fliplr(R_NZ_coat)], ...
              [0.95 0.9 0.6], 'EdgeColor', 'none');
    end

    % --- Nozzle structure (outer) ---
    patch([X_NZ_coat, fliplr(X_NZ_strut)], ...
          [R_NZ_coat, fliplr(R_NZ_strut)], ...
          [0.3 0.3 0.3], 'EdgeColor', 'k');

    % --- Formatting ---
    plot([0, x_end_nozzle], [0, 0], '-.', 'Color', [0.5 0.5 0.5]);
    xline(x_t_line, ':', 'Throat');

    axis equal; grid on; box on;
    xlabel('Axial Position x [m]');
    ylabel('Radius r [m]');
    title('SRM Section with Internal Nozzle Coating');

    xlim([0, x_end_nozzle * 1.05]);
    ylim([0, (R_cc + t_iso + t_cc) * 1.3]);

    hold off;
end

end