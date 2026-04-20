function hg = bartz_curv(Dt, mu0, cp, Pr_g, G_t, rc, AR, sigma_B)
% BARTZ_CURV  Bartz heat transfer coefficient with throat curvature correction.
%
%   hg = bartz_curv(Dt, mu0, cp, Pr_g, G_t, rc, AR, sigma_B)
%
%   Computes the gas-side convective heat transfer coefficient hg [W/(m²·K)]
%   according to the classical Bartz (1957) correlation, including the
%   throat curvature term (Dt/rc)^0.1 and the sigma correction factor.
%
%   -----------------------------------------------------------------------
%   Inputs
%   -----------------------------------------------------------------------
%   Dt        throat diameter                                   [m]
%   mu0       gas dynamic viscosity at chamber temperature      [Pa·s]
%   cp        gas specific heat at constant pressure            [J/(kg·K)]
%   Pr_g      gas Prandtl number                                [-]
%   G_t       mass flux at throat                               [kg/(m²·s)]
%   rc        throat radius of curvature (upstream)             [m]
%   AR        local area ratio A/A_t                            [-]
%   sigma_B   Bartz sigma correction factor (from sigma_bartz)  [-]
%
%   -----------------------------------------------------------------------
%   Output
%   -----------------------------------------------------------------------
%   hg        gas-side convective heat transfer coefficient     [W/(m²·K)]
%
%   -----------------------------------------------------------------------
%   Reference
%   -----------------------------------------------------------------------
%   Bartz D.R. (1957), "A Simple Equation for Rapid Estimation of
%   Rocket Nozzle Convective Heat Transfer Coefficients", Jet Propulsion
%   27(7):558-566.
%
%   See also: sigma_bartz, thermal_chain.

    hg = (0.026./Dt.^0.2) * (mu0.^0.2 * cp ./ Pr_g.^0.6) * G_t.^0.8 ...
         * (Dt./rc).^0.1 .* (1./AR).^0.9 .* sigma_B;
end