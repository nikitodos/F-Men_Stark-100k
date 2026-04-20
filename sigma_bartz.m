function sigma_B = sigma_bartz(T_wg, Tc, gamma, M, omega)
% SIGMA_BARTZ  Bartz sigma correction factor for property variation.
%
%   sigma_B = sigma_bartz(T_wg, Tc, gamma, M, omega)
%
%   Computes the sigma correction factor that accounts for the variation
%   of gas properties across the boundary layer in the Bartz heat transfer
%   correlation.  The exponent ω is the viscosity temperature exponent
%   (μ ∝ T^ω).  Standard value for combustion gases is ω ≈ 0.6.
%
%   -----------------------------------------------------------------------
%   Inputs
%   -----------------------------------------------------------------------
%   T_wg    wall temperature (gas‑side)                           [K]
%   Tc      chamber stagnation temperature                        [K]
%   gamma   specific heat ratio                                   [-]
%   M       local Mach number                                     [-]
%   omega   viscosity temperature exponent                        [-]
%
%   -----------------------------------------------------------------------
%   Output
%   -----------------------------------------------------------------------
%   sigma_B   sigma correction factor (dimensionless)             [-]
%
%   -----------------------------------------------------------------------
%   Reference
%   -----------------------------------------------------------------------
%   Bartz D.R. (1957), Jet Propulsion 27(7):558-566.
%   Huzel D.K. & Huang D.H. (1992), Modern Engineering for Design of
%   Liquid‑Propellant Rocket Engines, AIAA, Eq. (4-19).

    fac     = 1 + (gamma-1)/2 * M.^2;
    ratio   = T_wg / Tc;
    sigma_B = 1 ./ ((0.5*ratio.*fac + 0.5).^(0.8 - omega*0.2) .* fac.^0.12);
end