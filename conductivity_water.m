function k = conductivity_water(T_K)
% CONDUCTIVITY_WATER  Thermal conductivity of liquid water.
%
%   k = conductivity_water(T_K)
%
%   Returns the thermal conductivity of pure liquid water [W/(m·K)]
%   as a function of temperature, using a quadratic fit valid for
%   0–100 °C (273–373 K).
%
%   -----------------------------------------------------------------------
%   Input
%   -----------------------------------------------------------------------
%   T_K       temperature                                       [K]
%
%   -----------------------------------------------------------------------
%   Output
%   -----------------------------------------------------------------------
%   k         thermal conductivity                              [W/(m·K)]
%
%   -----------------------------------------------------------------------
%   Reference
%   -----------------------------------------------------------------------
%   Curve fit based on data from Incropera et al., Fundamentals of Heat
%   and Mass Transfer, 7th ed., Appendix A.

    Tc = T_K - 273.15;
    k = 0.5636 + 1.946e-3*Tc - 8.151e-6*Tc.^2;
end