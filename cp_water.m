function cp = cp_water(T_K)
% CP_WATER  Specific heat capacity of liquid water.
%
%   cp = cp_water(T_K)
%
%   Returns the specific heat capacity of pure liquid water [J/(kg·K)]
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
%   cp        specific heat capacity                            [J/(kg·K)]
%
%   -----------------------------------------------------------------------
%   Reference
%   -----------------------------------------------------------------------
%   Curve fit based on data from Incropera et al., Fundamentals of Heat
%   and Mass Transfer, 7th ed., Appendix A.

    Tc = T_K - 273.15;
    cp = 4217.6 - 3.453*Tc + 0.01379*Tc.^2;
end