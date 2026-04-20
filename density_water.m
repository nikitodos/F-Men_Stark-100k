function rho = density_water(T_K)
% DENSITY_WATER  Density of liquid water.
%
%   rho = density_water(T_K)
%
%   Returns the density of pure liquid water [kg/m³] as a function of
%   temperature, using a cubic polynomial fit valid for 0–100 °C.
%
%   -----------------------------------------------------------------------
%   Input
%   -----------------------------------------------------------------------
%   T_K       temperature                                       [K]
%
%   -----------------------------------------------------------------------
%   Output
%   -----------------------------------------------------------------------
%   rho       density                                           [kg/m³]

    Tc = T_K - 273.15;
    rho = 999.84 + 0.0673*Tc - 0.00894*Tc.^2 + 2.57e-5*Tc.^3;
end
