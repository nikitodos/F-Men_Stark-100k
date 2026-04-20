function mu = viscosity_water(T_K)
% VISCOSITY_WATER  Dynamic viscosity of liquid water.
%
%   mu = viscosity_water(T_K)
%
%   Returns the dynamic viscosity of pure liquid water [Pa·s] as a function
%   of temperature, using the correlation μ = 2.414e-5 * 10^(247.8/(T[°C]+133.15)).
%   Valid for 0–100 °C.
%
%   -----------------------------------------------------------------------
%   Input
%   -----------------------------------------------------------------------
%   T_K       temperature                                       [K]
%
%   -----------------------------------------------------------------------
%   Output
%   -----------------------------------------------------------------------
%   mu        dynamic viscosity                                 [Pa·s]

    Tc = T_K - 273.15;
    mu = 2.414e-5 * 10.^(247.8./(Tc + 133.15));
end