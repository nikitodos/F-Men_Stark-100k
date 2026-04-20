function P = P_sat_water(T_K)
% P_SAT_WATER  Saturation pressure of water (Antoine equation).
%
%   P = P_sat_water(T_K)
%
%   Returns the saturation vapour pressure of pure water [Pa] as a function
%   of temperature, using the Antoine equation with coefficients valid for
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
%   P         saturation pressure                               [Pa]

    A = 8.07131; B = 1730.63; C = 233.426;
    Tc = T_K - 273.15;
    P_mmHg = 10^(A - B/(Tc+C));
    P = P_mmHg * 133.322;
end