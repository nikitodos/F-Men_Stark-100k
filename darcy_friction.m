function f = darcy_friction(Re, roughness, D_h)
% DARCY_FRICTION  Darcy-Weisbach friction factor for pipe flow.
%
%   f = darcy_friction(Re, roughness, D_h)
%
%   Computes the Darcy friction factor f [-] for fully developed internal
%   flow.  For Re < 2300, laminar f = 64/Re is returned.  For turbulent
%   flow (Re ≥ 2300), the Colebrook-White equation is solved iteratively
%   if the relative roughness roughness/D_h > 0; otherwise the smooth
%   Blasius correlation f = 0.316·Re^(-0.25) is used.
%
%   -----------------------------------------------------------------------
%   Inputs
%   -----------------------------------------------------------------------
%   Re          Reynolds number based on hydraulic diameter        [-]
%   roughness   absolute pipe roughness                            [m]
%   D_h         hydraulic diameter                                 [m]
%
%   -----------------------------------------------------------------------
%   Output
%   -----------------------------------------------------------------------
%   f           Darcy friction factor                              [-]

    if Re < 2300; f = 64/Re; return; end
    if roughness == 0
        f = 0.316 * Re^(-0.25); return;
    end
    eps_D = roughness / D_h;
    f = 0.316 * Re^(-0.25);
    for iter = 1:50
        f_new = 1/(-2*log10(eps_D/3.7 + 2.51/(Re*sqrt(f))))^2;
        if abs(f_new-f) < 1e-10; f = f_new; return; end
        f = f_new;
    end
end