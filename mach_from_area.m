function M = mach_from_area(AR, gamma, supersonic)
% MACH_FROM_AREA  Mach number from area ratio for isentropic flow.
%
%   M = mach_from_area(AR, gamma, supersonic)
%
%   Solves the area‑Mach relation A/A* = (1/M)*((2/(γ+1))*(1+(γ-1)/2*M²))^((γ+1)/(2(γ-1)))
%   for the Mach number M.  The user must specify whether the solution lies
%   in the subsonic (M < 1) or supersonic (M > 1) branch.
%
%   -----------------------------------------------------------------------
%   Inputs
%   -----------------------------------------------------------------------
%   AR          area ratio A/A* (≥ 1)                          [-]
%   gamma       specific heat ratio                            [-]
%   supersonic  logical: true  → supersonic branch (M > 1)
%                       false → subsonic branch   (M < 1)
%
%   -----------------------------------------------------------------------
%   Output
%   -----------------------------------------------------------------------
%   M           Mach number                                    [-]
%
%   -----------------------------------------------------------------------
%   Algorithm
%   -----------------------------------------------------------------------
%   Bisection method with 300 iterations and tolerance 1e-10.

    if supersonic; lo = 1.0; hi = 50.0;
    else;          lo = 1e-4; hi = 1.0 - 1e-10; end
    f = @(M) (1/M)*((2/(gamma+1))*(1+(gamma-1)/2*M^2))^((gamma+1)/(2*(gamma-1))) - AR;
    % FIX: Replaced '~' with 'iter'
    for iter = 1:300
        mid = 0.5*(lo+hi);
        if f(mid)*f(lo) < 0; hi = mid; else; lo = mid; end
        if hi-lo < 1e-10; break; end
    end
    M = 0.5*(lo+hi);
end
