function Solution = solve_geometry_fzero2(V_sp, r_b, t_b, n)
% SOLVE_GEOMETRY_FZERO2  Solve for BATES grain internal geometry.
%
%   Solution = solve_geometry_fzero2(V_sp, r_b, t_b, n)
%
%   Determines the initial port radius r_p and chamber radius r_cc for a
%   cylindrical BATES grain that burns from the inner surface and both
%   ends.  The condition enforced is that the initial burning area A_b,i
%   equals the final burning area A_b,f (neutral burning profile).
%
%   -----------------------------------------------------------------------
%   Inputs
%   -----------------------------------------------------------------------
%   V_sp    total propellant volume                              [m³]
%   r_b     burning rate                                          [m/s]
%   t_b     burn time                                             [s]
%   n       number of grains (default = 1)                        [-]
%
%   -----------------------------------------------------------------------
%   Output
%   -----------------------------------------------------------------------
%   Solution   struct with fields:
%                .sol      = 1 if solution found, 0 otherwise
%                .r_p      initial port radius                    [m]
%                .r_cc     chamber inner radius                   [m]
%                .L_g      initial grain length                   [m]
%                .web      web thickness (= r_b * t_b)            [m]
%                .Ab_i     initial burning area                   [m²]
%                .Ab_f     final burning area                     [m²]

if nargin < 4, n = 1; end

web = r_b * t_b;

% residual: A_b,i - A_b,f
obj = @(rp) deal( ...
    2*pi*rp*(V_sp/(pi*((rp+web)^2-rp^2)*n))*n + 2*n*pi*((rp+web)^2-rp^2) - ...
    2*pi*(rp+web)*(V_sp/(pi*((rp+web)^2-rp^2)*n) - 2*n*web)*n );

% bracket search
rp_vec = linspace(1e-4, 0.5, 1000);
res    = arrayfun(obj, rp_vec);
idx    = find(diff(sign(res)) ~= 0, 1);

if isempty(idx), Solution.sol = 0; return; end

% refine with fzero
rp   = fzero(obj, [rp_vec(idx), rp_vec(idx+1)]);
rcc  = rp + web;
Lg   = V_sp / (pi*(rcc^2 - rp^2)*n);
Lg_f = Lg - 2*n*web;

Solution.sol  = 1;
Solution.r_p  = rp;
Solution.r_cc = rcc;
Solution.L_g  = Lg;
Solution.web  = web;
Solution.Ab_i = 2*pi*rp*Lg*n + 2*n*pi*(rcc^2 - rp^2);
Solution.Ab_f = 2*pi*rcc*Lg_f*n;

end