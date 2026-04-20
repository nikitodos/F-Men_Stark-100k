function meop_out = run_MC(k_scale, N_s, a_nom, n_nom, sa_rel, sn, ...
                           Geo, c_st, rho_p, m_p, Gas_in, lam, T_n, rb_n)
% RUN_MC  Helper function for Monte Carlo uncertainty calibration.
%
%   meop_out = run_MC(k_scale, N_s, a_nom, n_nom, sa_rel, sn, ...
%                     Geo, c_st, rho_p, m_p, Gas_in, lam, T_n, rb_n)
%
%   Runs N_s Monte Carlo samples of the Firing simulation, scaling the
%   raw uncertainties sa_rel and sn by the factor k_scale.  Returns the
%   99th percentile of the mean chamber pressure (MEOP) in bar.
%
%   This function is used inside a bisection loop to find the scale factor
%   k_unc that forces the simulated MEOP to match the design chamber
%   pressure (70 bar).  It is not intended for standalone use.
%
%   -----------------------------------------------------------------------
%   Inputs
%   -----------------------------------------------------------------------
%   k_scale   scale factor applied to raw uncertainties            [-]
%   N_s       number of Monte Carlo samples                        [-]
%   a_nom     nominal Vieille coefficient a                        [m/s/Pa^n]
%   n_nom     nominal Vieille exponent n                           [-]
%   sa_rel    raw relative uncertainty on a (from regression)      [-]
%   sn        raw absolute uncertainty on n (from regression)      [-]
%   Geo       geometry struct (see Firing.m)                       [-]
%   c_st      nominal characteristic velocity c*                   [m/s]
%   rho_p     propellant density                                   [kg/m³]
%   m_p       propellant mass                                      [kg]
%   Gas_in    gas properties struct (see Firing.m)                 [-]
%   lam       divergence loss factor                               [-]
%   T_n       nominal thrust                                       [N]
%   rb_n      nominal burning rate                                 [m/s]
%
%   -----------------------------------------------------------------------
%   Output
%   -----------------------------------------------------------------------
%   meop_out  99th percentile of mean chamber pressure distribution [bar]

    rng(42);
    pm = zeros(N_s,1);
    for ii = 1:N_s
        a_i = a_nom * exp(k_scale * sa_rel * randn());
        n_i = n_nom + k_scale * sn * randn();
        try
            En = Firing(Geo, a_i, n_i, c_st, rho_p, 1e-3, 1e-3, ...
                        m_p, Gas_in, lam, T_n, rb_n, false);
            pm(ii) = En.p_mean;
        catch
            pm(ii) = NaN;
        end
    end
    pm(isnan(pm)) = [];
    meop_out = prctile(pm, 99);
end