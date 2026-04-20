function s = ternary_str(condition, s_true, s_false)
% TERNARY_STR  Return s_true if condition is true, s_false otherwise.
%
%   s = ternary_str(condition, s_true, s_false)
%
%   Convenience one-liner for conditional string selection in fprintf calls.
%   Avoids multi-line if/else blocks when building formatted output.
%
%   Example:
%     fprintf('Flow regime: %s\n', ternary_str(Re > 4000, 'turbulent', 'laminar'));
%
%   Inputs
%   ------
%   condition   scalar logical (or numeric; treated as logical)
%   s_true      string returned when condition is true
%   s_false     string returned when condition is false
%
%   Output
%   ------
%   s           char row vector

if condition
    s = s_true;
else
    s = s_false;
end
end