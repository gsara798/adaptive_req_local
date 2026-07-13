function T_out = add_effective_window_metrics(T_in, varargin)
%ADD_EFFECTIVE_WINDOW_METRICS Add effective window-size diagnostics.
%
% This function adds diagnostic quantities related to the number of true
% wavelengths contained in the analysis window.
%
% These quantities use cs_bg and are meant for simulation diagnostics, not
% necessarily for the final experimental model.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.add_effective_window_metrics';

addRequired(p, 'T_in', @istable);

addParameter(p, 'MVar', 'REQ_M', @(x) ischar(x) || isstring(x));
addParameter(p, 'F0Var', 'SIM_f0', @(x) ischar(x) || isstring(x));
addParameter(p, 'CsTrueVar', 'SIM_cs_bg', @(x) ischar(x) || isstring(x));
addParameter(p, 'CsGuess', 3.0, @(x) isnumeric(x) && isscalar(x) && x > 0);

parse(p, T_in, varargin{:});

M_var = char(p.Results.MVar);
f0_var = char(p.Results.F0Var);
cs_true_var = char(p.Results.CsTrueVar);
cs_guess = p.Results.CsGuess;

required_vars = string({M_var, f0_var, cs_true_var});

for i = 1:numel(required_vars)
    if ~ismember(required_vars(i), string(T_in.Properties.VariableNames))
        error('Required variable not found: %s', required_vars(i));
    end
end

T_out = T_in;

M = T_out.(M_var);
f0 = T_out.(f0_var);
cs_true = T_out.(cs_true_var);

T_out.cs_guess_diag = repmat(cs_guess, height(T_out), 1);

T_out.lambda_true_diag = cs_true ./ f0;
T_out.lambda_guess_diag = cs_guess ./ f0;

T_out.window_length_guess_diag = M .* T_out.lambda_guess_diag;

T_out.M_eff_true_diag = ...
    T_out.window_length_guess_diag ./ T_out.lambda_true_diag;

T_out.M_eff_true_diag_alt = M .* cs_guess ./ cs_true;

T_out.cs_ratio_guess_to_true_diag = cs_guess ./ cs_true;

end