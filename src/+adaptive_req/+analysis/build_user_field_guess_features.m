function T_out = build_user_field_guess_features(T, varargin)
%BUILD_USER_FIELD_GUESS_FEATURES Add synthetic user field guess priors.
%
% This helper can either add a single "unknown" guess to the input table or
% expand the table four-fold with all supported synthetic guesses.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.build_user_field_guess_features';
addRequired(p, 'T', @istable);
addParameter(p, 'ExpandGuesses', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'DefaultGuess', "unknown", @(x) ischar(x) || isstring(x));
parse(p, T, varargin{:});

required = ["q_theory_dir2D", "q_theory_diffuse2D", "q_theory_projected3D", ...
    "q_theory_mean_dir2D_projected3D", "q_theory_mean_all"];
missing = required(~ismember(required, string(T.Properties.VariableNames)));
assert(isempty(missing), 'Missing theory-q candidate variables: %s', ...
    strjoin(missing, ', '));

if logical(p.Results.ExpandGuesses)
    guesses = ["directional_like"; "partially_diffuse"; "diffuse_like"; "unknown"];
    T_out = table();
    for i = 1:numel(guesses)
        Ti = T;
        Ti.user_field_guess = categorical(repmat(guesses(i), height(Ti), 1), guesses);
        Ti.q_user_guess_prior = prior_from_guess(Ti, guesses(i));
        T_out = adaptive_req.analysis.Test12Analysis.concatTables(T_out, Ti);
    end
else
    guess = string(p.Results.DefaultGuess);
    guesses = ["directional_like"; "partially_diffuse"; "diffuse_like"; "unknown"];
    T_out = T;
    T_out.user_field_guess = categorical(repmat(guess, height(T), 1), guesses);
    T_out.q_user_guess_prior = prior_from_guess(T_out, guess);
end

end

function q = prior_from_guess(T, guess)

switch string(guess)
    case "directional_like"
        q = T.q_theory_dir2D;
    case "partially_diffuse"
        q = T.q_theory_mean_dir2D_projected3D;
    case "diffuse_like"
        q = T.q_theory_projected3D;
    otherwise
        q = T.q_theory_mean_all;
end

end
