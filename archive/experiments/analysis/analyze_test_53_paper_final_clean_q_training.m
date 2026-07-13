%% analyze_test_53_paper_final_clean_q_training.m
% Test 53: paper-final clean q/SWS training orchestration.
%
% This script intentionally reuses the mature Test 38 training engine, but
% pins the paper-facing design to common frequencies, M values, and output
% folders. It trains clean spectral q models only. No confidence detector,
% prior frozen q model, true patch purity, or readout/noise variable is used
% as an operational predictor.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST53_MODE                  = quick | medium | full
%   ADAPTIVE_REQ_TEST53_VALIDATE_ONLY         = true | false
%   ADAPTIVE_REQ_TEST53_SAVE_ALL_MAPS         = true | false
%   ADAPTIVE_REQ_TEST53_USE_PARFOR            = true | false
%   ADAPTIVE_REQ_TEST53_USE_PARALLEL_TRAINING = true | false
%   ADAPTIVE_REQ_TEST53_MAX_PATCHES_PER_CONDITION
%   ADAPTIVE_REQ_TEST53_MAX_MODEL_TRAIN_ROWS
%   ADAPTIVE_REQ_TEST53_MAX_MODEL_EVAL_ROWS
%
% Defaults for paper v1:
%   frequencies = [200 300 400 500 600] Hz
%   M           = [2 3]
%   dx = dz     = 0.2 mm
%   clean simulations only

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));

mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST53_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","medium","full"]), ...
    'ADAPTIVE_REQ_TEST53_MODE must be quick, medium, or full.');

fprintf('\nTest 53: paper-final clean q/SWS training\n');
fprintf('This delegates to Test 38 with paper-final controlled overrides.\n');
fprintf('Frequencies: [200 300 400 500 600] Hz | M: [2 3] | dx=0.2 mm\n');
fprintf('No realistic readout/SNR is included in this training set.\n\n');

forward_env('ADAPTIVE_REQ_TEST53_MODE', 'ADAPTIVE_REQ_TEST38_MODE', char(mode));
forward_env('ADAPTIVE_REQ_TEST53_VALIDATE_ONLY', 'ADAPTIVE_REQ_TEST38_VALIDATE_ONLY', '');
forward_env('ADAPTIVE_REQ_TEST53_SAVE_ALL_MAPS', 'ADAPTIVE_REQ_TEST38_SAVE_ALL_MAPS', '');
forward_env('ADAPTIVE_REQ_TEST53_USE_PARFOR', 'ADAPTIVE_REQ_TEST38_USE_PARFOR', '');
forward_env('ADAPTIVE_REQ_TEST53_USE_PARALLEL_TRAINING', ...
    'ADAPTIVE_REQ_TEST38_USE_PARALLEL_TRAINING', '');
forward_env('ADAPTIVE_REQ_TEST53_MAX_PATCHES_PER_CONDITION', ...
    'ADAPTIVE_REQ_TEST38_MAX_PATCHES_PER_CONDITION', '');
forward_env('ADAPTIVE_REQ_TEST53_MAX_MODEL_TRAIN_ROWS', ...
    'ADAPTIVE_REQ_TEST38_MAX_MODEL_TRAIN_ROWS', '');
forward_env('ADAPTIVE_REQ_TEST53_MAX_MODEL_EVAL_ROWS', ...
    'ADAPTIVE_REQ_TEST38_MAX_MODEL_EVAL_ROWS', '');

setenv('ADAPTIVE_REQ_TEST38_OUTPUT_NAME', 'test_53_paper_final_clean_q_training');
setenv('ADAPTIVE_REQ_TEST38_FREQUENCIES_HZ', '200 300 400 500 600');
setenv('ADAPTIVE_REQ_TEST38_M_LIST', '2 3');
setenv('ADAPTIVE_REQ_TEST38_DX_LIST_MM', '0.2');

run(fullfile(root_dir, 'experiments', 'analysis', ...
    'analyze_test_38_velocity_field_diverse_q_training.m'));

function forward_env(src_name, dst_name, default_value)
val = getenv(src_name);
if isempty(val)
    val = default_value;
end
if ~isempty(val)
    setenv(dst_name, val);
end
end
