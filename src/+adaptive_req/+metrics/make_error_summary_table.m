function T_error = make_error_summary_table(method_names, q_hat_cell, cs_hat_cell, ape_cell, q_true, cs_true)

nM = numel(method_names);

rmse_q = nan(nM,1);
mae_q = nan(nM,1);
rho_q = nan(nM,1);

rmse_cs = nan(nM,1);
mae_cs = nan(nM,1);

mean_ape = nan(nM,1);
median_ape = nan(nM,1);
p75_ape = nan(nM,1);
p90_ape = nan(nM,1);
max_ape = nan(nM,1);

for m = 1:nM

    qhat = q_hat_cell{m};
    cshat = cs_hat_cell{m};
    ape = ape_cell{m};

    valid_q = isfinite(qhat) & isfinite(q_true);
    valid_cs = isfinite(cshat) & isfinite(cs_true);
    valid_ape = isfinite(ape);

    rmse_q(m) = sqrt(mean((qhat(valid_q) - q_true(valid_q)).^2, 'omitnan'));
    mae_q(m) = mean(abs(qhat(valid_q) - q_true(valid_q)), 'omitnan');

    if nnz(valid_q) >= 3
        rho_q(m) = corr(qhat(valid_q), q_true(valid_q), 'Type', 'Spearman');
    end

    rmse_cs(m) = sqrt(mean((cshat(valid_cs) - cs_true(valid_cs)).^2, 'omitnan'));
    mae_cs(m) = mean(abs(cshat(valid_cs) - cs_true(valid_cs)), 'omitnan');

    mean_ape(m) = mean(ape(valid_ape), 'omitnan');
    median_ape(m) = median(ape(valid_ape), 'omitnan');
    p75_ape(m) = prctile(ape(valid_ape), 75);
    p90_ape(m) = prctile(ape(valid_ape), 90);
    max_ape(m) = max(ape(valid_ape));
end

T_error = table( ...
    string(method_names), ...
    rmse_q, mae_q, rho_q, ...
    rmse_cs, mae_cs, ...
    mean_ape, median_ape, p75_ape, p90_ape, max_ape, ...
    'VariableNames', { ...
    'method', ...
    'rmse_q', 'mae_q', 'rho_q', ...
    'rmse_cs_mps', 'mae_cs_mps', ...
    'mean_APE_percent', 'median_APE_percent', ...
    'p75_APE_percent', 'p90_APE_percent', 'max_APE_percent'});

end