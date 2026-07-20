function [predictions, augmented_features, composition] = ...
    predict_frozen_models(feature_table, bundle)
%PREDICT_FROZEN_MODELS Apply the frozen REQ quantile model family.

arguments
    feature_table table
    bundle struct
end

if ~isfield(bundle, "MODELS")
    error( ...
        "adaptive_req:InvalidModelBundle", ...
        "Bundle is missing MODELS.");
end

models = bundle.MODELS;

composition = ...
    adaptive_req.inference.predict_composition( ...
        feature_table, ...
        models.composition);

augmented_features = feature_table;

augmented_features.predicted_patch_purity = ...
    composition.predicted_patch_purity;

augmented_features.p_mixed = ...
    composition.p_mixed;

augmented_features.p_strong_mixed = ...
    composition.p_strong_mixed;

model_names = ...
    string(models.model_names(:));

parts = cell(numel(model_names), 1);

for model_index = 1:numel(model_names)
    public_name = model_names(model_index);

    switch public_name
        case "q_spectrum_only"
            model_container = ...
                models.q.spectrum_only;

        case "q_spectrum_plus_composition"
            model_container = ...
                models.q.spectrum_plus_composition;

        otherwise
            error( ...
                "adaptive_req:UnknownFrozenModel", ...
                "Unknown frozen model: %s", ...
                public_name);
    end

    predictor_names = ...
        string(model_container.features(:));

    predictors = ...
        adaptive_req.inference.validate_predictor_table( ...
            augmented_features, ...
            predictor_names, ...
            public_name);

    q_raw = ...
        double(predict( ...
            model_container.model, ...
            predictors));

    q_raw = q_raw(:);

    if any(~isfinite(q_raw))
        error( ...
            "adaptive_req:InvalidQuantilePrediction", ...
            "%s produced non-finite quantiles.", ...
            public_name);
    end

    q_clamped = ...
        min(max(q_raw, 0), 1);

    number_of_rows = ...
        height(feature_table);

    result = table();

    result.source_row = ...
        (1:number_of_rows).';

    result.model_name = ...
        repmat(public_name, number_of_rows, 1);

    result.q_raw = q_raw;
    result.q_pred = q_clamped;

    result.predicted_patch_purity = ...
        composition.predicted_patch_purity;

    result.p_mixed = ...
        composition.p_mixed;

    result.p_strong_mixed = ...
        composition.p_strong_mixed;

    coordinate_names = [
        "map_iz"
        "map_ix"
        "cx"
        "cz"
        "x_center_m"
        "z_center_m"
    ];

    for coordinate_name = coordinate_names.'
        if ismember( ...
                coordinate_name, ...
                string(feature_table.Properties.VariableNames))
            result.(coordinate_name) = ...
                feature_table.(coordinate_name);
        end
    end

    parts{model_index} = result;
end

predictions = ...
    vertcat(parts{:});

end
