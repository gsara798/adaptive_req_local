function sws_m_s = quantiles_to_sws( ...
    req_mappings, quantiles, frequency_hz)
%QUANTILES_TO_SWS Convert predicted q values using local REQ mappings.

arguments
    req_mappings
    quantiles (:,1) double
    frequency_hz double
end

req_mappings = req_mappings(:);
quantiles = double(quantiles(:));

if numel(req_mappings) ~= numel(quantiles)
    error( ...
        "adaptive_req:MappingQuantileSizeMismatch", ...
        "REQ mappings and quantiles must have the same length.");
end

if isscalar(frequency_hz)
    frequency_hz = ...
        repmat(double(frequency_hz), numel(quantiles), 1);
else
    frequency_hz = ...
        double(frequency_hz(:));
end

if numel(frequency_hz) ~= numel(quantiles)
    error( ...
        "adaptive_req:FrequencyQuantileSizeMismatch", ...
        "Frequency and quantile vectors must have the same length.");
end

sws_m_s = ...
    nan(numel(quantiles), 1);

for index = 1:numel(quantiles)
    mapping = req_mappings{index};

    sws_m_s(index) = ...
        adaptive_req.quantile.quantile_to_cs( ...
            mapping, ...
            quantiles(index), ...
            frequency_hz(index));
end

end
