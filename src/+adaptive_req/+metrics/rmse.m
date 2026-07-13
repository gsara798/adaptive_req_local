function val = rmse(yhat, y)
    val = sqrt(mean((yhat - y).^2, 'omitnan'));
end