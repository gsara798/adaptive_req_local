function format_colorbar_ticks(cb, nDecimals)
%FORMAT_COLORBAR_TICKS Format colorbar tick labels.

if nargin < 2 || isempty(nDecimals)
    nDecimals = 1;
end

fmt = sprintf('%%.%df', nDecimals);
cb.TickLabels = compose(fmt, cb.Ticks);

end