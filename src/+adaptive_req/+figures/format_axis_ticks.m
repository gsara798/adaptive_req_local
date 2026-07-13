function format_axis_ticks(ax, nDecimals)
%FORMAT_AXIS_TICKS Format x and y axis tick labels.

if nargin < 2 || isempty(nDecimals)
    nDecimals = 1;
end

fmt = sprintf('%%.%df', nDecimals);

xtickformat(ax, fmt);

if numel(ax.YAxis) == 1
    ytickformat(ax, fmt);
else
    yyaxis(ax, 'left');
    ytickformat(ax, fmt);

    yyaxis(ax, 'right');
    ytickformat(ax, fmt);

    yyaxis(ax, 'left');
end

end