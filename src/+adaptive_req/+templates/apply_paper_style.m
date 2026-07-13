function apply_paper_style(fig, ax, font_name, font_size)
%APPLY_PAPER_STYLE Apply clean paper-like style to figure axes.

if nargin < 3 || isempty(font_name)
    font_name = 'Times New Roman';
end

if nargin < 4 || isempty(font_size)
    font_size = 18;
end

if nargin < 2 || isempty(ax)
    ax = findall(fig, 'Type', 'axes');
end

ax = ax(isgraphics(ax));
ax = ax(~arrayfun(@(h) strcmp(get(h, 'Tag'), 'Colorbar'), ax));

set(findall(fig, '-property', 'FontName'), 'FontName', font_name);
set(findall(fig, '-property', 'FontSize'), 'FontSize', font_size);

set(ax, ...
    'Box', 'on', ...
    'TickDir', 'in', ...
    'LineWidth', 1.2, ...
    'FontName', font_name, ...
    'FontSize', font_size);

lgd = findall(fig, 'Type', 'Legend');
if ~isempty(lgd)
    set(lgd, 'FontName', font_name, 'FontSize', font_size, 'Box', 'on');
end

end
