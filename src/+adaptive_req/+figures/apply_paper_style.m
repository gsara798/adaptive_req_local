function apply_paper_style(fig, ax, fontName, fontSize)

    if nargin < 3 || isempty(fontName)
        fontName = 'Times New Roman';
    end
    if nargin < 4 || isempty(fontSize)
        fontSize = 15;  %era 15
    end

    if nargin < 2 || isempty(ax)
        ax = findall(fig, 'Type', 'axes');
    end

    ax = ax(isgraphics(ax));
    ax = ax(~arrayfun(@(h) strcmp(get(h, 'Tag'), 'Colorbar'), ax));

    set(findall(fig, '-property', 'FontName'), 'FontName', fontName);
    set(findall(fig, '-property', 'FontSize'), 'FontSize', fontSize);

    set(ax, ...
        'Box', 'on', ...
        'TickDir', 'in', ...
        'LineWidth', 1, ...
        'FontName', fontName, ...
        'FontSize', fontSize);

    lgd = findall(fig, 'Type', 'Legend');
    if ~isempty(lgd)
        set(lgd, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'on');
    end
end