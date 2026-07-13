function [map, xc, zc] = table_to_map_grid(T, value_var)
%TABLE_TO_MAP_GRID Convert patch-center table values to image grid.

value_var = char(value_var);
xc = sort(unique(T.x_center_m(:))).';
zc = sort(unique(T.z_center_m(:)));

map = nan(numel(zc), numel(xc));
[~, ix] = ismember(T.x_center_m, xc);
[~, iz] = ismember(T.z_center_m, zc);

for i = 1:height(T)
    if ix(i) > 0 && iz(i) > 0
        map(iz(i), ix(i)) = T.(value_var)(i);
    end
end

end
