function ids = find_realisations(cfg)
% Numbers of the realisations of cfg.material with a tessellation
% data/<MAT><id>.tess, in ascending order.
f = dir(fullfile(cfg.data_dir, [cfg.material '*.tess']));
ids = zeros(1,0);
for i = 1:numel(f)
    t = regexp(f(i).name, ['^' cfg.material '(\d+)\.tess$'], 'tokens', 'once');
    if ~isempty(t), ids(end+1) = str2double(t{1}); end %#ok<AGROW>
end
ids = sort(ids);
end
