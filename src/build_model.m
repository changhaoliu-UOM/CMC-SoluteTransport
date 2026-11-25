function M = build_model(cfg, rid)
% Builds the model of realisation rid of the material in cfg (config_material.m):
% the physical complex M from data/<MAT><rid>.tess, its Forman subdivision K,
% the pore assignment, the conductances of the 1-cells of K and the pore
% storage of the vertices of K.
%
% M.fcmplx, M.vcmplx  complexes K and M
% M.fBO_1             boundary operator of the 1-cells of K (n1 x n0)
% M.fvol0             dual volume of every vertex of K
% M.vol_cells, M.h_face   voluminous pores (3-cells of M), sheet apertures (2-cells of M)
% M.G, M.Dop          diagonal hydraulic and diffusive conductances (n1 x n1)
% M.node_active       vertices of K that are monitored
% M.S_pore            pore volume of every vertex of K (um^3)
% M.pores             porosities and counts of the pore assignment
tess  = fullfile(cfg.data_dir, sprintf('%s%d.tess', cfg.material, rid));
f_vol = fullfile(cfg.data_dir, [cfg.pore_data '_voluminous_pores.csv']);
f_sh  = fullfile(cfg.data_dir, [cfg.pore_data '_sheet_pores_cdf.csv']);
for f = {tess, f_vol, f_sh}
    if ~isfile(f{1}), error('build_model: missing input file %s', f{1}); end
end

[vcmplx, fcmplx, fBO_1, fs1] = build_complexes(tess);

V_vol   = readmatrix(f_vol);                       % volumes of the voluminous pores (um^3)
V_vol   = repmat(V_vol(:,1), cfg.void_copies, 1);
V_sheet = readmatrix(f_sh);                        % [volume (um^3), cumulative probability]
[vol_cells, h_face, pores] = assign_pores(vcmplx, fcmplx, cfg.vf, V_vol, V_sheet(:,1:2));

M = struct('fcmplx', fcmplx, 'vcmplx', vcmplx, 'fBO_1', fBO_1, 'fvol0', full(diag(fs1)), ...
           'vol_cells', vol_cells, 'h_face', h_face, 'pores', pores);

[gH, gD, node_active] = pore_conductances(M, cfg.mu, cfg.Dm);
n1 = numel(gH);
M.G   = spdiags(gH, 0, n1, n1);
M.Dop = spdiags(gD, 0, n1, n1);
M.node_active = node_active;
M.S_pore = pore_storage(fcmplx, vcmplx, vol_cells, h_face);
end
