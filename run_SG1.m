% Transient solute transport in the reconstruction SG1 (sintered glass),
% computed with combinatorial mesh calculus for the exponential and
% sigmoidal (erf) inlet concentration histories at the two flow rates.
% Results: results/SG1_low.mat and results/SG1_high.mat (struct R).
clear; clc;

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));

tess_file = fullfile(root, 'data', 'SG1.tess');
V_vol   = readmatrix(fullfile(root, 'data', 'SG_voluminous_pores.csv'));
V_sheet = readmatrix(fullfile(root, 'data', 'SG_sheet_pores_cdf.csv'));
vf      = 0.278;                  % target porosity

Q_low = 1.42e8;                   % um^3/s (0.142 uL/s); high flow = 2 Q_low
flows = {'low', 'high'};
Tf_in = [182 92];                 % s, inlet 10-90% rise time [low high]

opt.dt    = 0.1;                  % s
opt.t_end = 800;                  % s
opt.z_fraction     = [0.1 0.3 0.5 0.7 0.9];
opt.slab_halfwidth = 0.02;        % fraction of the length

out_dir = fullfile(root, 'results');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

[vcmplx, fcmplx, fBO_1, fs1, fs2] = build_complexes(tess_file);
[k_H, D, vol_cells, h_face] = assign_pores(vcmplx, fcmplx, vf, V_vol, V_sheet);
W = struct('fcmplx', fcmplx, 'vcmplx', vcmplx, 'fBO_1', fBO_1, 'fs1', fs1, ...
           'fs2', fs2, 'k_H', k_H, 'D', D, 'vol_cells', vol_cells, 'h_face', h_face);

for iq = 1:numel(flows)
    opt.Q_in  = iq * Q_low;
    opt.Tf_in = Tf_in(iq);
    fprintf('SG1, %s flow ...\n', flows{iq});
    R = solve_transport(W, opt);
    save(fullfile(out_dir, sprintf('SG1_%s.mat', flows{iq})), 'R', '-v7.3');
end
