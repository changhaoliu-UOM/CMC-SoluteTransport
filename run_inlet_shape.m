%% run_inlet_shape.m
% Sensitivity to the shape of the inlet concentration history (Appendix A of
% the paper): exponential and sigmoidal (erf) histories with the same 10-90%
% rise time, SG1 and BS1, both flow rates, dt = 0.1 s, t = 800 s (to allow for
% the delay t0 of the sigmoidal history). Arrival velocity: regression of t10
% on z with an intercept; dTf: median t90 - t10 in the last minus the first of
% 20 axial intervals; both over the nodes that reach 0.9 for both histories.
%
% Output  results/inlet_shape/<MAT>1_<flow>_<history>.mat
%         results/inlet_shape/inlet_shape.csv
clear; clc;
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));

materials = {'SG', 'BS'};
rid   = 1;
hists = {'exp', 'erf'};
out_dir = fullfile(root, 'results', 'inlet_shape');
if ~isfolder(out_dir), mkdir(out_dir); end

rows = {};
for im = 1:numel(materials)
    cfg = config_material(materials{im});
    M = [];
    for iq = 1:numel(cfg.flows)
        fl = cfg.flows{iq};
        Rh = cell(1,2);
        for ih = 1:2
            f = fullfile(out_dir, sprintf('%s%d_%s_%s.mat', cfg.material, rid, fl, hists{ih}));
            if isfile(f)
                S = load(f, 'R');  Rh{ih} = S.R;  continue
            end
            if isempty(M), M = build_model(cfg, rid); end
            opt = struct('Q_in', iq * cfg.Q_low, 'Tf_in', cfg.Tf_in(iq), 'dt', cfg.dt, 't_end', 800, ...
                         'z_fraction', cfg.z_fraction, 'slab_halfwidth', cfg.slab_halfwidth, ...
                         'history', hists{ih});
            fprintf('\n---- %s%d, %s flow, %s history ----\n', cfg.material, rid, fl, hists{ih});
            R = solve_transport(M, opt);
            save(f, 'R', '-v7.3');
            Rh{ih} = R;
        end
        both = isfinite(Rh{1}.T_conv(:,3)) & isfinite(Rh{2}.T_conv(:,3));   % 0.9 reached for both
        Me = threshold_metrics(Rh{1}, 'conv', both);
        Ms = threshold_metrics(Rh{2}, 'conv', both);
        rows(end+1,:) = {cfg.material, rid, fl, Rh{1}.opt.Q_in, Me.v_arr, Ms.v_arr, ...
            100*(Ms.v_arr/Me.v_arr - 1), Me.dTf, Ms.dTf, 100*(Ms.dTf/Me.dTf - 1), ...
            nnz(both), numel(both), Rh{2}.mass_res(2), Rh{2}.cmin(2), Rh{2}.cmax(2)}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'material','id','flow','Q_um3_s','v_arr_exp','v_arr_erf', ...
    'v_change_pct','dTf_exp','dTf_erf','dTf_change_pct','n_both_reach90','n_monitored', ...
    'mass_res_erf','c_min_erf','c_max_erf'});
writetable(T, fullfile(out_dir, 'inlet_shape.csv'));
disp(T);
fprintf('\n%s\n', fullfile(out_dir, 'inlet_shape.csv'));
