%% run_inlet_shape.m
% Sensitivity to the shape of the inlet concentration history (Appendix A of
% the paper): exponential and sigmoidal (erf) histories with the same 10-90%
% rise time, SG1 and BS1, both flow rates, dt = 0.1 s, t = 800 s (to allow for
% the delay t0 of the sigmoidal history). Arrival velocity and increase of the
% filling duration as in Sections 4.3 and 4.4 (analyse_transport.m), each
% history over its own nodes.
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
        Ae = analyse_transport(Rh{1}, cfg.Dm);
        As = analyse_transport(Rh{2}, cfg.Dm);
        rows(end+1,:) = {cfg.material, rid, fl, Rh{1}.opt.Q_in, Ae.v_arr_conv, As.v_arr_conv, ...
            100*(As.v_arr_conv/Ae.v_arr_conv - 1), Ae.dTf_conv, As.dTf_conv, ...
            100*(As.dTf_conv/Ae.dTf_conv - 1), Rh{2}.mass_res(2), Rh{2}.cmin(2), Rh{2}.cmax(2)}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'material','id','flow','Q_um3_s','v_arr_exp','v_arr_erf', ...
    'v_change_pct','dTf_exp','dTf_erf','dTf_change_pct','mass_res_erf','c_min_erf','c_max_erf'});
writetable(T, fullfile(out_dir, 'inlet_shape.csv'));
disp(T);
fprintf('\n%s\n', fullfile(out_dir, 'inlet_shape.csv'));
