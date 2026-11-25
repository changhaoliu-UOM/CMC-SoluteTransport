%% run_transport.m
% Transport calculations of the paper for every realisation in data/
% (data/SG<id>.tess, data/BS<id>.tess and the extended model data/SGL1.tess):
% pore assignment, bulk permeability, and steady flow with transient
% advection-diffusion for the step inlet and the inlet concentration history
% at the two flow rates (composition model, diagonal Hodge star, exponential
% fitting, dt = 0.1 s).
%
% Output  results/runs/<MAT><id>_<flow>.mat   struct R (solve_transport.m),
%                                             R.summary (analyse_transport.m)
%         results/summary.csv                 one row per run
% Saved runs are read, not recomputed, so the script can be stopped and
% started again.
clear; clc;
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));

materials = {'SG', 'BS', 'SGL'};
out_dir = fullfile(root, 'results', 'runs');
if ~isfolder(out_dir), mkdir(out_dir); end

rows = {};
for im = 1:numel(materials)
    cfg = config_material(materials{im});
    ids = find_realisations(cfg);
    if isempty(ids)
        fprintf('%s: no tessellation in data/, skipped\n', cfg.material);
        continue
    end
    for rid = ids
        M = [];
        k = NaN;
        for iq = 1:numel(cfg.flows)
            fl = cfg.flows{iq};
            f = fullfile(out_dir, sprintf('%s%d_%s.mat', cfg.material, rid, fl));
            if isfile(f)
                S = load(f, 'R');
                R = S.R;
                fprintf('>> %s read\n', f);
            else
                if isempty(M)
                    fprintf('\n==================== %s%d ====================\n', cfg.material, rid);
                    M = build_model(cfg, rid);
                    k = bulk_permeability(M, cfg.mu);
                    fprintf('bulk permeability %.4e m^2\n', k);
                end
                opt = struct('Q_in', iq * cfg.Q_low, 'Tf_in', cfg.Tf_in(iq), 'dt', cfg.dt, ...
                             't_end', cfg.t_end, 'z_fraction', cfg.z_fraction, ...
                             'slab_halfwidth', cfg.slab_halfwidth);
                fprintf('\n---- %s%d, %s flow (Q = %.3e um^3/s, Tf_in = %g s) ----\n', ...
                    cfg.material, rid, fl, opt.Q_in, opt.Tf_in);
                R = solve_transport(M, opt);
                R.summary  = analyse_transport(R, cfg.Dm);
                R.k_bulk   = k;
                R.pores    = M.pores;
                R.material = cfg.material;
                R.id       = rid;
                R.flow     = fl;
                save(f, 'R', '-v7.3');
            end
            A = R.summary;
            rows(end+1,:) = {cfg.material, rid, fl, R.Q_final, R.k_bulk, ...
                R.pores.vf_vol + R.pores.vf_sheet, R.v_storage, A.v_arr_step, A.v_arr_conv, ...
                cfg.v_meas(iq), 100*(A.v_arr_conv/cfg.v_meas(iq) - 1), ...
                A.Tf_conv_in, A.Tf_conv_out, A.dTf_conv, cfg.dTf_meas(iq), ...
                A.Tf_step_out, A.DL_over_Dm, min(R.cmin), max(R.cmax), max(R.mass_res), ...
                R.frac90(1), R.frac90(2)}; %#ok<SAGROW>
            fprintf('%s%d %s: v_arrival %.2f um/s (step %.2f, v_p %.2f, measured %.1f); Tf %.1f -> %.1f s; D_L,eq/D_m %.2f\n', ...
                cfg.material, rid, fl, A.v_arr_conv, A.v_arr_step, R.v_storage, cfg.v_meas(iq), ...
                A.Tf_conv_in, A.Tf_conv_out, A.DL_over_Dm);
        end
    end
end
if isempty(rows), return; end

T = cell2table(rows, 'VariableNames', {'material','id','flow','Q_um3_s','k_bulk_m2', ...
    'porosity','v_p','v_arr_step','v_arr_history','v_meas','v_diff_pct', ...
    'Tf_history_in','Tf_history_out','dTf_history','dTf_meas','Tf_step_out','DL_over_Dm', ...
    'c_min','c_max','mass_residual','reached90_step','reached90_history'});
writetable(T, fullfile(root, 'results', 'summary.csv'));

% ensemble mean and standard deviation (Table 6 and Sections 4.1-4.5)
fprintf('\n%-4s %-5s %3s %11s %16s %16s %16s %16s %14s\n', 'mat', 'flow', 'n', 'k (m^2)', ...
    'v_step (um/s)', 'v_hist (um/s)', 'Tf_out,hist (s)', 'Tf_out,step (s)', 'D_L,eq/D_m');
for im = 1:numel(materials)
    for iq = 1:2
        m = strcmp(T.material, materials{im}) & strcmp(T.flow, cfg.flows{iq});
        if ~any(m), continue; end
        fprintf('%-4s %-5s %3d %11.3e %8.2f +- %5.2f %8.2f +- %5.2f %8.1f +- %5.1f %8.1f +- %5.1f %7.2f +- %4.2f\n', ...
            materials{im}, cfg.flows{iq}, nnz(m), mean(T.k_bulk_m2(m)), ...
            mean(T.v_arr_step(m)), std(T.v_arr_step(m)), mean(T.v_arr_history(m)), std(T.v_arr_history(m)), ...
            mean(T.Tf_history_out(m)), std(T.Tf_history_out(m)), mean(T.Tf_step_out(m)), std(T.Tf_step_out(m)), ...
            mean(T.DL_over_Dm(m)), std(T.DL_over_Dm(m)));
    end
end
fprintf('\n%s\n', fullfile(root, 'results', 'summary.csv'));
