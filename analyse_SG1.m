% Arrival velocity and filling duration for the exponential and sigmoidal
% inlet histories, from the results of run_SG1.m. Nodes at which both
% histories reach 0.9 are used.
%   v_arrival  ordinary linear regression of t10 on z
%   Tf         median of t90 - t10 in 20 axial intervals; dTf = last - first
clear; clc;

root  = fileparts(mfilename('fullpath'));
flows = {'low', 'high'};
nBins = 20;
edges = linspace(0, 1, nBins + 1);
zmid  = (edges(1:end-1) + edges(2:end))' / 2;

summary = table();
profile = table(zmid, 'VariableNames', {'z_over_L'});
for iq = 1:numel(flows)
    S = load(fullfile(root, 'results', sprintf('SG1_%s.mat', flows{iq})), 'R');
    R = S.R;
    zL = R.z / R.L;
    ok = all(isfinite(R.T_exp), 2) & all(isfinite(R.T_erf), 2) & zL > 0 & zL < 1;

    T  = {R.T_exp, R.T_erf};
    v  = zeros(1,2);
    Tf = nan(nBins, 2);
    for s = 1:2
        p = polyfit(R.z(ok), T{s}(ok,1), 1);
        v(s) = 1 / p(1);
        tf = T{s}(:,3) - T{s}(:,1);
        for b = 1:nBins
            if b < nBins
                in = ok & zL >= edges(b) & zL < edges(b+1);
            else
                in = ok & zL >= edges(b) & zL <= edges(b+1);
            end
            Tf(b,s) = median(tf(in));
        end
    end

    summary = [summary; table(flows(iq), R.Q, R.v_p, v(1), v(2), ...
        Tf(end,1) - Tf(1,1), Tf(end,2) - Tf(1,2), ...
        'VariableNames', {'flow', 'Q_um3_s', 'v_p', 'v_arrival_exp', 'v_arrival_erf', ...
                          'dTf_exp', 'dTf_erf'})]; %#ok<AGROW>
    profile.(['Tf_exp_' flows{iq}]) = Tf(:,1);
    profile.(['Tf_erf_' flows{iq}]) = Tf(:,2);
end

disp(summary);
writetable(summary, fullfile(root, 'results', 'SG1_summary.csv'));
writetable(profile, fullfile(root, 'results', 'SG1_filling_duration.csv'));
