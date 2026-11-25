function A = analyse_transport(R, Dm)
% Arrival velocity, filling duration and equivalent dispersion coefficient of
% one run of solve_transport.m (Sections 4.3-4.5 of the paper).
%   Arrival velocity: regression of t10 on z, through the origin for the step
%     inlet (A.v_arr_step) and with an intercept for the inlet history (A.v_arr_conv).
%   Filling duration Tf = t90 - t10: median and 10th/90th percentiles in 20
%     axial intervals (A.bins_step, A.bins_conv); A.dTf_conv = increase from the
%     first to the last interval for the inlet history.
%   D_L,eq: constant dispersion coefficient for which the Ogata-Banks solution
%     with the interstitial velocity gives the median step-inlet Tf of the nodes
%     with z/L > 0.95 (A.DL, A.DL_over_Dm = A.DL / Dm); A.OB: Tf(z) of the full
%     Ogata-Banks solution and of the asymptotic Fickian relation
%     3.6248 sqrt(D_L z / v^3) with that D_L,eq.
[A.v_arr_step, ~] = arrival_fit(R, R.T_step(:,1));
[~, A.v_arr_conv] = arrival_fit(R, R.T_conv(:,1));
edges = linspace(0, 1, 21);
A.bins_step = bin_filling(R, R.T_step, edges, [10 90]);
A.bins_conv = bin_filling(R, R.T_conv, edges, [10 90]);
A.Tf_step_in  = A.bins_step.median(1);   A.Tf_step_out = A.bins_step.median(end);
A.Tf_conv_in  = A.bins_conv.median(1);   A.Tf_conv_out = A.bins_conv.median(end);
A.dTf_conv    = A.Tf_conv_out - A.Tf_conv_in;
A.band_step_out = A.bins_step.p_hi(end) - A.bins_step.p_lo(end);
A.band_conv_out = A.bins_conv.p_hi(end) - A.bins_conv.p_lo(end);
Tf = R.T_step(:,3) - R.T_step(:,1);
zL = (R.z_act - R.minz) / R.L;
band = isfinite(Tf) & Tf >= 0 & zL > 0.95 & zL < 1;
Tf_out = median(Tf(band));  z_star = median(R.z_act(band) - R.minz);
v = R.v_storage;
DL0 = v^3 * Tf_out^2 / (3.6248^2 * z_star);
A.DL = fit_ogata_banks_DL(z_star, v, Tf_out, DL0, R.L/v);
A.DL_over_Dm = A.DL / Dm;
A.v_p = v;
A.z_star = z_star;  A.Tf_step_target = Tf_out;
% full Ogata-Banks solution and asymptotic Fickian relation with v_p and D_L,eq
zg = linspace(0.01, 1, 80)' * R.L;
Tfull = arrayfun(@(zz) ogata_banks_filling_duration(zz, v, A.DL, R.L/v), zg);
A.OB = table(zg / R.L, Tfull, 3.6248 * sqrt(A.DL * zg / v^3), ...
             'VariableNames', {'z_over_L','Tf_full','Tf_asymptotic'});
end

function [v_zero, v_free] = arrival_fit(R, t10)
ok = isfinite(t10);
x = R.z_act(ok) - R.minz;  y = t10(ok);
v_zero = 1 / (sum(x.*y) / sum(x.^2));
p = polyfit(x, y, 1);  v_free = 1 / p(1);
end

function B = bin_filling(R, T, edges, pct)
nB = numel(edges) - 1;
Tf = T(:,3) - T(:,1);
zL = (R.z_act - R.minz) / R.L;
ok = isfinite(Tf) & Tf >= 0 & zL > 0 & zL < 1;
n = zeros(nB,1);  med = nan(nB,1);  lo = med;  hi = med;  p25 = med;  p75 = med;
for j = 1:nB
    if j < nB, in = ok & zL >= edges(j) & zL < edges(j+1);
    else,      in = ok & zL >= edges(j) & zL <= edges(j+1); end
    n(j) = nnz(in);
    if n(j) == 0, continue; end
    med(j) = median(Tf(in));
    lo(j)  = prctile(Tf(in), pct(1));  hi(j)  = prctile(Tf(in), pct(2));
    p25(j) = prctile(Tf(in), 25);      p75(j) = prctile(Tf(in), 75);
end
zc = 0.5*(edges(1:end-1) + edges(2:end))';
B = table(zc, n, med, lo, hi, p25, p75, 'VariableNames', ...
    {'z_over_L','n_nodes','median','p_lo','p_hi','p25','p75'});
end

function C = ogata_banks_concentration_fit(z,t,v,DL)
a = (z-v.*t) ./ (2.*sqrt(DL.*t));
b = (z+v.*t) ./ (2.*sqrt(DL.*t));
C = 0.5 .* (erfc(a) + exp(-a.^2) .* erfcx(b));
end

function tp = ogata_banks_threshold_time_fit(z,v,DL,threshold,time_scale)
t_lower = 1e-10;
f = @(tt) ogata_banks_concentration_fit(z,tt,v,DL) - threshold;
if f(t_lower) >= 0, tp = t_lower; return; end
t_upper = max(10.*time_scale, 2.*z./v);
f_upper = f(t_upper);  n_expand = 0;
while (~isfinite(f_upper) || f_upper < 0) && n_expand < 50
    t_upper = 2.*t_upper;  f_upper = f(t_upper);  n_expand = n_expand+1;
end
tp = fzero(f, [t_lower,t_upper]);
end

function Tf = ogata_banks_filling_duration(z,v,DL,time_scale)
Tf = ogata_banks_threshold_time_fit(z,v,DL,0.90,time_scale) - ...
     ogata_banks_threshold_time_fit(z,v,DL,0.10,time_scale);
end

function DL_fit = fit_ogata_banks_DL(z,v,Tf_target,DL_initial,time_scale)
objective = @(lD) ogata_banks_filling_duration(z,v,10.^lD,time_scale) - Tf_target;
logDL0 = log10(DL_initial);  span = 1;  found = false;
for k = 1:12
    lo = logDL0-span;  hi = logDL0+span;
    flo = objective(lo);  fhi = objective(hi);
    if isfinite(flo) && isfinite(fhi) && flo*fhi <= 0, found = true; break; end
    span = 1.5*span;
end
if ~found, error('Unable to bracket the Ogata-Banks DL solution.'); end
DL_fit = 10.^fzero(objective, [lo,hi]);
end
