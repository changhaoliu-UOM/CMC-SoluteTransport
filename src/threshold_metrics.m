function Mt = threshold_metrics(R, which, mask)
% Arrival velocity, increase of the filling duration and median t50 of one run
% of solve_transport.m, over the monitored nodes with t10, t50 and t90 all
% reached and 0 < z < L (the definitions of cmc-py; Appendix A of the paper):
%   v_arr  1 / slope of the least-squares line t10 = a z + b
%   dTf    median of t90 - t10 in the last minus the first of 20 intervals of z/L
%   t50_median  median of t50 over the nodes that reached 0.5
% which = 'conv' (inlet history, default) or 'step'.
% mask (optional, logical over the monitored nodes): further restriction,
%   e.g. the nodes that reach 0.9 for two inlet histories.
if nargin < 2, which = 'conv'; end
if strcmp(which, 'step'), T = R.T_step; else, T = R.T_conv; end
zz = R.z_act - R.minz;  L = R.L;
ok = all(isfinite(T), 2) & zz > 0 & zz < L;
if nargin > 2 && ~isempty(mask), ok = ok & mask(:); end
c  = polyfit(zz(ok), T(ok,1), 1);
Mt.v_arr = 1 / c(1);
zl = zz / L;  tf = T(:,3) - T(:,1);
edges = linspace(0, 1, 21);  med = nan(20,1);
for b = 1:20
    if b < 20, inb = ok & zl >= edges(b) & zl < edges(b+1);
    else,      inb = ok & zl >= edges(b) & zl <= edges(b+1); end
    if any(inb), med(b) = median(tf(inb)); end
end
Mt.Tf_first = med(1);  Mt.Tf_last = med(20);  Mt.dTf = med(20) - med(1);
Mt.t50_median = median(T(isfinite(T(:,2)),2));
Mt.reached90 = mean(isfinite(T(:,3)));
Mt.n_monitored = numel(zz);
end
