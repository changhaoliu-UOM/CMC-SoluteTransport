function R = solve_transport(M, opt)
% Steady flow and transient advection-diffusion on K (Section 3.3 of the paper).
%   M    model of build_model.m
%   opt  Q_in            imposed flow rate (um^3/s)
%        Tf_in           10-90% rise time of the inlet history (s)
%        dt, t_end       time step and end time (s)
%        z_fraction      monitoring planes (z/L) for the slab medians
%        slab_halfwidth  half width of the monitoring slabs (fraction of L)
%        art_scale       (optional, default 1) factor on the diffusion added by
%                        the exponential fitting (1 = exponential fitting;
%                        other values only to estimate its numerical diffusion)
%        history         (optional) 'exp' (default) or 'erf'
%
% Flow: steady_flow.m with a unit pressure difference, scaled to Q_in (the
% problem is linear).
% Transport: backward Euler; advection by exponential fitting
% (Scharfetter-Gummel): on a 1-cell e with flow rate Q_e from node i to node j
% and diffusive conductance G_e, the transport rate is
%     J_e = Q_e (c_i + c_j)/2 - G_e alpha(P_e) (c_j - c_i),
%     alpha(P) = (P/2) coth(P/2),  P_e = |Q_e| / G_e,
% i.e. equal (cup-product) weights plus the added diffusion G_e (alpha - 1).
% The outlet carries the advective outflow; the matrix is factorised once.
%
% Two inlet conditions are solved together as two right-hand sides:
%   1 (step):    c_in = 1
%   2 (history): c_in = 1 - exp(-t/tau), tau = Tf_in/ln 9         ('exp')
%                or the sigmoidal history of the same 10-90% rise time
%                f = (g(t) - g(0))/(1 - g(0)), g = erfc(-(t - t0)/(sqrt(2) sigma))/2,
%                sigma = Tf_in/(2 sqrt(2) erfinv(0.8)), t0 = 3 sigma   ('erf', Appendix A)
% In R, fields ending in _step refer to inlet 1 and fields ending in _conv to
% inlet 2 (the convolution of the step response with the inlet history).
%
% R.T_step, R.T_conv  first times at which c reaches 0.1, 0.5, 0.9 at the
%                     active interior pore nodes (rows), at heights R.z_act
% R.cmed_step, R.cmed_conv  median c in the slabs at R.z_fraction, every step
% R.v_storage         interstitial velocity Q L / pore volume (um/s)
t_all = tic;
if ~isfield(opt, 'art_scale'), opt.art_scale = 1; end
if ~isfield(opt, 'history'),   opt.history = 'exp'; end
small = 1e-8;
fcmplx = M.fcmplx;
B  = M.fBO_1;
z  = fcmplx(1).bc(:,3);
minz = min(z);  maxz = max(z);  L = maxz - minz;
n0 = fcmplx(1).num(1).val;
n1 = fcmplx(2).num(2).val;
bottom = z < (minz + small);
top    = z > (maxz - small);

%% ---- steady flow scaled to Q_in -----------------------------------------
[~, F, finfo] = steady_flow(fcmplx, B, M.G);
zsf_e = abs(B) * double(top);
zsf_e(zsf_e > 1) = 0;                              % 1-cells with one end on the outlet
Q_trial = sum(zsf_e .* abs(F));
if Q_trial < small, error('solve_transport: no flow through the outlet'); end
p_in = abs(opt.Q_in) / Q_trial;
F = F * p_in;
Q_final = sum(zsf_e .* abs(F));

%% ---- storage ---------------------------------------------------------------
S_pore = M.S_pore;
pore_node = S_pore > 0;
S = M.fvol0;                                       % solid nodes are isolated (no conductance)
S(pore_node) = S_pore(pore_node);                  % and stay at c = 0
v_storage = Q_final * L / sum(S_pore);

%% ---- transport operator -------------------------------------------------
AD_T = -B' * M.Dop * B;                            % diffusion
dD   = full(diag(M.Dop));                          % diffusive conductance of every 1-cell
Pe_e = abs(F) ./ max(dD, realmin);                 % Peclet number of every 1-cell
pore_e = dD > 0;
Pe_stats = [median(Pe_e(pore_e)), prctile(Pe_e(pore_e), 95), mean(Pe_e(pore_e) > 2)];

art = zeros(n1,1);                                 % diffusion added by the exponential fitting
art(pore_e) = opt.art_scale * dD(pore_e) .* (bernoulli_alpha(Pe_e(pore_e)) - 1);
AD_T = AD_T - B' * spdiags(art, 0, n1, n1) * B;

bnet  = B' * F;
q_out = zeros(n0,1);
q_out(top) = bnet(top);                            % advective outflow at the outlet
[ie, jn] = find(B);
Wc = sparse(ie, jn, 0.5 * ones(numel(ie),1), n1, n0);   % equal weights of the two ends
AC = B' * spdiags(F, 0, n1, n1) * Wc;
AA_T = -(AC - spdiags(q_out, 0, n0, n0));
art_stats = [sum(art(pore_e)) / sum(dD(pore_e)), mean(art(pore_e) > 0.01 * dD(pore_e))];
Lop = AD_T - AA_T;

Dn = find(bottom);                                 % inlet nodes (Dirichlet)
In = setdiff((1:n0)', Dn);
dt = opt.dt;
nt = round(opt.t_end / dt);
A  = spdiags(S(In), 0, numel(In), numel(In)) - dt * Lop(In,In);
dA = decomposition(A, 'lu');
B_in = dt * Lop(In,Dn) * ones(numel(Dn),1);        % inlet coupling for c_in = 1

time = (1:nt)' * dt;
if strcmp(opt.history, 'erf')
    sig = opt.Tf_in / (2*sqrt(2)*erfinv(0.8));
    t0  = 3*sig;
    g   = @(t) 0.5 * erfc(-(t - t0) / (sqrt(2)*sig));
    c_in = [ones(nt,1), (g(time) - g(0)) / (1 - g(0))];
else
    tau_in = opt.Tf_in / log(9);
    c_in = [ones(nt,1), 1 - exp(-time/tau_in)];
end

%% ---- monitored quantities ------------------------------------------------
act = find(pore_node & M.node_active(:) & ~bottom & ~top);   % active interior pore nodes
nA  = numel(act);
thr = [0.1 0.5 0.9];
T_step = nan(nA,3);
T_conv = nan(nA,3);

nP = numel(opt.z_fraction);
slab = cell(nP,1);
for j = 1:nP
    zp = minz + opt.z_fraction(j)*L;
    slab{j} = find(pore_node & abs(z - zp) <= opt.slab_halfwidth*L);
end
cmed = nan(nt, nP, 2, 'single');

cI   = zeros(numel(In), 2);                        % interior concentrations
iIn  = zeros(n0,1);
iIn(In) = 1:numel(In);
pa   = iIn(act);
c_prev = zeros(nA,2);
cmin = [0 0];  cmax = [0 0];  cmin_all = [0 0];  cmax_all = [0 0];
cumIn = [0 0];  cumR = [0 0];  res_max = [0 0];
mask_I = pore_node(In);
LD = Lop(Dn,:);
colsum = full(sum(Lop,1))';
cfull = zeros(n0,2);

%% ---- time loop -----------------------------------------------------------
for n = 1:nt
    rhs = S(In).*cI + B_in * c_in(n,:);
    cI  = dA \ rhs;
    cfull(In,:) = cI;
    cfull(Dn,:) = repmat(c_in(n,:), numel(Dn), 1);

    % threshold times, by linear interpolation between steps
    ca = cI(pa,:);
    for k = 1:3
        for s = 1:2
            if s == 1, Tk = T_step(:,k); else, Tk = T_conv(:,k); end
            hit = isnan(Tk) & ca(:,s) >= thr(k);
            if any(hit)
                c0 = c_prev(hit,s);  c1 = ca(hit,s);
                Tk(hit) = (n-1)*dt + (thr(k)-c0)./max(c1-c0,eps)*dt;
            end
            if s == 1, T_step(:,k) = Tk; else, T_conv(:,k) = Tk; end
        end
    end
    c_prev = ca;

    % slab medians (pore nodes with c > 0)
    for s = 1:2
        for j = 1:nP
            cj = cfull(slab{j},s);  cj = cj(cj > 0);
            if ~isempty(cj), cmed(n,j,s) = median(cj); end
        end
    end

    % concentration range and solute mass balance
    cmin = min(cmin, min(cI(mask_I,:),[],1));
    cmax = max(cmax, max(cI(mask_I,:),[],1));
    cmin_all = min(cmin_all, min(cI,[],1));
    cmax_all = max(cmax_all, max(cI,[],1));
    q_in  = -sum(LD * cfull, 1);
    q_ou  =  sum(q_out .* cfull, 1);
    gen   =  colsum' * cfull;                      % zero for a conservative operator
    cumIn = cumIn + dt*q_in;
    cumR  = cumR  + dt*(gen + q_ou);
    res_max = max(res_max, abs(cumR) ./ max(cumIn, eps));

    if mod(n, round(50/dt)) == 0
        fprintf('   t = %6.1f / %g s\n', n*dt, opt.t_end);
    end
end

%% ---- output ----------------------------------------------------------------
R.z_act      = z(act);
R.minz       = minz;
R.L          = L;
R.T_step     = T_step;          % [t10 t50 t90], step inlet
R.T_conv     = T_conv;          % [t10 t50 t90], inlet history
R.time       = time;
R.cmed_step  = cmed(:,:,1);
R.cmed_conv  = cmed(:,:,2);
R.z_fraction = opt.z_fraction;
R.Q_final    = Q_final;         % um^3/s
R.p_in       = p_in;            % Pa
R.v_storage  = v_storage;       % um/s
R.pore_volume = sum(S_pore);    % um^3
R.n_inactive_flow = finfo.n_inactive;
R.cmin       = cmin;            % pore nodes [step history]
R.cmax       = cmax;
R.cmin_all   = cmin_all;        % all interior nodes
R.cmax_all   = cmax_all;
R.Pe_stats   = Pe_stats;        % pore 1-cells: [median, 95th percentile, fraction > 2]
R.art_stats  = art_stats;       % added / physical diffusion, fraction of 1-cells with > 1 %
R.mass_res   = res_max;         % largest relative mass-balance residual [step history]
R.frac90     = [mean(isfinite(T_step(:,3))) mean(isfinite(T_conv(:,3)))];
R.opt        = opt;

fprintf('   Q = %.4e um^3/s, v_p = %.3f um/s, %d monitored nodes\n', Q_final, v_storage, nA);
fprintf('   Peclet number of the pore 1-cells: median %.2f, 95th percentile %.2f\n', Pe_stats(1:2));
fprintf('   c range: step [%.2e, %.5f], history [%.2e, %.5f]; mass-balance residual %.1e, %.1e\n', ...
    cmin(1), cmax(1), cmin(2), cmax(2), res_max(1), res_max(2));
fprintf('   nodes reaching 0.9: step %.2f %%, history %.2f %% (%.0f s)\n', 100*R.frac90, toc(t_all));
end

% -------------------------------------------------------------------------
function a = bernoulli_alpha(Pe)
% alpha = (Pe/2) coth(Pe/2), even in Pe; 1 + Pe^2/12 for small Pe
y = 0.5 * abs(Pe);
a = ones(size(y));
big = y > 1e-4;
a(big)  = y(big) ./ tanh(y(big));
a(~big) = 1 + y(~big).^2 / 3;
end
