function R = solve_transport(W, opt)
% Steady flow scaled to the flow rate opt.Q_in, followed by transient
% advection-diffusion (backward Euler, one LU factorisation). The
% exponential and sigmoidal (erf) inlet histories, with the same 10-90%
% rise time opt.Tf_in, are solved as two right-hand sides.
%   W    fcmplx, vcmplx, fBO_1, fs1, fs2, k_H, D, vol_cells, h_face
%   opt  Q_in, Tf_in, dt, t_end, z_fraction, slab_halfwidth
%   R    threshold times [t10 t50 t90] at the active pore nodes (T_exp,
%        T_erf), median concentration in axial slabs (cmed), flow rate Q
%        and interstitial velocity v_p

fcmplx = W.fcmplx;
B   = W.fBO_1;
fs2 = W.fs2;
n0  = fcmplx(1).num(1).val;
n1  = fcmplx(2).num(2).val;

z = fcmplx(1).bc(:,3);
minz = min(z);
maxz = max(z);
L = maxz - minz;
tol = 1e-8;
bottom = z < (minz + tol);
top    = z > (maxz - tol);

% steady flow: unit inlet pressure, then scaled to Q_in
In_H = setdiff((1:n0)', find(top | bottom));
kd = full(spdiags(W.k_H, 0));
kd(kd == 0) = 1e-5;
A  = -B' * fs2 * spdiags(kd, 0, n1, n1) * B;
dA = decomposition(A(In_H,In_H), 'lu');

p = zeros(n0,1);
p(bottom) = 1;
p(In_H) = dA \ (-A(In_H,:) * p);
q = W.k_H * B * p;
e_top = abs(B) * double(top);
e_top(e_top > 1) = 0;
p_in = abs(opt.Q_in) / sum(e_top .* abs(fs2 * q));

p = zeros(n0,1);
p(bottom) = p_in;
p(In_H) = dA \ (-A(In_H,:) * p);
q = W.k_H * B * p;
Q = sum(e_top .* abs(fs2 * q));

% storage
S_pore = pore_storage(fcmplx, W.vcmplx, W.vol_cells, W.h_face);
pore = S_pore > 0;
S = full(diag(W.fs1));
S(pore) = S_pore(pore);
on_active = full(abs(B)' * double(full(spdiags(W.D, 0)) > 0)) > 0;

% transport operator: diffusion, advection (cup product), advective outflow
F = -(fs2 * q);
AC = B' * spdiags(F, 0, n1, n1) * (0.5 * abs(B));
bnet = B' * F;
q_out = zeros(n0,1);
q_out(top) = bnet(top);
Lop = -B' * fs2 * W.D * B + (AC - spdiags(q_out, 0, n0, n0));

Dn = find(bottom);
In = setdiff((1:n0)', Dn);
dt = opt.dt;
nt = round(opt.t_end / dt);
dM = decomposition(spdiags(S(In), 0, numel(In), numel(In)) - dt * Lop(In,In), 'lu');
b_in = dt * Lop(In,Dn) * ones(numel(Dn), 1);

% inlet histories
time = (1:nt)' * dt;
tau = opt.Tf_in / log(9);
sig = opt.Tf_in / (2 * sqrt(2) * erfinv(0.8));
t0  = 3 * sig;
g = @(t) 0.5 * erfc(-(t - t0) / (sqrt(2) * sig));
c_in = [1 - exp(-time / tau), (g(time) - g(0)) / (1 - g(0))];

% monitored nodes
act = find(pore & on_active & ~bottom & ~top);
pa = zeros(n0,1);
pa(In) = 1:numel(In);
pa = pa(act);
thr = [0.1 0.5 0.9];
T = nan(numel(act), 3, 2);

nZ = numel(opt.z_fraction);
slab = cell(nZ,1);
for j = 1:nZ
    zj = minz + opt.z_fraction(j) * L;
    slab{j} = find(pore & abs(z - zj) <= opt.slab_halfwidth * L);
end
cmed = nan(nt, nZ, 2, 'single');

c  = zeros(numel(In), 2);
cf = zeros(n0, 2);
c_prev = zeros(numel(act), 2);
for n = 1:nt
    c = dM \ (S(In) .* c + b_in * c_in(n,:));
    cf(In,:) = c;
    cf(Dn,:) = repmat(c_in(n,:), numel(Dn), 1);
    ca = c(pa,:);
    for s = 1:2
        for k = 1:3
            hit = isnan(T(:,k,s)) & ca(:,s) >= thr(k);
            if any(hit)
                c0 = c_prev(hit,s);
                c1 = ca(hit,s);
                T(hit,k,s) = (n - 1) * dt + (thr(k) - c0) ./ max(c1 - c0, eps) * dt;
            end
        end
        for j = 1:nZ
            cj = cf(slab{j},s);
            cj = cj(cj > 0);
            if ~isempty(cj), cmed(n,j,s) = median(cj); end
        end
    end
    c_prev = ca;
end

R.z     = z(act) - minz;
R.L     = L;
R.T_exp = T(:,:,1);
R.T_erf = T(:,:,2);
R.time  = time;
R.c_in  = c_in;
R.cmed  = cmed;
R.z_fraction = opt.z_fraction;
R.Q   = Q;
R.v_p = Q * L / sum(S_pore);
R.opt = opt;
end
