function [p, F, info] = steady_flow(fcmplx, B, G)
% Steady flow on the vertices of K with p = 1 at z = min(z), p = 0 at
% z = max(z) and impermeable lateral faces. The flow rates on the 1-cells are
% F = -G B p, with the same conductances G as the pressure equation, so the
% flow balances at every vertex to round-off. Vertices not connected to the
% inlet or outlet through 1-cells of nonzero conductance carry no flow
% (p = 0 there, info.active = false).
%   B  boundary operator of the 1-cells of K (n1 x n0)
%   G  diagonal hydraulic conductances (n1 x n1)
small = 1e-8;
z  = fcmplx(1).bc(:,3);
n0 = numel(z);
bottom = z < min(z) + small;
top    = z > max(z) - small;
fixed  = bottom | top;
p = zeros(n0,1);
p(bottom) = 1;

A   = -B' * G * B;
Adj = A - spdiags(diag(A), 0, n0, n0);
[i, j] = find(Adj);
up  = i < j;
comp = conncomp(graph(i(up), j(up), ones(nnz(up),1), n0));
comp = comp(:);
anchored = false(max(comp),1);
anchored(comp(fixed)) = true;
active = anchored(comp) & full(sum(abs(A), 2)) > 0;
active(fixed) = true;

free = find(active & ~fixed);
p(free) = decomposition(A(free,free), 'lu') \ (-A(free,fixed) * p(fixed));
F = -(G * (B * p));
info = struct('active', active, 'n_inactive', nnz(~active), 'n_unknown', numel(free), ...
              'bottom', bottom, 'top', top);
end
