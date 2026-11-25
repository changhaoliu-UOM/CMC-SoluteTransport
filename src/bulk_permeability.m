function k = bulk_permeability(M, mu)
% Bulk permeability (m^2) of the model M (build_model.m) from the steady flow
% with a unit pressure difference between the inlet (z = min) and the outlet
% (z = max) and impermeable lateral faces (Section 4.1 of the paper).
%   mu  dynamic viscosity (Pa s)
bc = M.fcmplx(1).bc;
[~, F, fi] = steady_flow(M.fcmplx, M.fBO_1, M.G);
zsf = abs(M.fBO_1) * double(fi.top);
zsf(zsf > 1) = 0;                                       % 1-cells with one end on the outlet
Q   = sum(zsf .* abs(F));                               % um^3/s per Pa
L   = max(bc(:,3)) - min(bc(:,3));
Axy = (max(bc(:,1)) - min(bc(:,1))) * (max(bc(:,2)) - min(bc(:,2)));
k   = (Q*1e-18) * mu * (L*1e-6) / (Axy*1e-12);          % m^2
end
