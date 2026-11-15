function [k_H, D, vol_cells, h_face] = assign_pores(vcmplx, fcmplx, vf, V_vol, V_sheet)
% Assigns voluminous pores to 3-cells and sheet-like pores to 2-cells of M
% and builds the hydraulic (k_H) and diffusion (D) coefficients on the
% 1-cells of K.
%   vf       target porosity
%   V_vol    volumes of the voluminous pores (um^3)
%   V_sheet  [volume (um^3), cumulative probability] of the sheet-like pores
%   vol_cells  3-cells of M assigned as voluminous pores
%   h_face     aperture of each 2-cell of M (zero if not a sheet-like pore)

mu     = 1e-3;                    % Pa s
D_pore = 1770;                    % um^2/s

bc = fcmplx(1).bc;
V_box = max(bc(:,1)) * max(bc(:,2)) * max(bc(:,3));
nC = vcmplx(4).num(4).val;
nF = vcmplx(3).num(3).val;
nE = vcmplx(2).num(2).val;
nmax = max(arrayfun(@(i) nnz(vcmplx(4).bndop(3).indx(i,:)), 1:nC));

grain = vcmplx(4).vvol;
area  = vcmplx(3).vvol;
thick  = zeros(nF,1);
h_face = zeros(nF,1);
h_edge = zeros(nE,1);
k_cell = zeros(nC, nmax);

% voluminous pores: 3-cells with the closest volumes
vol_cells = zeros(size(V_vol));
used = false(size(grain));
for i = 1:numel(V_vol)
    d = abs(grain - V_vol(i));
    d(used) = Inf;
    [~, j] = min(d);
    vol_cells(i) = j;
    used(j) = true;
end
vf_vol = sum(grain(vol_cells)) / V_box;
vol_cells = sort(vol_cells);

for cell_id = vol_cells(:)'
    x3 = vcmplx(4).bc(cell_id,1:3);
    f = nonzeros(vcmplx(4).bndop(3).indx(cell_id,:));
    for j = 1:numel(f)
        k_cell(cell_id,j) = (area(f(j)) / pi) / (8 * mu);
        thick(f(j)) = norm(x3 - vcmplx(3).bc(f(j),1:3));
    end
    e = nonzeros(vcmplx(4).bndop(2).indx(cell_id,:));
    for j = 1:numel(e)
        h_edge(e(j)) = norm(x3 - vcmplx(2).bc(e(j),1:3));
    end
end

% sheet-like pores: volumes sampled from the cumulative distribution
b = area;
for cell_id = vol_cells(:)'
    b(nonzeros(vcmplx(4).bndop(3).indx(cell_id,:))) = 0;
end
V_s = [];
V_acc = 0;
vf_s = 0;
rng(1);
while vf_s < vf - vf_vol
    v = interp1(V_sheet(:,2), V_sheet(:,1), rand, 'linear', 'extrap');
    V_acc = V_acc + v;
    V_s = [V_s; v]; %#ok<AGROW>
    vf_s = V_acc / V_box;
end

% 2-cells for the sheet-like pores
nS = numel(V_s);
sel_area = nan(nS,1);
sel_face = nan(nS,1);
for i = 1:nS
    lo = (pi/4) * (20 * V_s(i) / pi)^(2/3);
    hi = (pi/4) * (50 * V_s(i) / pi)^(2/3);
    idx = find(isfinite(b) & b > lo & b < hi);
    if ~isempty(idx)
        [a, k] = min(b(idx));
    else
        idx = find(isfinite(b) & b > 0);
        if isempty(idx), continue; end
        [a, k] = max(b(idx));
    end
    sel_area(i) = a;
    sel_face(i) = idx(k);
    b(idx(k)) = -Inf;
end
for i = 1:nS
    if isnan(sel_face(i)) || isnan(sel_area(i)), continue; end
    f = sel_face(i);
    thick(f)  = V_s(i) / sel_area(i);
    h_face(f) = V_s(i) / sel_area(i);
end

% apertures of the 1-cells of M bounding pores
for i = 1:nF
    e = nonzeros(vcmplx(3).bndop(2).indx(i,:));
    if thick(i) ~= 0
        h_edge(e) = max(h_edge(e), h_face(i));
    end
end
h_sum = zeros(nE,1);
h_cnt = zeros(nE,1);
m = h_edge > 0;
h_sum(m) = h_edge(m);
h_cnt(m) = 1;
for i = 1:nF
    if h_face(i) > 0
        e = nonzeros(vcmplx(3).bndop(2).indx(i,:));
        h_sum(e) = h_sum(e) + h_face(i);
        h_cnt(e) = h_cnt(e) + 1;
    end
end
m = h_cnt > 0;
h_edge(m) = h_sum(m) ./ h_cnt(m);

% hydraulic coefficients: (a0<a1), (a1<a2) or (a1<b2), (a2<b3)
k1 = ((h_edge / 2).^2) / (8 * mu);
k1 = reshape([k1 k1]', [], 1);

c2 = thick.^2 / (12 * mu);
k2 = ones(sum(vcmplx(3).num(2).val), 1);
ct = 0;
for i = 1:nF
    for j = 1:size(vcmplx(3).bndop(2).indx, 2)
        if vcmplx(3).bndop(2).indx(i,j) ~= 0
            ct = ct + 1;
            k2(ct) = c2(i);
        end
    end
end

k3 = ones(sum(vcmplx(4).num(3).val), 1);
ct = 0;
for i = 1:nC
    for j = 1:nnz(vcmplx(4).bndop(3).indx(i,:))
        if vcmplx(4).bndop(3).indx(i,j) ~= 0
            ct = ct + 1;
            k3(ct) = k_cell(i,j);
        end
    end
end

k = [k1; k2; k3];
n1 = fcmplx(2).num(2).val;
k_H = spdiags(k, 0, n1, n1);
D   = spdiags(D_pore * double(k > 0), 0, n1, n1);
end
