function [vol_cells, h_face, info] = assign_pores(vcmplx, fcmplx, vf, V_vol, V_sheet)
% Assigns the voluminous pores to 3-cells and the sheet-like pores to 2-cells
% of the physical complex M (Section 3.2 of the paper).
%   vf       target porosity
%   V_vol    volumes of the voluminous pores (um^3)
%   V_sheet  [volume (um^3), cumulative probability] of the sheet-like pores
%
%   vol_cells  3-cells of M assigned as voluminous pores (sorted)
%   h_face     aperture of every 2-cell of M (zero if not a sheet-like pore)
%   info       porosities and counts
%
% Voluminous pores: every measured volume takes the unused 3-cell with the
% closest volume. Sheet-like pores: volumes are drawn from the measured
% cumulative distribution (rng(1)) until the target porosity is reached; each
% takes the smallest free 2-cell with area between (pi/4)(20 V/pi)^(2/3) and
% (pi/4)(50 V/pi)^(2/3), or the largest free 2-cell if there is none, and its
% aperture is V/A. 2-cells bounding voluminous pores are not used.

bc = fcmplx(1).bc;
V_box = max(bc(:,1)) * max(bc(:,2)) * max(bc(:,3));
nF = vcmplx(3).num(3).val;

grain = vcmplx(4).vvol;               % volumes of the 3-cells of M
area  = vcmplx(3).vvol;               % areas of the 2-cells of M

% ---- voluminous pores ---------------------------------------------------
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

% ---- sheet-like pores: volumes ------------------------------------------
b = area;
for c = vol_cells(:)'
    b(nonzeros(vcmplx(4).bndop(3).indx(c,:))) = 0;   % faces of voluminous pores
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

% ---- sheet-like pores: 2-cells ------------------------------------------
nS = numel(V_s);
sel_area = nan(nS,1);
sel_face = nan(nS,1);
n_fallback = 0;
n_failed = 0;
for i = 1:nS
    lo = (pi/4) * (20 * V_s(i) / pi)^(2/3);
    hi = (pi/4) * (50 * V_s(i) / pi)^(2/3);
    idx = find(isfinite(b) & b > lo & b < hi);
    if ~isempty(idx)
        [a, k] = min(b(idx));              % smallest admissible free face
    else
        idx = find(isfinite(b) & b > 0);
        if isempty(idx)
            warning('assign_pores: no free 2-cell for sheet-like pore %d', i);
            n_failed = n_failed + 1;
            continue
        end
        [a, k] = max(b(idx));              % fallback: largest free face
        n_fallback = n_fallback + 1;
    end
    sel_area(i) = a;
    sel_face(i) = idx(k);
    b(idx(k)) = -Inf;                      % each face is used once
end

h_face = zeros(nF,1);
for i = 1:nS
    if isnan(sel_face(i)) || isnan(sel_area(i)), continue; end
    h_face(sel_face(i)) = V_s(i) / sel_area(i);
end

info = struct('vf_vol', vf_vol, 'vf_sheet', vf_s, 'n_vol', numel(vol_cells), ...
              'n_sheet', nnz(h_face > 0), 'n_fallback', n_fallback, 'n_failed', n_failed);
fprintf('pores: %d voluminous (porosity %.4f), %d sheet-like (porosity %.4f, median aperture %.1f um), %d fallback faces\n', ...
    info.n_vol, vf_vol, info.n_sheet, vf_s, median(h_face(h_face > 0)), n_fallback);
end
