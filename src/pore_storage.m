function [S_pore, info] = pore_storage(fcmplx, vcmplx, vol_cells, h_face)
% Pore volume stored at every vertex (0-cell) of K, in um^3.
%   Voluminous pores (3-cells vol_cells of M): every cube of K inside the pore
%     gives 1/8 of its volume to each of its 8 vertices.
%   Sheet-like pores (2-cells of M with aperture h_face > 0): every quadrilateral
%     of K in the face gives 1/4 of (area x aperture) to each of its 4 vertices.
% The total equals the reconstructed pore volume.

n0 = fcmplx(1).num(1).val;
S_pore = zeros(n0,1);

% ---- voluminous pores ----------------------------------------------------
% The cubes of K are numbered polyhedron by polyhedron of M, one per vertex.
nP = vcmplx(4).num(4).val;
cube_owner = repelem((1:nP)', vcmplx(4).num(1).val(:));
if numel(cube_owner) ~= fcmplx(4).num(4).val
    error('pore_storage: cube count of K does not match the polyhedron-vertex count of M');
end
is_vol = false(nP,1);
is_vol(vol_cells(:)) = true;
kc = find(is_vol(cube_owner));
cube_nodes = abs(fcmplx(4).bndop(1).indx(kc,1:8));
if any(cube_nodes(:) == 0)
    error('pore_storage: a cube of K without 8 vertices');
end
if any(max(cube_nodes,[],2) - fcmplx(1).break(3) ~= cube_owner(kc))
    error('pore_storage: cube ownership check failed');
end
w_cube = fcmplx(4).bfvol(kc) / 8;
S_pore = S_pore + accumarray(cube_nodes(:), repmat(w_cube,8,1), [n0 1]);

% ---- sheet-like pores ----------------------------------------------------
% The first fcmplx(3).break(1) 2-cells of K lie in the faces of M, numbered
% face by face, one per vertex of the face.
nF = vcmplx(3).num(3).val;
nq = fcmplx(3).break(1);
quad_owner = repelem((1:nF)', vcmplx(3).num(1).val(:));
if numel(quad_owner) ~= nq
    error('pore_storage: in-face 2-cell count of K does not match the face-vertex count of M');
end
h = zeros(nF,1);
h(1:numel(h_face)) = h_face(:);
kq = find(h(quad_owner) > 0);
quad_nodes = abs(fcmplx(3).bndop(1).indx(kq,1:4));
if any(quad_nodes(:) == 0)
    error('pore_storage: an in-face 2-cell of K without 4 vertices');
end
if any(max(quad_nodes,[],2) - fcmplx(1).break(2) ~= quad_owner(kq))
    error('pore_storage: quadrilateral ownership check failed');
end
w_quad = fcmplx(3).bfvol(kq) .* h(quad_owner(kq)) / 4;
S_pore = S_pore + accumarray(quad_nodes(:), repmat(w_quad,4,1), [n0 1]);

% ---- checks ----------------------------------------------------------------
sheet_faces = find(h > 0);
bc = fcmplx(1).bc;
info.V_vol    = sum(fcmplx(4).bfvol(kc));
info.V_sheet  = sum(fcmplx(3).bfvol(kq) .* h(quad_owner(kq)));
info.V_total  = sum(S_pore);
info.porosity = info.V_total / prod(max(bc) - min(bc));
info.n_pore   = nnz(S_pore > 0);
fprintf('storage: voluminous %.4e um^3 (target %.4e), sheet-like %.4e um^3 (target %.4e), porosity %.4f\n', ...
    info.V_vol, sum(vcmplx(4).vvol(vol_cells)), info.V_sheet, ...
    sum(vcmplx(3).vvol(sheet_faces) .* h(sheet_faces)), info.porosity);
end
