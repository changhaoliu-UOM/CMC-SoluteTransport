function S = pore_storage(fcmplx, vcmplx, vol_cells, h_face)
% Pore storage volume of each node of K. Each 3-cell of K in a voluminous
% pore gives 1/8 of its volume to each of its 8 vertices; each 2-cell of K
% in a sheet-like pore gives 1/4 of (area x aperture) to each of its 4
% vertices.
n0 = fcmplx(1).num(1).val;
S  = zeros(n0,1);

nP = vcmplx(4).num(4).val;
owner3 = repelem((1:nP)', vcmplx(4).num(1).val(:));
is_vol = false(nP,1);
is_vol(vol_cells(:)) = true;
kc = find(is_vol(owner3));
nodes = abs(fcmplx(4).bndop(1).indx(kc,1:8));
S = S + accumarray(nodes(:), repmat(fcmplx(4).bfvol(kc) / 8, 8, 1), [n0 1]);

nF = vcmplx(3).num(3).val;
owner2 = repelem((1:nF)', vcmplx(3).num(1).val(:));
h = zeros(nF,1);
h(1:numel(h_face)) = h_face(:);
kq = find(h(owner2) > 0);
nodes = abs(fcmplx(3).bndop(1).indx(kq,1:4));
S = S + accumarray(nodes(:), repmat(fcmplx(3).bfvol(kq) .* h(owner2(kq)) / 4, 4, 1), [n0 1]);
end
