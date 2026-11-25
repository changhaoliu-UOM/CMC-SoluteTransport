function [vcmplx, fcmplx, fBO_1, fs1, fs2] = build_complexes(tess_file)
% Physical complex M from a Neper tessellation, its Forman subdivision K,
% the boundary operator of the 1-cells of K and the Hodge stars s0 and s1.
names = {'vertex', 'edge', 'face', 'polyhedron'};

vcmplx = setup_complex(4);
for k = 1:4, vcmplx(k).name = names{k}; end
vcmplx = read_tess(tess_file, vcmplx);
vcmplx = vgeo(vcmplx);
vcmplx = vcon(vcmplx);

fcmplx = setup_complex(4);
for k = 1:4, fcmplx(k).name = names{k}; end
fcmplx = v2f(fcmplx, vcmplx);
fcmplx = fcon(fcmplx);
fBO_1  = build_BO(fcmplx, 1);
fcmplx = fgeo2026(fcmplx);

n0 = fcmplx(1).num(1).val;
n1 = fcmplx(2).num(2).val;
fs1 = spdiags(fcmplx(1).fvol, 0, n0, n0);
fs2 = spdiags(fcmplx(2).fvol, 0, n1, n1);
end
