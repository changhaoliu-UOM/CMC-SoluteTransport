function [gH, gD, node_active, info] = pore_conductances(M, mu, D_pore)
% Hydraulic and diffusive conductances of the 1-cells of K by the composition
% model (Section 3.3 of the paper; cmc-py, model = "composition"), with the
% diagonal Hodge star.
%   M       struct with fcmplx (K), vcmplx (M), vol_cells and h_face (assign_pores.m)
%   mu      dynamic viscosity (Pa s);  D_pore  molecular diffusion coefficient (um^2/s)
%
%   gH, gD       (n1 x 1) hydraulic (um^3/(Pa s)) and diffusive (um^3/s) conductances
%   node_active  (n0 x 1) vertices of K with a positive diffusive pair value at
%                one of their corners (the monitored nodes)
%
% Every 1-cell e of K collects one pair (e, F) per cube of K at each of its
% ends: F is the face of the cube orthogonal to e at that end, of area mu(F),
% and the pair contributes kappa(e,F) mu(F) / (8 mu(e)) to the conductance of e.
% The pair value kappa is set by the pore features that the cross-section F
% meets, with voluminous pores in the role of grains (value D_g) and
% sheet-like pores of aperture h in the role of grain boundaries (value D_f);
% the solid has the value zero:
%   parallel pair (F = Q(s,g), s a 1-cell): every edge of F at the vertex that
%       lies in a sheet adds a strip of width h/2 with the value D_f; the rest of
%       F has the value D_g; the parts act in parallel (area-weighted mean);
%   exchange pair (e leaves a sheet into a 3-cell): the half layer h/2 and the
%       3-cell act in series over the half length of e (length-weighted
%       harmonic mean); a zero value blocks the pair;
%   junction pair (all other pairs): D_g.
% Values: hydraulic  D_g = k_b/mu in voluminous pores, k_b the area-weighted mean
%                    over the faces f of the pore of (A_f/pi)/8 (Hagen-Poiseuille);
%                    D_f = h^2/(12 mu) in sheet-like pores (parallel plates)
%         diffusive  D_g = D_f = D_pore
% There are no separate channels: sheet-like pores sharing an edge of M are
% connected through their layers along that edge.
fc = M.fcmplx;  vc = M.vcmplx;
b1 = fc(1).break(1);  b2 = fc(1).break(2);  b3 = fc(1).break(3);
n0 = fc(1).num(1).val;  n1 = fc(2).num(2).val;
n2 = fc(3).num(3).val;  n3 = fc(4).num(4).val;

% cell of M represented by every vertex of K: dimension and index
vdim = zeros(n0,1);
vdim(b1+1:b2) = 1;  vdim(b2+1:b3) = 2;  vdim(b3+1:n0) = 3;
off  = [0; b1; b2; b3];
vid  = (1:n0)' - off(vdim + 1);

% 1-cells: end vertices, carrier (cell of M containing the 1-cell)
EV = fc(2).bndop(1).indx(:,1:2);
[etdim, ie] = max(vdim(EV), [], 2);
etid = vid(EV(sub2ind([n1 2], (1:n1)', ie)));
mu1  = fc(2).bfvol(:);

% 2-cells Q(sigma,tau): vertices, edges, lowest and highest cell of M
FV = fc(3).bndop(1).indx(:,1:4);
FE = fc(3).bndop(2).indx(:,1:4);
fsdim = min(vdim(FV), [], 2);
[ftdim, jf] = max(vdim(FV), [], 2);
ftid = vid(FV(sub2ind([n2 4], (1:n2)', jf)));
mu2  = fc(3).bfvol(:);

% 3-cells Q(v,g): faces, edges, cell g of M
CF = fc(4).bndop(3).indx(:,1:6);
CE = fc(4).bndop(2).indx(:,1:12);
CV = fc(4).bndop(1).indx(:,1:8);
if any(CF(:) == 0) || any(CE(:) == 0) || any(CV(:) == 0)
    error('pore_conductances: a 3-cell of K without 6 faces, 12 edges and 8 vertices');
end
gv = max(CV .* (vdim(CV) == 3), [], 2);
cg = gv - b3;                                        % cell of M of every cube

% pore features
nC = vc(4).num(4).val;  nF = vc(3).num(3).val;
h = zeros(nF,1);  h(1:numel(M.h_face)) = M.h_face(:);           % sheet aperture (0: none)
isbody = false(nC,1);  isbody(M.vol_cells) = true;
kb = zeros(nC,1);
for c = M.vol_cells(:)'
    f  = nonzeros(vc(4).bndop(3).indx(c,:));
    Af = vc(3).vvol(f);  Af = Af(:);
    kb(c) = sum(Af .* (Af / pi) / 8) / sum(Af);
end
DgH = isbody(cg) .* kb(cg) / mu;     DgD = D_pore * double(isbody(cg));
DfH = h.^2 / (12 * mu);              DfD = D_pore * double(h > 0);

gH = zeros(n1,1);  gD = zeros(n1,1);  node_active = false(n0,1);
n_par = 0;  n_exch = 0;  n_junc = 0;  n_resc = 0;  n_exch_open = 0;
for j = 1:12
    e = CE(:,j);
    for k = 1:2
        a  = EV(e,k);  bb = EV(e,3-k);
        % the face of the cube at a orthogonal to e: contains a, not the other end
        F = zeros(n3,1);  nhit = zeros(n3,1);
        for m = 1:6
            Fm = CF(:,m);
            fv = FV(Fm,:);
            pick = any(fv == a, 2) & ~any(fv == bb, 2);
            F(pick) = Fm(pick);  nhit = nhit + pick;
        end
        if any(nhit ~= 1)
            error('pore_conductances: orthogonal face not unique in %d cubes', nnz(nhit ~= 1));
        end
        w = mu2(F) ./ (8 * mu1(e));

        par  = fsdim(F) == 1;                         % F = Q(s,g)
        exch = ~par & ftdim(F) == 2 & (a == b2 + ftid(F));   % F = Q(v,f), a = vertex of f
        junc = ~par & ~exch;
        kH = DgH;  kD = DgD;                          % junction pairs: grain value

        % parallel pairs: strips of the sheets along the edges of F at a
        Aq = mu2(F) / 4;
        SA = zeros(n3,1);  SAH = SA;  SAD = SA;
        for m = 1:4
            Phi = FE(F,m);
            inA = EV(Phi,1) == a | EV(Phi,2) == a;
            sh  = par & inA & etdim(Phi) == 2;
            fi  = etid(Phi(sh));
            hs  = h(fi);
            Ab  = zeros(n3,1);  Ab(sh) = 0.5 * hs .* mu1(Phi(sh)) / 2;
            dH  = zeros(n3,1);  dH(sh) = DfH(fi);
            dD  = zeros(n3,1);  dD(sh) = DfD(fi);
            SA  = SA + Ab;  SAH = SAH + Ab .* dH;  SAD = SAD + Ab .* dD;
        end
        tot = SA ./ Aq;
        sc  = ones(n3,1);  over = tot > 1;  sc(over) = 1 ./ tot(over);
        gfr = min(max(1 - tot .* sc, 0), 1);
        kH(par) = SAH(par) ./ Aq(par) .* sc(par) + gfr(par) .* DgH(par);
        kD(par) = SAD(par) ./ Aq(par) .* sc(par) + gfr(par) .* DgD(par);
        n_resc = n_resc + nnz(par & over);

        % exchange pairs: half layer and grain in series over the half 1-cell
        ix = find(exch);
        if ~isempty(ix)
            fi = ftid(F(ix));
            L  = mu1(e(ix)) / 2;
            ef = min(0.5 * h(fi), L) ./ L;  gf = 1 - ef;
            kH(ix) = series_value(ef, DfH(fi), gf, DgH(ix));
            kD(ix) = series_value(ef, DfD(fi), gf, DgD(ix));
            n_exch_open = n_exch_open + nnz(ef > 0 & gf == 0);
        end

        gH = gH + accumarray(e, kH .* w, [n1 1]);
        gD = gD + accumarray(e, kD .* w, [n1 1]);
        node_active(a(kD > 0)) = true;
        n_par = n_par + nnz(par);  n_exch = n_exch + nnz(exch);  n_junc = n_junc + nnz(junc);
    end
end
info = struct('kb', kb(M.vol_cells), 'n_pairs', [n_par n_exch n_junc], ...
              'n_rescaled', n_resc, 'n_exchange_through_sheet', n_exch_open, ...
              'n_pore_edges', nnz(gD > 0));
end

% -------------------------------------------------------------------------
function v = series_value(ef, Df, gf, Dg)
% value of a pair made of the fractions ef (value Df) and gf (value Dg) in
% series; a zero value with a positive fraction blocks the pair
v = zeros(size(ef));
blocked = (ef > 0 & Df <= 0) | (gf > 0 & Dg <= 0);
rs = zeros(size(ef));                                % resistance per unit length
t = ef > 0 & ~blocked;  rs(t) = rs(t) + ef(t) ./ Df(t);
t = gf > 0 & ~blocked;  rs(t) = rs(t) + gf(t) ./ Dg(t);
ok = ~blocked & rs > 0;
v(ok) = 1 ./ rs(ok);
end
