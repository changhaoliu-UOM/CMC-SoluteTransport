function cfg = config_material(mat)
% Material data and run settings.
%   mat  'SG' (sintered glass), 'BS' (Bentheimer sandstone) or
%        'SGL' (extended SG model, 4 x 4 x 10 mm)
% Lengths are in um, times in s and pressures in Pa.
switch upper(mat)
    case 'SG'
        cfg.pore_data = 'SG';         % prefix of the pore-data files in data/
        cfg.vf        = 0.278;        % target (resolved) porosity
        cfg.Tf_in     = [182 92];     % s, 10-90% rise time of the inlet history [low high]
        cfg.v_meas    = [37.1 75.5];  % um/s, measured arrival velocity [low high]
        cfg.dTf_meas  = [28 18];      % s, measured increase of the filling duration [low high]
        cfg.t_end     = 650;          % s
        cfg.void_copies = 1;
    case 'BS'
        cfg.pore_data = 'BS';
        cfg.vf        = 0.193;
        cfg.Tf_in     = [187 102];
        cfg.v_meas    = [52.1 96.9];
        cfg.dTf_meas  = [23 13];
        cfg.t_end     = 650;
        cfg.void_copies = 1;
    case 'SGL'
        cfg.pore_data = 'SG';         % pore data of SG
        cfg.vf        = 0.278;
        cfg.Tf_in     = [182 92];
        cfg.v_meas    = [37.1 75.5];
        cfg.dTf_meas  = [NaN NaN];    % measured on the 5 mm sample only
        cfg.t_end     = 1300;
        cfg.void_copies = 2;          % voluminous pores assigned twice: keeps the
                                      % voluminous-pore porosity of SG
    otherwise
        error('config_material: unknown material %s', mat);
end
cfg.material = upper(mat);
root = fileparts(fileparts(mfilename('fullpath')));   % this file is in <root>/src
cfg.root     = root;
cfg.data_dir = fullfile(root, 'data');

cfg.Q_low = 1.42e8;            % um^3/s (0.142 uL/s); higher flow rate = 2 Q_low
cfg.flows = {'low', 'high'};
cfg.Dm    = 1770;              % um^2/s, molecular diffusion coefficient
cfg.mu    = 1e-3;              % Pa s, dynamic viscosity
cfg.dt    = 0.1;               % s
cfg.z_fraction     = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];   % monitoring planes, z/L
cfg.slab_halfwidth = 0.02;          % half width of the monitoring slabs, fraction of L
end
