% Setup_parameters for raytomo
%
% https://github.com/jbrussell
addpath('./tomo_functions');
% parameters.qpath = './MINEOS_qfile/';

parameters.workingdir = '../RESULTS/';

% Path to station file
parameters.station_list = './stations.txt'; % text file containing stations [staname, lat, lon, elev]
[stalist, stlat, stlon, stel] = textread(parameters.station_list,'%s %f %f %f');

parameters.periods = [5 10 20]; % periods of interest
parameters.frange = [1/max(parameters.periods) 1/min(parameters.periods)]; % [Hz]
parameters.per_ind = [1:3]; % index of periods to consider

% Lat Lon
parameters.lalim = [5 12] ;
parameters.lolim = [-150 -142];
parameters.gridsize = 0.1; %0.25; % degrees?
parameters.agebins = [165:5:175];
parameters.bathybins = [-9000 :5000: 1000];
parameters.gridsize_azi = 0.1; %0.25; %3; %1.5; % gridsize for 2D azimuthal anisotropy (degrees)
parameters.r = 0.05; %0.01; % controls color bar [avgv(1-r) avgv(1+r)]

% QC parameters
parameters.snr_tol = 3; % minimum signal-to-noise
parameters.r_tol_min = 0; % [km] minimum station separation
parameters.r_tol_max = 600; % [km] maximum station separation
parameters.err_tol = inf; % maximum misfit of bessel fit between observed and synthetic
parameters.is_rtolmin_wavelength = 1; parameters.wl_fac = 1.0; % determine distance tolerance by wavelength?
parameters.dep_tol = [0 0]; % [sta1 sta2] OBS depth tolerance
parameters.is_raydensity_thresh = 0; % Apply raydensity threshold to wipe out poorly constrained grid cells?
parameters.min_dep = 9999; %-3500 for min station depth to use

% Kernel Parameters
parameters.kernel_path = ['./SEM2D_FFK_save/']; % path to saved kernels
parameters.kmode = -1; % <0 : bandwidth; =0 instantaneous; >0 Fresnel zones
parameters.nfreqs = 30; % number of frequencies to include in bandwidth average
parameters.bw_frac = 0.25; % fraction of frequency to include in bandwidht
parameters.tphase_in = []; % Path to SEM traveltime surfaces. If left blank, will use analytical kernels by default.

% Smoothing parameters
parameters.damp0 = 1e-4; % Norm damping of isotropic phase velocity
parameters.smweight0 = 10; % isotropic second derivative smoothing
parameters.smweight0_azi = 1e3; %1000; % anisotropic second derivative smoothing
parameters.flweight0_azi = 0; %1000; % anisotropic first derivative flatness
parameters.damp0_azi = 1e3; %0; % anisotropic norm damping (damp to zero)
parameters.is_wlsmooth = 0; %1; % weight smoothing by wavelength 0:no, 1:yes

% parameters for the tomography (QC)
parameters.raydensetol=deg2km(parameters.gridsize)*0.25; %deg2km(parameters.gridsize); %deg2km(parameters.gridsize)*2;
parameters.raydensetol_azi=deg2km(parameters.gridsize_azi)*0.25; %deg2km(parameters.gridsize)*2;
parameters.fiterrtol = 2;   % error allowed in the wavelet fitting
parameters.dterrtol = 4;    % largest variance of the inversion error allowed
parameters.maxerrweight = 5; % Maximum error weight
parameters.polyfit_dt_err = 2; % (s) dt error greater than this, weighted 0
parameters.stderrfac = 2; % remove measurements with err(i) > std(err)*stderrfac

if ~exist(parameters.workingdir)
    mkdir(parameters.workingdir)
end