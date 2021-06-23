% Setup_parameters for raytomo
%
% https://github.com/jbrussell
addpath('./tomo_functions');
% parameters.qpath = './MINEOS_qfile/';

parameters.workingdir = './RESULTS/';

% Path to station file
parameters.station_list = './stations.txt'; % text file containing stations [staname, lat, lon, elev]
[stalist, stlat, stlon, stel] = textread(parameters.station_list,'%s %f %f %f');

parameters.periods = [5 10 20]; % periods of interest

% Lat Lon
parameters.lalim = [5 12] ;
parameters.lolim = [-150 -142];
parameters.gridsize = 0.25; % degrees?
parameters.agebins = [165:5:175];
parameters.bathybins = [-9000 :5000: 1000];
parameters.gridsize_azi = 0.25; %3; %1.5; % gridsize for 2D azimuthal anisotropy (degrees)
parameters.r = 0.05; %0.01; % controls color bar [avgv(1-r) avgv(1+r)]

% Smoothing parameters
parameters.smweight0 = [2 5 10]; % isotropic second derivative smoothing
parameters.smweight0_azi = 1e3 * ones(size(parameters.periods)); %1000; % anisotropic second derivative smoothing
parameters.flweight0_azi = 1000 * ones(size(parameters.periods)); %1000; % anisotropic first derivative flatness
parameters.damp0_azi = 1000 * ones(size(parameters.periods)); %1000; % anisotropic norm damping

% parameters for the tomography (QC)
parameters.raydensetol=deg2km(parameters.gridsize)*0.25; %deg2km(parameters.gridsize); %deg2km(parameters.gridsize)*2;
parameters.raydensetol_azi=deg2km(parameters.gridsize_azi)*0.25; %deg2km(parameters.gridsize)*2;
parameters.fiterrtol = 2;   % error allowed in the wavelet fitting
parameters.dterrtol = 4;    % largest variance of the inversion error allowed
parameters.maxerrweight = 5; % Maximum error weight
parameters.polyfit_dt_err = 2; % (s) dt error greater than this, weighted 0

if ~exist(parameters.workingdir)
    mkdir(parameters.workingdir)
end