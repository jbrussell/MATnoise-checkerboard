% Trace travel times through synthetic phase velocity map creating synthetic
% CSmeasure dataset with mostly dummy variables.
%
clear;

setup_parameters_tomo;

fname_checker = [parameters.workingdir,'/checker.mat']; % path to checkerboard map

workingdir = parameters.workingdir;
Xspoutputpath = [workingdir,'/Xsp/'];
if ~exist(Xspoutputpath)
    mkdir(Xspoutputpath)
end
lalim=parameters.lalim;
lolim=parameters.lolim;
periods = parameters.periods;

% Read checkboard map
temp = load(fname_checker);
checker = temp.checker;
xi = checker(1).xi;
yi = checker(1).yi;
xnode = xi(1:end,1)';
ynode = yi(1,1:end);

% Trace rays for each station
ipair = 0;
for ista1 = 1:length(stalist)
    sta1 = stalist{ista1};
    for ista2 = 1:length(stalist)
        sta2 = stalist{ista2};
        if strcmp(sta1,sta2)
            continue
        end
        ipair = ipair + 1;
        disp([sta1,'_',sta2]);

        % Trace traveltime through checkerboard map
        lat1 = stlat(ista1);
        lon1 = stlon(ista1);
        lat2 = stlat(ista2);
        lon2 = stlon(ista2);
        r = distance(lat1,lon1,lat2,lon2,referenceEllipsoid('GRS80'))/1000;
        dr = deg2km(mean(diff(xnode)));
        Nr = floor(r/dr);
        [lat_way,lon_way] = gcwaypts(lat1,lon1,lat2,lon2,Nr);
        dtp = [];
        for ip = 1:length(periods)
            phv = checker(ip).phv;
            phv_path = interp2(yi,xi,phv,lon_way,lat_way);
            dtp(ip) = r ./ mean(phv_path(:));
        end


        xspinfo.sta1 = sta1;
        xspinfo.sta2 = sta2;
        xspinfo.lat1 = lat1;
        xspinfo.lon1 = lon1;
        xspinfo.lat2 = lat2;
        xspinfo.lon2 = lon2;
        xspinfo.r = r;
        xspinfo.tw = dtp;
        xspinfo.coherenum = 1;
        xspinfo.sumerr = 0;
        xspinfo.snr = 999;
        
        twloc = 1./periods*2*pi;

        save([Xspoutputpath,'/',sta1,'_',sta2,'_xsp.mat'],'xspinfo','twloc');
    end
end

