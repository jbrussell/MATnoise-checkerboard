% mk_stageom
% quick script to make different station configuations
%% set up
clear all;

nsta = 50;
distint = 50/1000; % in km

option = 1;
outpath = './';
outfile = sprintf('stations_%d.txt',option);

% region info
lat = [21 22];
lon = [-158 -157];
deplim = 0; % depth limit to avoid plotting crustal seismicity
latstart = 21.45;
lonstart = -157.45;

% load map grid to extract elevations
[elev,longs,lats]=m_etopo2([lon,lat]);
[a,b] = size(elev);
elvvec = reshape(elev,[1,a*b]);
latvec = reshape(lats,[1,a*b]);
lonvec = reshape(longs,[1,a*b]);

%% testing possible station locations

if option ==1 % square grid
    ninst = round(sqrt(nsta)); % width and height of grid
    az = 30;
    [latline1,lonline1] = reckon(latstart,lonstart,distint/111.1.*[0:ninst-1],az);
    az = az+90;
    for ii = 1:length(latline1)
        
    [latline(:,ii),lonline(:,ii)] = reckon(latline1(ii),lonline1(ii),distint/111.1.*[0:ninst-1],az);
    end
    [a,b] = size(latline);
    stalats = reshape(latline,[a*b,1]);
    stalons = reshape(lonline,[a*b,1]);
end


%% put together vectors for new stations

for ii = 1:length(stalats)
    [arclen,az] = distance(stalats(ii),stalons(ii),latvec,lonvec);
    [val,idx] = min(arclen);
    staelvs(ii) = elvvec(idx);
    staelvs(ii) = 0; %for now just set default to sealevel
end
stanmss = [1:length(stalats)];

%% organizing and writing text file

formatSpec = ('%s %4.5f %5.5f %6.5f\n');
fileID = fopen([outpath,outfile],'w');

for ii = 1:length(stalats)
    stalat = stalats(ii);
    stalon = stalons(ii);
    staelv = staelvs(ii);
    stanms = num2str(stanmss(ii));
    fprintf(fileID,formatSpec,stanms,stalat,stalon,staelv);
end


fclose(fileID);

