% Script to do the ray theory tomography based on the ambient noise measurement. 
% Solves for azimuthal anisotropy and isotropic velocity on a 2D grid.
%
% This version uses a wavelength dependent smoothing factor
%
% phv(theta,freq) = phv_iso(freq) + Ac2(freq)*cos(2*theta) + As2(freq)*sin(2*theta)
%                                 + Ac4(freq)*cos(4*theta) + As4(freq)*sin(4*theta)
%
% Written by Ge Jin, jinwar@gmail.com
% Nov 2012
% JBR: Modified in 2018 to include azimuthal anisotropy
% https://github.com/jbrussell
clear; close all;

setup_parameters_tomo;

workingdir = parameters.workingdir;
Xsp_path = [workingdir,'/Xsp/'];
periods = parameters.periods;

%%
%======================= PARAMETERS =======================%

% Save results?
isoutput = 1;
savefile = [workingdir,'/raytomo.mat'];

frange = [1/max(periods) 1/min(periods)]; % [Hz]

average_vel = 3.8; % [km/s] For calculating wavelength for determining r_tol_min

% QC parameters
snr_tol = 0; % minimum signal-to-noise
is_rtolmin_wavelength = 1; wl_fac = 1.0; % determine distance tolerance by wavelength?
r_tol_min = 0; % [km] minimum station separation
r_tol_max = 600; % [km] maximum station separation
err_tol = 999; % maximum misfit of bessel fit between observed and synthetic

fastdir = 78; % Fast direction for azimuthal anisotropy (only for plotting purposes);
iscompare_aniso = 0; % compare to old anisotropic measurements
comp = {'ZZ'};
%==========================================================%
%%
% Set up geometry parameters
setup_parameters_tomo;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
gridsize_azi = parameters.gridsize_azi;
station_list = parameters.station_list;

% Load station info
[sta.name, sta.lat, sta.lon, sta.dep] = textread(station_list,'%s %f %f %f');

fiterrtol = parameters.fiterrtol;
maxerrweight = parameters.maxerrweight;
polyfit_dt_err = parameters.polyfit_dt_err;
is_wl_smoothing = parameters.is_wl_smoothing; % wavelength dependent smoothing
smweight0 = parameters.smweight0;
smweight0_azi = parameters.smweight0_azi;
flweight0_azi = parameters.flweight0_azi;
damp0_azi = parameters.damp0_azi;
dterrtol = parameters.dterrtol;
raydensetol = parameters.raydensetol;
raydensetol_azi = parameters.raydensetol_azi;
r = parameters.r;

xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
[xi yi] = ndgrid(xnode,ynode);
Nx = length(xnode);
Ny = length(ynode);

xnode_azi=lalim(1):gridsize_azi:lalim(2);
ynode_azi=lolim(1):gridsize_azi:lolim(2);
[xi_azi yi_azi] = ndgrid(xnode_azi,ynode_azi);
Nx_azi = length(xnode_azi);
Ny_azi = length(ynode_azi);

% figure output path
phv_fig_path = ['./figs/'];
if ~exist(phv_fig_path)    
    mkdir(phv_fig_path);
end

% read in bad station list, if existed
if exist('badsta.lst')
    badstnms = textread('badsta.lst','%s');
    badstaids = find(ismember({stainfo.staname},badstnms));
    disp('Found Bad stations:')
    disp(badstnms)
end

%% Set up initial smoothing kernels (second derivative)
% Isotropic smoothing kernels
F_iso = smooth_kernel_build(xnode, ynode, Nx*Ny);
F = sparse(Nx*Ny*2,Nx*Ny+Nx_azi*Ny_azi*2);
F(1:end,1:Nx*Ny) = F_iso;
% Azimuthal smoothing kernels
F_azi = smooth_kernel_build(xnode_azi, ynode_azi, Nx_azi*Ny_azi);
F_azic = sparse(Nx_azi*Ny_azi*2,Nx*Ny+Nx_azi*Ny_azi*2);
F_azic(1:end,Nx*Ny+[1:Nx_azi*Ny_azi]) = F_azi;
F_azis = sparse(Nx_azi*Ny_azi*2,Nx*Ny+Nx_azi*Ny_azi*2);
F_azis(1:end,Nx*Ny+Nx_azi*Ny_azi+[1:Nx_azi*Ny_azi]) = F_azi;

%% Set up initial flattening kernels (first derivative)
% Azimuthal flattening kernels
J_azi = flat_kernel_build(xnode_azi, ynode_azi, Nx_azi*Ny_azi);
J_azic = sparse(Nx_azi*Ny_azi*2,Nx*Ny+Nx_azi*Ny_azi*2);
J_azic(1:end,Nx*Ny+[1:Nx_azi*Ny_azi]) = J_azi;
J_azis = sparse(Nx_azi*Ny_azi*2,Nx*Ny+Nx_azi*Ny_azi*2);
J_azis(1:end,Nx*Ny+Nx_azi*Ny_azi+[1:Nx_azi*Ny_azi]) = J_azi;

%% Anisotropy norm damping
Nxyazi = Nx_azi*Ny_azi*2;
Nxy = Nx*Ny;
Areg_azi = zeros(Nxyazi,Nxy+Nxyazi);
azipart = eye(Nxyazi,Nxyazi); %.* diag(damp_azi);
Areg_azi(1:Nxyazi,Nxy+1:Nxy+Nxyazi) = azipart;
F_azi_damp = Areg_azi;

%%
% Initialize the xsp structure
xspfiles = dir([Xsp_path,'*_xsp.mat']);

disp('Looking at Xsp Files')
for ixsp = 1:length(xspfiles)
    
    temp = load([Xsp_path,xspfiles(ixsp).name]);
    xspinfo = temp.xspinfo;
    
    if ixsp ==1
        Tperiods = (2*pi)./temp.twloc;
%         waxis = temp.waxis;
%         twloc = temp.twloc;
        xspinfo.isgood = zeros(size(Tperiods));
        xspsum = xspinfo;
        wavelength = average_vel*Tperiods;
    else
        xspinfo.isgood = zeros(size(Tperiods));
        xspsum = [xspsum;xspinfo];
    end
    clear temp

    
    % 	xspinfo(ixsp).isgood = 0;
%     if xspsum(ixsp).sumerr < errlevel ...
%             && xspsum(ixsp).snr > snrtol && xspsum(ixsp).coherenum > mincoherenum
%         xspsum(ixsp).isgood = 1;
%     end

    for ip = 1:length(Tperiods)
        if ~is_rtolmin_wavelength && xspinfo.snr >= snr_tol && xspinfo.r >= r_tol_min && xspinfo.r <= r_tol_max && xspinfo.sumerr <= err_tol
            xspsum(ixsp).isgood(ip) = 1;
        elseif  is_rtolmin_wavelength && xspinfo.snr >= snr_tol && xspinfo.r >= wavelength(ip)*wl_fac && xspinfo.r <= r_tol_max && xspinfo.sumerr <= err_tol
            xspsum(ixsp).isgood(ip) = 1;
        end

        if isempty(find(strcmp(sta.name,xspsum(ixsp).sta1))) || isempty(find(strcmp(sta.name,xspsum(ixsp).sta2)))
            xspsum(ixsp).isgood(ip) = 0;
        end
    end
    
    if rem(ixsp,500)==0
        disp(['Looking at #',num2str(ixsp),' of ',num2str(length(xspfiles))])
    end
end % end of loop ixsp'


% Loop through periods
for ip=1:length(Tperiods)
    disp(' ');
    disp(['Inverting Period: ',num2str(Tperiods(ip))]);
    clear rays dt fiterr mat phaseg err raydense dist azi mat_azi phv
    raynum = 0;

    for ixsp = 1:length(xspsum)
        if xspsum(ixsp).isgood(ip) ==0
            continue;
        end
%         if xspsum(ixsp).r > refv*Tperiods(ip)*distrange(2)...
%                 || xspsum(ixsp).r < refv*Tperiods(ip)*distrange(1)
%             continue;
%         end
        
        raynum = raynum+1;
        rays(raynum,1) = xspsum(ixsp).lat1;
        rays(raynum,2) = xspsum(ixsp).lon1;
        rays(raynum,3) = xspsum(ixsp).lat2;
        rays(raynum,4) = xspsum(ixsp).lon2;
        
        % dist(raynum) = deg2km(distance(rays(raynum,1),rays(raynum,2),rays(raynum,3),rays(raynum,4)));
        dist(raynum) = distance(rays(raynum,1),rays(raynum,2),rays(raynum,3),rays(raynum,4),referenceEllipsoid('GRS80'))/1000;
        dt(raynum) = xspsum(ixsp).tw(ip);
        phv(raynum) = dist(raynum)./dt(raynum);
        
        dep1 = sta.dep(strcmp(xspsum(raynum).sta1,sta.name));
        dep2 = sta.dep(strcmp(xspsum(raynum).sta2,sta.name));
        dep(raynum) = mean([dep1 dep2]);
        
        %JRB - load azimuthal anisotropy
        [~,azi(raynum)]=distance(xspsum(ixsp).lat1,xspsum(ixsp).lon1,xspsum(ixsp).lat2,xspsum(ixsp).lon2);
        if azi(raynum) > 180
            azi(raynum) = azi(raynum) - 360;
        end
        if iscompare_aniso
            A2 = aniso.A2(ip);
            A4 = aniso.A4(ip);
            phi2 = aniso.phi2(ip);
            phi4 = aniso.phi4(ip);
            c_iso = aniso.c_iso(ip);       
        end
%         if comp{1}(1) == 'Z'
%             phv_cor = dist(raynum)./dt(raynum) - A2*c_iso*cosd(2*(azi - phi2));
%         elseif comp{1}(1) == 'T'
%             phv_cor = dist(raynum)./dt(raynum) - A2*c_iso*cosd(2*(azi-phi2)) - A4*c_iso*cosd(4*(azi-phi4)); 
%         end
%         dt(raynum) = dist(raynum)./phv_cor;
        
%         err = smooth((abs(xspsum(ixsp).err)./mean(abs(xspsum(ixsp).xsp))).^2,round(length(waxis)/length(twloc)));
%         fiterr(raynum) = interp1(waxis(:),err(:),twloc(ip)); 
        % Fix the fact that last period always breaks (JBR 9/29/17)
%         if isnan(fiterr(raynum))
%             [~,I] = min(abs(twloc(ip)-waxis(:)));
%             fiterr(raynum) = err(I);
%         end
        fiterr(raynum) = 1;
        csnum(raynum) = xspsum(ixsp).coherenum;
        snr(raynum) = xspsum(ixsp).snr;

        
        % JBR - Build azimuthal part of data kernel
%         mat_azi(raynum,:) = dist(raynum) * [cosd(2*azi(raynum)), sind(2*azi(raynum)), cosd(4*azi(raynum)), sind(4*azi(raynum)) ];
   
    end
    if size(dt,1) ~=raynum
        dt = dt';
    end
    
    % Building the isotropic data kernel
    disp('Start building the kernel');
    tic
    mat_iso=ray_kernel_build(rays,xnode,ynode);   
    toc
    [mat_azi, mat_azi_hits] = ray_kernel_build_azi(rays,xnode_azi,ynode_azi);
    
    % JBR - Combine isotropic and anisotropic
    mat = [mat_iso, mat_azi];
    
    % Calculate the weighting matrix
    W = sparse(length(dt),length(dt));
    for i=1:length(dt)
        W(i,i)=1./fiterr(i);
    end
    ind = find(W > maxerrweight);
    W(ind) = maxerrweight;
    ind = find(W < 1/fiterrtol);
    W(ind) = 0;
    for i=1:length(dt)
        W(i,i)=W(i,i).*(csnum(i).^0.5);
    end
    para = polyfit(dist(:),dt,1);
    polyerr = polyval(para,dist(:)) - dt;
    errind = find(abs(polyerr) > polyfit_dt_err);
    for i = errind
        W(i,i) = 0;
    end
    
    if is_wl_smoothing
        wavelength = nanmean(phv) .* Tperiods(ip);
        wl_smooth_factor = wavelength ./ deg2km(gridsize);
    else
        wl_smooth_factor = 1;
    end
    
    % calculate the smoothing weight
    NR=norm(F,1);
    NA=norm(W*mat,1);
    smweight = smweight0(ip)*NA/NR * wl_smooth_factor;
    
    NR=norm(F_azi_damp,1);
    NA=norm(W*mat,1);
    damp_azi = damp0_azi(ip)*NA/NR;
    
    NR=norm(F_azic,1);
    NA=norm(W*mat,1);
    smweight_azi = smweight0_azi(ip)*NA/NR * wl_smooth_factor;
    
    NR=norm(J_azic,1);
    NA=norm(W*mat,1);
    flweight_azi = flweight0_azi(ip)*NA/NR * wl_smooth_factor;
    
    disp('start inverse');
    A=[W*mat; smweight*F; damp_azi*F_azi_damp; smweight_azi*F_azic; smweight_azi*F_azis; flweight_azi*J_azic; flweight_azi*J_azis];
    rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azi_damp,1),1); zeros(size(F_azic,1),1); zeros(size(F_azis,1),1); zeros(size(J_azic,1),1); zeros(size(J_azis,1),1)];

    phaseg=(A'*A)\(A'*rhs);
    %        toc
    %        disp('Done');
    
    
    % Iteratively down weight the measurement with high error
    niter=1;
    
    while niter < 2
        niter=niter+1;
        err = mat*phaseg - dt;

        stderr=std(err);
        if stderr > dterrtol
            stderr = dterrtol;
        end
        ind = find(diag(W)==0);
        disp('Before iter:');
        disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
        disp(['Bad Measurement Number: ', num2str(length(ind))]);
        for i=1:length(err)
            if abs(err(i)) > 2*stderr
                W(i,i)=0;
            end
        end
        ind = find(diag(W)==0);
        disp('After iter:');
        disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
        disp(['Bad Measurement Number: ', num2str(length(ind))]);
        
        if is_wl_smoothing
            wavelength = 1./nanmean(phaseg(1:Nx*Ny)) .* Tperiods(ip);
            wl_smooth_factor = wavelength ./ deg2km(gridsize);
        else
            wl_smooth_factor = 1;
        end
        
        % Rescale the smooth kernel
        NR=norm(F,1);
        NA=norm(W*mat,1);
        smweight = smweight0(ip)*NA/NR * wl_smooth_factor;
        
        NR=norm(F_azi_damp,1);
        NA=norm(W*mat,1);
        damp_azi = damp0_azi(ip)*NA/NR;
        
        NR=norm(F_azic,1);
        NA=norm(W*mat,1);
        smweight_azi = smweight0_azi(ip)*NA/NR * wl_smooth_factor;
        
        NR=norm(J_azic,1);
        NA=norm(W*mat,1);
        flweight_azi = flweight0_azi(ip)*NA/NR * wl_smooth_factor;
        
        % Invert
        A=[W*mat; smweight*F; damp_azi*F_azi_damp; smweight_azi*F_azic; smweight_azi*F_azis; flweight_azi*J_azic; flweight_azi*J_azis];
        rhs=[W*dt; zeros(size(F,1),1); zeros(size(F_azi_damp,1),1); zeros(size(F_azic,1),1); zeros(size(F_azis,1),1); zeros(size(J_azic,1),1); zeros(size(J_azis,1),1)];

        phaseg=(A'*A)\(A'*rhs);

        
    end
    
    % Anisotropic terms from model vector
    phaseg_azic = phaseg(Nx*Ny+1 : Nx*Ny+Nx_azi*Ny_azi);
    phaseg_azis = phaseg(Nx*Ny+Nx_azi*Ny_azi+1 : Nx*Ny+Nx_azi*Ny_azi*2);
    
    % Isotropic phase velocity
    phv_iso = dist'./(mat_iso*phaseg(1:Nx*Ny));
    
    %        disp(' Get rid of uncertainty area');
    
    Igood = find(diag(W)~=0);
    mat_good = mat(Igood,:);
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            %raydense(i,j) = sum(mat(:,n));
            raydense(i,j) = sum(mat_good(:,n));
            if raydense(i,j) < raydensetol
                phaseg(n)=NaN;
            end
        end
    end
    for i=1:Nx_azi
        for j=1:Ny_azi 
            n=Ny_azi*(i-1)+j;
            raydense_azi(i,j) = sum(mat_azi_hits(:,n));
            if raydense_azi(i,j) < raydensetol_azi
                phaseg_azic(n)=NaN;
                phaseg_azis(n)=NaN;
            end
        end
    end
    
    % Convert into phase velocity
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            GV(i,j)= 1./phaseg(n);
        end
    end
    for i=1:Nx_azi
        for j=1:Ny_azi
            n=Ny_azi*(i-1)+j;
            Ac(i,j) = -phaseg_azic(n);
            As(i,j) = -phaseg_azis(n);
        end
    end
    
    % JBR - Get Azimuthal coefficients from phaseg (s/km)
    phv_av = nanmean(GV(:));
%     phv_av = nanmedian(GV(:));
    phv_avstd = nanstd(GV(:));
    slow_av = 1/phv_av;
%     slow_av = 1/mean(phv_iso);
    Ac2 = Ac./slow_av;
    As2 = As./slow_av;
    A2 = sqrt(Ac2.^2+As2.^2);
    phi2 = 1/2*atan2d(As2,Ac2);

    raytomo(ip).GV = GV;
    raytomo(ip).mat = mat;
    raytomo(ip).raydense = raydense;
    raytomo(ip).period = Tperiods(ip);
    raytomo(ip).w = diag(W);
    raytomo(ip).err = err;
    raytomo(ip).rays = rays;
    raytomo(ip).fiterr = fiterr;
    raytomo(ip).dt = dt;
    raytomo(ip).smweight0 = smweight0(ip) * wl_smooth_factor;
    raytomo(ip).smweight0_azi = smweight0_azi(ip) * wl_smooth_factor;
    raytomo(ip).flweight0_azi = flweight0_azi(ip) * wl_smooth_factor;
    %JBR    
    raytomo(ip).phv_iso = phv_iso;    
    raytomo(ip).phv_av = phv_av;
    raytomo(ip).phv_avstd = phv_avstd;
    raytomo(ip).Ac2 = Ac2;
    raytomo(ip).As2 = As2;
%     raytomo(ip).Ac4 = Ac4;
%     raytomo(ip).As4 = As4;
    raytomo(ip).A2 = A2;
%     raytomo(ip).A4 = A4;
    raytomo(ip).phi2 = phi2;
%     raytomo(ip).phi4 = phi4;
    raytomo(ip).phv = phv;
    raytomo(ip).azi = azi;
    raytomo(ip).dist = dist;
    
    
    if 0
        figure(1)
        clf
        ax = worldmap(lalim, lolim);
        set(ax, 'Visible', 'off')
        surfacem(xi,yi,raytomo(ip).GV);
%         drawlocal
        title([num2str(Tperiods(ip))],'fontsize',15)
        avgv = nanmean(raytomo(ip).GV(:));
        caxis([avgv*(1-r) avgv*(1+r)])
        colorbar
        colormap(tomo_cmap(100))
        
%         pause;
    end
    
end % end of period loop

lalim = [min(xnode) max(xnode)];
lolim = [min(ynode) max(ynode)];
[xi yi] = ndgrid(xnode,ynode);
% isoutput = 1;
if isoutput
    save(savefile,'raytomo','xnode','ynode');
    save('coor.mat','xi','yi','xnode','ynode','gridsize','lalim','lolim');
end

%% Azimuthal anisotropy (%)

% Mp = 3; Np = 3;
% fig16 = figure(16);
% set(gcf,'position',[1    1   1244   704]);
% clf
% vperc = [-r r];
% for ip=1:length(Tperiods)
%     subplot(Mp,Np,ip)
%     ax = worldmap(lalim+[-0.001 +0.001], lolim+[-0.001 +0.001]);
%     setm(gca,'MapProjection','mercator','FLineWidth',1.5,'FontSize',13)
%     tightmap
%     set(ax, 'Visible', 'on')
%     set(gcf,'color','w')
%     setm(gca,'FFaceColor',[0.9 0.9 0.9])
%     A2 = raytomo(ip).A2;
% %     surfacem(xi,yi,resid);
%     levels = linspace(0,0.03,10)*100;
%     surfacem(xi_azi,yi_azi,A2*100,'Linestyle','none');
% %     drawlocal
%     title([num2str(Tperiods(ip))],'fontsize',15)
%     caxis([min(levels) max(levels)])
%     colorbar
% %     colormap(seiscmap)
%     rbc = flip(redbluecmap);
% %     rbc = rbc([1 2 3 4 5 7 8 9 10 11],:);
% %     colormap(rbc);
%     colormap('parula');
%     
%     u=raytomo(ip).A2 .* cosd(raytomo(ip).phi2)*20;
% 	v=raytomo(ip).A2 .* sind(raytomo(ip).phi2)*20./cosd(mean(lalim));
% 	[m n]=size(xi_azi);
%     hold on;
%     xpts = [];
%     ypts = [];
% 	for ix=1:m
% 		for iy=1:n
% % 			if avgphv_aniso(ip).aniso_azi_std(ix,iy) < 40 && avgphv_aniso(ip).aniso_strength(ix,iy)>0.02
% %                 [xi_azi(ix,iy)-u(ix,iy)/2 xi_azi(ix,iy)+u(ix,iy)/2]
% %                 [yi_azi(ix,iy)-v(ix,iy)/2 yi_azi(ix,iy)+v(ix,iy)/2]
% % % 			geoshow([xi_azi(ix,iy)-u(ix,iy)/2 xi_azi(ix,iy)+u(ix,iy)/2]+gridsize_azi/2,...
% % % 					[yi_azi(ix,iy)-v(ix,iy)/2 yi_azi(ix,iy)+v(ix,iy)/2]+gridsize_azi/2,'Color','k','linewidth',2);
% %             plotm([yi_azi(ix,iy)-v(ix,iy)/2 yi_azi(ix,iy)+v(ix,iy)/2],...
% % 					[xi_azi(ix,iy)-u(ix,iy)/2 xi_azi(ix,iy)+u(ix,iy)/2],'k-','linewidth',2);
% % 			end
%             xpts = [xpts, [xi_azi(ix,iy)-u(ix,iy)/2 xi_azi(ix,iy)+u(ix,iy)/2]+gridsize_azi/2, nan];
%             ypts = [ypts, [yi_azi(ix,iy)-v(ix,iy)/2 yi_azi(ix,iy)+v(ix,iy)/2]+gridsize_azi/2, nan];
%         end
% 	end
%     plotm(xpts,ypts,'Color','k','linewidth',2);
%     hold on;
%     plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0 0 0]);
% end

%% Phase Velocity Maps (km/s)

Mp = 3; Np = 3;
fig17 = figure(17);
set(gcf,'position',[1    1   1244   704]);
clf
colormap(tomo_cmap(100));
for ip=1:length(Tperiods)
    subplot(Mp,Np,ip)
    ax = worldmap(lalim+[-0.001 +0.001], lolim+[-0.001 +0.001]);
    setm(gca,'MapProjection','mercator','FLineWidth',1.5,'FontSize',13)
    tightmap
    set(ax, 'Visible', 'on')
    set(gcf,'color','w')
    setm(gca,'FFaceColor',[0.9 0.9 0.9])
    avgv = nanmean(raytomo(ip).GV(:));
    levels = linspace(avgv*(1-r), avgv*(1+r),100);
    contourfm(xi,yi,raytomo(ip).GV,levels,'LineStyle','none');
    title([num2str(Tperiods(ip)),' s'],'fontsize',15)
    caxis([avgv*(1-r) avgv*(1+r)])
    cb = colorbar;
    ylabel(cb,'Phase Velocity','fontsize',15);
    set(cb,'linewidth',1.5,'fontsize',15);
    set(gca,'fontsize',15);
    hold on;
%     plotm(stlat,stlon,'ok','markerfacecolor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
end

%%
% RAY DENSITY
fig18 = figure(18);
set(gcf,'position',[1    1   1244   704]);
clf

for ip=1:length(Tperiods)
subplot(Mp,Np,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    surfacem(xi,yi,raytomo(ip).raydense);
%     drawlocal
    title([num2str(Tperiods(ip))],'fontsize',15)
    colorbar
    colormap(flip(hot));
    caxis([0 500])
end


%% PLOT AZIMUTHAL MODEL

fig4 = figure(4); clf;
% set(gcf,'position',[4         325        1239         380],'color','w');
set(gcf,'position',[6   220   490   485]);
periods = Tperiods;
% Old Values
if iscompare_aniso
    wRMS_2A = aniso.wRMS_2A;
    wRMS_4A = aniso.wRMS_4A;
    err_phi2 = aniso.err_phi2;
    err_phi4 = aniso.err_phi4;
    A2_2 = aniso.A2;
    A4_2 = aniso.A4;
    phi2_2 = aniso.phi2;
    phi4_2 = aniso.phi4;
end
for iper = 1:length(periods)
    % New values
    A2_rt(iper) = nanmean(raytomo(iper).A2(:));
    phi2_rt(iper) = nanmean(raytomo(iper).phi2(:));
    phv_av_rt(iper) = raytomo(iper).phv_av;
    phv_avstd_rt(iper) = raytomo(iper).phv_avstd;
end


% peak-to-peak
subplot(2,1,1); hold on;
if iscompare_aniso
    h3(2) = errorbar(periods,A4_2*2*100,wRMS_4A*100,'--ob','linewidth',2);
    h3(1) = errorbar(periods,A2_2*2*100,wRMS_2A*100,'--o','color',[0 0.7 0],'linewidth',2);
end
% h3(2) = plot(periods,A4_rt*2*100,'-ob','linewidth',2);
h3(1) = plot(periods,A2_rt*2*100,'-o','color',[0 0.7 0],'linewidth',2);
xlim(flip(1./frange));
ylim([0 5]);
set(gca,'linewidth',1.5,'xminortick','on','yminortick','on','fontsize',18);
xlabel('Period (s)','fontsize',18);
ylabel('Peak-to-peak amp (%)','fontsize',18);
legend(h3,{'2\theta','4\theta'},'location','northwest');

% Azimuth
subplot(2,1,2); hold on;
for iper = 1:length(periods)
    % OLD AZIMUTHAL PARAMTERS
    if iscompare_aniso
        phi2_vec(1) = phi2_2(iper);
        phi2_vec(2) = phi2_2(iper)+180;
        phi2_vec(3) = phi2_2(iper)-180;
        if comp{1}(1) == 'Z'
            [dif, I] = min(abs(phi2_vec-fastdir));
            phi2_2(iper) = phi2_vec(I);
            if phi2_2(iper) < fastdir && dif > 10
                phi2_2(iper) = phi2_2(iper)+180;
            end
        elseif comp{1}(1) == 'T'
            [dif, I] = min(abs(phi2_vec-fastdir+90));
            phi2_2(iper) = phi2_vec(I);
            if phi2_2(iper) < fastdir && dif
                phi2_2(iper) = phi2_2(iper)+180;
            end
        end


        phi4_vec(1) = phi4_2(iper);
        phi4_vec(2) = phi4_2(iper)+90;
        phi4_vec(3) = phi4_2(iper)+180;
        phi4_vec(4) = phi4_2(iper)+270;
        phi4_vec(5) = phi4_2(iper)-90;
        phi4_vec(6) = phi4_2(iper)-180;
        phi4_vec(7) = phi4_2(iper)-270;
        [~, I] = min(abs(phi4_vec-fastdir+45));
        phi4_2(iper) = phi4_vec(I);
        if phi4_2(iper) < fastdir
            phi4_2(iper) = phi4_2(iper)+90;
        end
    end
    
    % NEW AZIMUTHAL PARAMTERS
    phi2_vec(1) = phi2_rt(iper);
    phi2_vec(2) = phi2_rt(iper)+180;
    phi2_vec(3) = phi2_rt(iper)-180;
    if comp{1}(1) == 'Z'
        [dif, I] = min(abs(phi2_vec-fastdir));
        phi2_rt(iper) = phi2_vec(I);
        if phi2_rt(iper) < fastdir && dif > 10
            phi2_rt(iper) = phi2_rt(iper)+180;
        end
    elseif comp{1}(1) == 'T'
        [dif, I] = min(abs(phi2_vec-fastdir+90));
        phi2_rt(iper) = phi2_vec(I);
        if phi2_rt(iper) < fastdir
            phi2_rt(iper) = phi2_rt(iper)+180;
        end
    end
    phi2_rt(phi2_rt>160) = phi2_rt(phi2_rt>160)-180;
end
plot(periods,ones(size(periods))*fastdir,'--k','linewidth',2);
plot(periods,ones(size(periods))*fastdir-90,'--','color',[0.5 0.5 0.5],'linewidth',2);
plot(periods,ones(size(periods))*fastdir+90,'--','color',[0.5 0.5 0.5],'linewidth',2);
if iscompare_aniso
    errorbar(periods,phi4_2,err_phi4,'--ob','linewidth',2);
    errorbar(periods,phi2_2,err_phi2,'--o','color',[0 0.7 0],'linewidth',2);
end
% plot(periods,phi4_rt,'-ob','linewidth',2);
plot(periods,phi2_rt,'-o','color',[0 0.7 0],'linewidth',2);
ylabel('Fast Direction (%)','fontsize',18);
ylim([fastdir-130 fastdir+130]);
xlim(flip(1./frange));
set(gca,'linewidth',1.5,'xminortick','on','yminortick','on','fontsize',18);
xlabel('Period (s)','fontsize',18);

%% Plot phase velocities
for ip = 1:length(Tperiods)
    avgv(ip) = nanmean(raytomo(ip).GV(:));
    avgv_std(ip) = nanstd(raytomo(ip).GV(:));
end

fig2 = figure(2);clf;
set(gcf,'position',[320     2   508   703]);

subplot(2,1,1);
hold on; box on;
if iscompare_aniso
    errorbar(Tperiods,aniso.c_iso,aniso.err_c_iso*2,'-k','linewidth',2);
end
try 
    plot(Tperiods,xspinfo.c_start,'ok','linewidth',2);
catch
    display('No starting data')
end
errorbar(Tperiods,phv_av_rt,phv_avstd_rt*2,'-r','linewidth',2);
% errorbar(Tperiods,mean([vertcat(raytomo(:).phv)],2),std([vertcat(raytomo(:).phv)],0,2)*2,'-b','linewidth',2);
title('Isotropic Phase Velocity');
xlabel('Period (s)','fontsize',16);
ylabel('Phase Velocity (km/s)','fontsize',16);
set(gca,'fontsize',16,'linewidth',1.5);
xlim(flip(1./frange));
if comp{1}(1) == 'Z'
    ylim([3.4 4.3]);
elseif comp{1}(1) == 'T'
    ylim([3.8 4.7]);
end


%% Plot Azimuthal Data
fig6 = figure(6); clf;
% set(gcf,'position',[10         248        1203         457]);
set(gcf,'position',[10          11        1203         695]);
for iper = 1:length(periods)
    azi = raytomo(iper).azi;
    azi(azi<0) = azi(azi<0) + 360;
    dphv = (raytomo(iper).phv' - raytomo(iper).phv_iso) ./ raytomo(iper).phv_iso;
    dist = raytomo(iper).dist;
    subplot(3,4,iper); hold on;
    x = [0:360];
    % PHV FIT = a*(1+d*cosd(c*(x-e)))
    h2(2) = plot(x,A2_rt(iper)*cosd(2*(x-phi2_rt(iper)))*100,'-','color',[0.5 0.5 0.5],'linewidth',3);
    scatter(azi,dphv*100,30,dist,'filled'); hold on;
    
    if iper == 4
        ax = get(gca);
        pos = ax.Position;
        colorbar;
        set(gca,'Position',pos);
    end
    title([num2str(periods(iper),'%0.1f'),' s'],'fontsize',30);
    if iper == length(periods)
        xlabel('Azimuth (degrees)','fontsize',15);
    end
    if iper == 1
        ylabel('\delta{c}/c (%)','fontsize',15);
    end
    set(gca,'fontsize',18);
    xlim([0 360]);
    ylim([-5 5]);
    
end
