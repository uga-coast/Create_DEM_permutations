clear
clc
close all
datadir = 'C:\Users\vhh04855\Documents\Idiealized_cases\DEM_create_\data';
fname = fullfile(datadir,'barrier_island_data.mat');
load(fname)

%%
Barri_x = xq ;
zqt     = zq(5,:);
Barri_z = smooth(zqt,50)';

figure
hold on ; grid on ; box on
plot(Barri_x,Barri_z)
title('select base depth for the BI and height location')

%%
inptval = input('Select base depth for the BI');
b_basez = inptval; % base of the BI

x_base  = linspace(0,15000,15000);
b_base  = repmat(b_basez,1,numel(x_base));

inptval     = input('Select BI height/mid X location ');
xbarr_mid   = inptval; % base of the BI

[~,hidx]    = min(abs(Barri_x-xbarr_mid));
h_bi        =  Barri_z(hidx);

%% Intersection points between base and barrier island 

L1       = [Barri_x;Barri_z];
L2       = [x_base; b_base];
intb     = InterX(L1,L2);

dif  = intb(1,:)- xbarr_mid;
idf  =  dif<0; % intersection locations at front
idb  =  dif>0; % intersection locations at back

if( sum(idf) ==0) % no intersections at the front
    disp('No intersections. Please choose another base value ')
else
    t_int    = intb(1,idf);
    fint     = max(t_int); % front intersection location
    [~, lfx] = min(abs(Barri_x-fint));
end

    if( sum(idb) ==0) % no intersections at the back
        % need aux slope for back
       disp('No intersections. Please choose another base value ')
    else
        t_int = intb(1,idb);
        bint  = min(t_int); % back intersection location
        [~, lbx] = min(abs(Barri_x-bint));
    end

Barri_x1 = Barri_x(lfx:lbx); 
Barri_z1 = Barri_z(lfx:lbx);


%% interpolate the BI and change x coord 
x_resl   = 0.5; % resolution of the BI qurey points 
Barri_x0 = Barri_x1-Barri_x1(1);
Barri_z0 = Barri_z1;
nx       = ceil(Barri_x0(end)/x_resl); 
Barri_xq = linspace(0,Barri_x0(end),nx);
Barri_zq = interp1(Barri_x0,Barri_z0,Barri_xq,'spline');

%% find the width and height of the BI

% width
% width is defined as the width at the MSL
x_base   = linspace(0,15000,15000);
b_base   = zeros(1,numel(x_base)); % at MSL
L1       = [Barri_xq;Barri_zq];
L2       = [x_base; b_base];

[~,hidx ]  = min(abs( Barri_zq-h_bi)); % find the index of the new mid/height loc
h_bi       = Barri_zq(hidx); % barrier height after interpolation 
xbarr_mid  = Barri_xq(hidx); % barrier mid after interpolation and coord change
intw       = InterX(L1,L2);

dif      = intw(1,:)- xbarr_mid;
idf      =  dif<0; % intersection locations at front
idb      =  dif>0; % intersection locations at back

t_int    = intw(1,idf);
fint     = max(t_int); % front intersection location
[~, lfx] = min(abs(Barri_xq-fint));


t_int    = intw(1,idb);
bint     = min(t_int); % back intersection location
[~, lbx] = min(abs(Barri_xq-bint));

w_bi     = Barri_xq(lbx)- Barri_xq(lfx);
% height : select an appropriate index from the plot

%%
figure
hold on; grid on; box on
plot(Barri_xq,Barri_zq,'k')
plot(Barri_xq(hidx),Barri_zq(hidx),'or')
plot(Barri_xq(lfx:lbx),zeros(1,numel(Barri_xq(lfx:lbx))),'r')
%% 

Barri.zOrg      = Barri_zq - b_basez;
Barri.xOrg      = Barri_xq;
Barri.b_basez   = b_basez;
Barri.h_bi      = h_bi;
Barri.hidx      = hidx;
Barri.w_bi      = w_bi;


fname = fullfile(datadir,'Origin_BI_data_3.mat');
load(fname)
save(fname,'Barri')










