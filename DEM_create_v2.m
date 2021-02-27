clear all
clc
close all

wrkdir = 'C:\Users\vhh04855\Documents\Idiealized_cases\DEM_create_';
%% base profile
prof = load('C:\Users\vhh04855\OneDrive - University of Florida\data_and_example_code\Bar_v9data_NEW.mat','x','z');
%
len_y   = 1000;
xres    = 5;
yres    = 5;
prtrans = 19;
prtim   = 55;
%
zsel    = cell2mat(prof.z(prtim,prtrans));
xsel    = cell2mat(prof.x(prtim,prtrans));
len_x   = xsel(end);
x1      = flip(-(xsel-xsel(end)));
z1      = flip(zsel);
xq      = linspace(x1(1),x1(end),ceil(len_x/xres)+1);
zq      = interp1(x1,z1,xq,'spline');
yq      = linspace(0,len_y,ceil(len_y/yres)+1);


[ymesh,xmesh]   = meshgrid(yq,xq);
zmesh           = repmat(zq',1,length(yq));

r_norm          = 0.1;  % norminal roughness of base profile
r_base          = 0.25; %  masked roughness of base profile
zthresh         = 1.3;  % threshold height for roughness
rmesh           = repmat(r_norm,numel(xq),numel(yq));
% roughness mask for base profile
base_mask       = zmesh > zthresh;
rmesh(base_mask)=  r_base;

barrier_island  = 1; % Boolien include barrier islands
marsh           = 1;% Boolien include barrier islands
%% Barrier island
if (barrier_island ==1)
    f_slopeb     = 1/10; % Barrier island front slope
    b_slopeb     = 1/16; % Barrier island back slope
    l_slopeb     = 1/8;  % Barrier island lateral slope
    widthb       = 34;   % Barrier island width
    barrier_hmax = 1.8;  % Barrier height at summit
    barrier_pati = 3;    % number of barrier island partitions
    xbarrier_loc = 300;  % Location of the barrier island comparative to shore
    ch_precent   = 0.01; % channel width as a percentage of individual barrier island partition
    r_barr       = 0.3;  % roughness barrier island
    
    barrier_l    = len_y/barrier_pati;
    barrier_op   = barrier_l*ch_precent;
    [~,MSL_idx]  = min(abs(zq));
    xbarr_mid    = xq(MSL_idx) - xbarrier_loc;
    xbarr_frt    = xbarr_mid-widthb/2;
    xbarr_bac    = xbarr_mid+widthb/2;
    [~,bmid_idx] = min(abs(xq-xbarr_mid));
    [~,bfrt_idx] = min(abs(xq-xbarr_frt));
    [~,bbac_idx] = min(abs(xq-xbarr_bac));
    
    channel_h2   = barrier_hmax:-yres*l_slopeb:zq(bmid_idx);
    chngrid      = floor(barrier_op/yres)+1;
    if(chngrid) <3
        chngrid = 3;
    end
    channel_h3   = repmat(zq(bmid_idx),1,chngrid);
    channel_h4   = flip(barrier_hmax:-yres*l_slopeb:zq(bmid_idx));
    
    channel_h1   = [channel_h2 channel_h3 channel_h4];
    channel_len  = yres*(numel(channel_h1)-1);
    nchannel     = barrier_pati-1;
    
    barrier_h1   = repmat(barrier_hmax,1,numel(yq));
    nid          = floor((numel(yq)-1)/barrier_pati);
    
    if (barrier_pati==1)
        % one long barrier island
    else
        for k1 = 1:barrier_pati-1
            
            cloc = nid*k1;
            i3 = (numel(channel_h1)-1)/2;
            idx  =  [(cloc-i3):cloc  (cloc+1):(cloc+i3)];
            
            barrier_h1(idx) = channel_h1;
            
        end
    end
    
    for j1 = 1:numel(yq)
        barrier_h    = barrier_h1(j1);
        xtempf       = xq(1:bfrt_idx);
        ztempf       = zq(1:bfrt_idx);
        cint         = barrier_h - xq(bfrt_idx)*f_slopeb;
        zbarr_frt    = xtempf*f_slopeb+cint;
        
        idt          = ztempf - zbarr_frt<0;
        xqbf         = xtempf(idt);
        zqbf         = zbarr_frt(idt);
        
        xtempb       = xq(bbac_idx:length(xq));
        ztempb       = zq(bbac_idx:length(xq));
        cint         = barrier_h - xq(bbac_idx)*-b_slopeb;
        zbarr_bac    = xtempb*-b_slopeb+cint;
        
        idt          = ztempb - zbarr_bac<0;
        xqbb         = xtempb(idt);
        zqbb         = zbarr_bac(idt);
        
        xqbm         = xq(bfrt_idx+1:bbac_idx-1);
        zqbm         = repmat(barrier_h,1,numel(xqbm));
        
        zbarrier     = [zqbf zqbm zqbb] ;
        xbarrier     = [xqbf xqbm xqbb] ;
        
        barr_mask    = zbarrier >0;
        rbarrier(barr_mask) =  r_barr;
        rbarrier(~barr_mask)=  r_norm;
        
        % superposition-barrier island- morphology
        for i1 =1:numel(xq)
            idx = xq(i1)== xbarrier;
            if any(idx)
                znew1(i1) = zbarrier(idx);
            else
                znew1(i1) = zq(i1);
            end
        end
        zmeshnew1(:,j1) = znew1;
        
        % superposition-barrier island- roughness
        for i1 = 1:numel(xq)
            idx = xq(i1)== xbarrier;
            if any(idx)
                rmeshnew1(i1,j1) = rbarrier(idx);
            else
                rmeshnew1(i1,j1) = rmesh(i1,j1);
            end
            
        end
    end
else
    zmeshnew1 = zmesh;
    rmeshnew1 = rmesh;
end

%% marsh
if (marsh ==1)
    f_slopem     = 1/10; % marsh front slope
    b_slopem     = 1/16; % marsh back slope
    width        = 67;   % marsh width
    marsh_h      = -0.1; % marsh height
    marsh_loc    = 100;  % Location of the marsh comparative to shore
    r_marsh      = 0.4;  % marsh roughness
    
    [~,MSL_idx]  = min(abs(zq));
    xmars_mid    = xq(MSL_idx) - marsh_loc;
    xmars_frt    = xmars_mid-width/2;
    xmars_bac    = xmars_mid+width/2;
    
    
    [~,mfrt_idx] = min(abs(xq-xmars_frt));
    [~,mbac_idx] = min(abs(xq-xmars_bac));
    
    xtempf       = xq(1:mfrt_idx);
    ztempf       = zq(1:mfrt_idx);
    cint         = marsh_h - xq(mfrt_idx)*f_slopem;
    zmars_frt    = xtempf*f_slopem+cint;
    
    idt          = ztempf - zmars_frt<0;
    xqmf         = xtempf(idt);
    zqmf         = zmars_frt(idt);
    
    xtempb       = xq(mbac_idx:length(xq));
    ztempb       = zq(mbac_idx:length(xq));
    cint         = marsh_h - xq(mbac_idx)*-b_slopem;
    zmars_bac    = xtempb*-b_slopem+cint;
    
    idt          = ztempb - zmars_bac<0;
    xqmb         = xtempb(idt);
    zqmb         = zmars_bac(idt);
    
    xqmm         = xq(mfrt_idx+1:mbac_idx-1);
    zqmm         = repmat(marsh_h,1,numel(xqmm));
    
    zmarsh       = [zqmf zqmm zqmb] ;
    xmarsh       = [xqmf xqmm xqmb] ;
    
    mars_mask    = zmarsh == zqmm(1);
    rmarsh(mars_mask) =  r_marsh;
    rmarsh(~mars_mask)=  r_norm;
    
    %superposition-marsh- morphology
    for j1 = 1:numel(yq)
        ztemp = zmeshnew1(:,j1);
        for i1 =1:numel(xq)
            idx = xq(i1)== xmarsh;
            if any(idx)
                znew2(i1) = zmarsh(idx);
            else
                znew2(i1) = ztemp(i1);
            end
        end
        zmeshnew2(:,j1) = znew2;
    end
    
    %superposition-marsh- roughness
    for j1 = 1:numel(yq)
        for i1 =1:numel(xq)
            idx = xq(i1)== xmarsh;
            if any(idx)
                rmeshnew2(i1,j1) = rmarsh(idx);
            else
                rmeshnew2(i1,j1) = rmeshnew1(i1,j1);
            end
        end
    end
    
else
    zmeshnew2 = zmeshnew1;
    rmeshnew2 = rmeshnew1;
end


%%
hFig=figure;
Y = 29.7; X = 21.0+15; %# A4 paper size
xSize = X ; ySize = Y ; %# figure size on paper (widht & height)
set(hFig, 'PaperUnits','centimeters');
set(hFig, 'PaperSize',[X Y/2]);
set(hFig, 'PaperPosition',[0.25 0.25 xSize ySize/2]);
set(hFig, 'PaperOrientation','portrait');

han(4) = subplot(2,3,[1 2 3])
hold on; grid on; box on
surf(ymesh,xmesh,zmeshnew2)
colormap winter
colorbar
view(-45,50)
xlabel('y [m]','FontSize',12)
ylabel('x [m]','FontSize',12)
zlabel('z [m]','FontSize',12)

han(3) = subplot(2,3,4)
hold on; grid on; box on
plot(xq,zq,'LineWidth',1,'Color','k','LineStyle','-')
xlabel('x [m]','FontSize',12)
ylabel('y [m]','FontSize',12)
text( 10 ,4,'base profile(cs)','FontSize',12)

han(2) = subplot(2,3,5)
hold on; grid on; box on
plot(xq,znew2,'LineWidth',1,'Color','k','LineStyle','-')
plot(xbarrier,zbarrier,'LineWidth',1,'Color','r','LineStyle','-')
plot(xmarsh,zmarsh,'LineWidth',1,'Color','g','LineStyle','-')
xlabel('x [m]','FontSize',12)
ylabel('y [m]','FontSize',12)
text( 10 ,4,'super imposed profile(cs)','FontSize',12)

han(1) = subplot(2,3,6)
contourf(ymesh,xmesh,rmeshnew2)
colormap(parula)
colorbar
xlabel('x [m]','FontSize',12)
ylabel('y [m]','FontSize',12)
text( 100 ,100,'roughness','FontSize',12,'Color','w')

for q1 = 1:4
    pos(q1,:) = get(han(q1),'Position');
    
end

for q1 = 1:3
    pos(q1,2)   = pos(q1,2) *.7;
    pos(q1,4)   = pos(q1,4) *.7;
end

pos(4,1)   = pos(4,1)*.7;
pos(4,2)   = pos(3,2) +pos(3,4)+0.06;
pos(4,3)   = pos(4,3)*1.3;
pos(4,4)   = pos(4,4)*1.7;

for q1 = 1:4
    set(han(q1),'Position',pos(q1,:));
end


