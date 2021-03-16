clear
clc
close all
%% set up base profile

len_y           = 2000; % length in y direction
xres            = 20;    % x grid resolution
yres            = 20;    % y grid resolution
r_norm          = 0.1;  % norminal roughness of base profile
r_base          = 0.25; % masked roughness of base profile
zthresh         = 1.3;  % threshold height for roughness
datadir          = 'C:\Users\vhh04855\Documents\Idiealized_cases\DEM_create_\data';
[base_data]     = Func_setup_base_profile(len_y,xres,yres,r_norm,r_base,zthresh,datadir);
base_data.yres  = yres; 
base_data.len_y = len_y;
%% main controls and parameter
barrier_island  = 1; % Boolien include barrier islands
marsh           = 1; % Boolien include barrier islands

%% Create barrier island
if((barrier_island ==1))
    w_req         = 600;   % BI width required
    h_req         = 2.6;  % BI height required
    xbarrier_loc  = 2000;  % Location of the barrier island comparative to shore
    aux_slopeb    = 1/10; % BI aux slope
    chn_slope     = 1/8;  % BI channel slope
    barrier_pati  = 3;    % number of barrier island partitions
    chn_bwidth    = 60;   % base width of the channel
    r_barr        = 0.3;  % roughness barrier island
    BI_data       = 'Origin_BI_data_3';
    [xbarrier, zbarrier, zmeshnew1, rmeshnew1] = Func_add_barrier_island (w_req, h_req, xbarrier_loc , barrier_pati, aux_slopeb , chn_bwidth, chn_slope, r_barr, r_norm, zthresh, base_data, datadir, BI_data);
    
else
    zmeshnew1 = base_data.zmesh;
    rmeshnew1 = base_data.rmesh;
end
%% create Marsh

if (marsh ==1)
    f_slopem     = 1/50; % marsh front slope
    b_slopem     = 1/660; % marsh back slope
    width        = 700;   % marsh width
    marsh_h      = -0.1; % marsh height
    marsh_loc    = 70;  % Location of the marsh comparative to shore
    r_marsh      = 0.4;  % marsh roughness
    [zmeshnew2, rmeshnew2, xmarsh, zmarsh] = Func_add_marsh (f_slopem, b_slopem, width, marsh_h, marsh_loc ,r_marsh, r_norm, zmeshnew1, rmeshnew1, base_data);
else
    zmeshnew2 = zmeshnew1;
    rmeshnew2 = rmeshnew1;
end

%% Plot 

xmesh   = base_data.xmesh;
ymesh   = base_data.ymesh;
zmesh   = base_data.zmesh;
rmesh   = base_data.rmesh;
xq      = base_data.xq ;
yq      = base_data.yq ;
zq      = base_data.zq;


hFig=figure;
Y = 29.7; X = 21.0+15; %# A4 paper size
xSize = X ; ySize = Y ; %# figure size on paper (widht & height)
set(hFig, 'PaperUnits','centimeters');
set(hFig, 'PaperSize',[X Y/2]);
set(hFig, 'PaperPosition',[0.1 0.1 xSize ySize/2]);
set(hFig, 'PaperOrientation','portrait');
set(hFig, 'visible','off');

han(4) = subplot(2,3,[1 2 3]);
hold on; grid on; box on
surf(ymesh,xmesh,zmeshnew2+0.003,'FaceAlpha',0.9)
% [cmap,climits] = demcmap(zmeshnew2);
% shading interp
c1 = colorbar;
colormap (jet)
ylabel(c1, 'depth [m]','FontSize',12)
view(-45,50)

zlim([-20 10])

xlabel('y [m]','FontSize',12)
ylabel('x [m]','FontSize',12)
zlabel('z [m]','FontSize',12)

han(3) = subplot(2,3,4);
hold on; grid on; box on
plot(xq,zq,'LineWidth',1,'Color','k','LineStyle','-')
xlabel('x [m]','FontSize',12)
ylabel('y [m]','FontSize',12)
text( 10 ,9,'base profile(cs)','FontSize',12)
ylim([-20 10])
xlim([0 max(xq)])

han(2) = subplot(2,3,5);
hold on; grid on; box on
plot(xq,zmeshnew2(:,1),'LineWidth',1,'Color','k','LineStyle','-')
plot(xbarrier,zbarrier,'LineWidth',1,'Color','r','LineStyle','-')
plot(xmarsh,zmarsh,'LineWidth',1,'Color','g','LineStyle','-')
xlabel('x [m]','FontSize',12)
ylabel('y [m]','FontSize',12)
text( 10 ,9,'super imposed profile(cs)','FontSize',12)
ylim([-20 10])
xlim([0 max(xq)])

han(1) = subplot(2,3,6);
contourf(ymesh,xmesh,rmeshnew2)
c2 = colorbar ;
colormap(han(1),parula)
ylabel(c2, 'Norminal roughness','FontSize',12)
xlabel('x [m]','FontSize',12)
ylabel('y [m]','FontSize',12)
text( 100 ,600,'roughness','FontSize',12,'Color','w')

for q1 = 1:4
    pos(q1,:) = get(han(q1),'Position');
    
end

for q1 = 1:3
    pos(q1,2)   = pos(q1,2) *.6;
    pos(q1,4)   = pos(q1,4) *.7;
end

pos(4,1)   = pos(4,1)*.7;
pos(4,2)   = pos(3,2) +pos(3,4)+0.1;
pos(4,3)   = pos(4,3)*1.3;
pos(4,4)   = pos(4,4)*1.7;

for q1 = 1:4
    set(han(q1),'Position',pos(q1,:));
end


wrkdir = 'C:\Users\vhh04855\OneDrive - University of Georgia\N-EWN\Notes_and_stuff\Meeting_stuff\DEM_arcmap';
ff1 = ('Fig_dem');
fname=fullfile(wrkdir, ff1);
print('-dpng','-r300', fname);

