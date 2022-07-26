clear
clc
close all
%% Main Controls and Parameters

ptype           = 2; % Sets profile type. (0 = Duck; 1 = Tyndall; 2 = Idealized)
barrier_island  = 1; % Boolean include barrier islands
btype           = 1; % Barrier Island profile type (1 = Real World Transect, 2 = Idealized Transect)
marsh           = 0; % Boolean include marsh
back_bay        = 1; % Boolean include back bay (Can only be 1 if barrier_island == 1 && ptype == 2)
mainland        = 1; % Boolean include mainland (Can only be 1 if barrier_island ==1 && ptype == 2)
write2file      = 0; % Specifies whether to write landscape to .xyz file
SMS             = 1; % Specifies whether output file is intended for use with SMS to interpolate to mesh
filename        = "Test"; %Output filename. .xyz extension will be added automatically
wrkdir          = 'C:\Users\sdm08233\Documents\MATLAB\DEM_create_permutations-main\Exp2'; %Working Directory.  Where files will be output.
makeplot        = 1; % Specifies whether to draw surface plot of resulting landscape
saveplot        = 0; % Specifies whether to save the surfaceplot from makeplot.  Can only be 1 if makeplot is 1.

%% Set up Base Profile

% Required base profile arameters
len_y           = 160000;  % Domain length in y direction (longshore)
xres            = 25;      % Target x grid resolution (Will change if x_len is not evenly divisible by xres)
yres            = 100;     % Target y grid resolution (Will change if y_len is not evenly divisible by yres)
r_norm          = 0.1;     % normal roughness of base profile
r_base          = 0.25;    % masked roughness of base profile
zthresh         = 1.3;     % threshold height for roughness

% Additional parameters for Idealized base profiles

len_x           = 80000;      % Domain/base profile length in x direction (cross shore)
DOC             = -2;         % Depth of Closure 
DOCSlope        = 5/800;      % Slope from land to Depth of Closure
OSSlope         = 195/137600; % Slope from DOC offshore
bay_Depth       = -6;         % Depth of back bay
len_Bay         = 5500;%, 500, 1000, 2000, 3000, 5000, 10000, 15000, 20000];      % Cross Shore length of back bay
startDepth      = 0;          % Depth to start profile
w_req           = 3000;       % BI width required
hmax_ml         = 8;          % Max mainland height - Along mainland boundary
mlslope         = .001;       % Mainland slope

datadir          = 'C:\Users\sdm08233\Documents\MATLAB\DEM_create_permutations-main\data';
for j = 1%:9
    %Checking that bay depth can be reached given chosen bay length and
    %slopes
    bay_override = 0;
    halfBay = abs(len_Bay/2);
    DOC_len = abs(DOC/DOCSlope);
    deepbay_len = [];
    OSoverride = 0;
    if bay_Depth < DOC
        deepbay_len = abs(bay_Depth - DOC)/OSSlope;
        if (deepbay_len + DOC_len) > halfBay
            warning("The chosen Bay Length is too short to reach desired bay depth at the chosen slopes. Consider revising. Landscape generation will continue by lengthening bay to accomadate chosen bay depth.")
            bay_override = 1;
            if DOC_len > halfBay
                OSoverride = 1;
            end
        end
    elseif DOC_len > halfBay
        warning("The chosen Bay Length is too short to reach desired bay depth at the chosen slopes. Consider revising. Landscape generation will continue by lengthening bay to accomadate chosen bay depth.")
        bay_override = 1;
    end
    
    %Calling the scipts to set up the base profile
    if ptype == 2
        [base_data, xbarr_mid, bay_check]     = Func_setup_base_profile(ptype,barrier_island,mainland,hmax_ml,mlslope,back_bay,len_x,DOC,DOCSlope,...
            OSSlope,bay_Depth,len_Bay(j),startDepth,w_req,len_y,xres,yres,r_norm,r_base,zthresh,datadir, bay_override, DOC_len, halfBay, deepbay_len, OSoverride);
    else
        [base_data, xbarr_mid, ~]     = Func_setup_base_profile(ptype,barrier_island,mainland,hmax_ml,mlslope,back_bay,len_x,DOC,DOCSlope,...
            OSSlope,bay_Depth,len_Bay(j),startDepth,w_req,len_y,xres,yres,r_norm,r_base,zthresh,datadir, bay_override, DOC_len, halfBay, deepbay_len, OSoverride);
    end
    base_data.yres  = yres; 
    base_data.len_y = len_y;
%% Set up Barrier Island    
    if((barrier_island == 1))
        % Required barrier island parameters
        h_req         =  3.29*1.03;      % BI height required
        xbarrier_loc  = len_Bay;  % Location of the barrier island comparative to shore
        aux_slopeb    = .003;     % BI aux slope
        chn_slope     = 1/10;    % BI channel slope
        barrier_pati  = 4;    % number of barrier island partitions
        inlet_depth   = -5.35;    % Depth of inlet
        chn_bwidth    = 2400;     % base width of the channel
        r_barr        = 0.3;      % roughness barrier island
        dist2edge     = 16000;    %Distance from end of BI to edge of mesh
        BI_data       = 'Origin_BI_data_4';
        
        % Calling Barrier Island Generation Script
        for i = 1
            close all
            filename = strcat('Bay',num2str(len_Bay(j)),'Dune',num2str(h_req(i),1)); %Output filename. .xyz extension will be added automatically
            [xbarrier, zbarrier, zmeshnew1, rmeshnew1] = Func_add_barrier_island (ptype, ...
                bay_Depth,back_bay,dist2edge,bay_check,bay_Depth,inlet_depth,xres,xbarr_mid, ...
                w_req,h_req(i),xbarrier_loc ,barrier_pati,aux_slopeb,chn_bwidth, ...
                chn_slope,r_barr,r_norm,zthresh,base_data,datadir,BI_data);

            zmeshnew2 = zmeshnew1;
            rmeshnew2 = rmeshnew1;
            %% Create Marsh
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
            %% Plotting
            if makeplot == 1
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
                set(hFig, 'visible','on');
                set(hFig, 'color', 'None');

                %han(4) = subplot(2,3,[1 2 3]);
                hold on; grid on; box on
                surf(ymesh,xmesh,zmeshnew2+0.003,'FaceAlpha',0.9, 'EdgeColor', 'none')
                % [cmap,climits] = demcmap(zmeshnew2);
                % shading interp
                c1 = colorbar;
                colormap (jet)
                caxis([-20 10])
                ylabel(c1, 'depth [m]','FontSize',12)
                view(-45,50)
                ylim([50000 80000])
                zlim([-40 5])

                xlabel('y [m]','FontSize',12)
                ylabel('x [m]','FontSize',12)
                zlabel('z [m]','FontSize',12)
                title("Generated Landscape")
            end
            %% Writing to file
            if write2file == 1 && SMS == 1
                zmeshnew2 = zmeshnew2.*(-1);
                zmeshnew2 = flipud(zmeshnew2);
                %Resize to fit global coordinates
                xmesh = xmesh .* (1.4876056/len_y);
                ymesh = ymesh .* (1.4876056/len_y);
                %Translate to fit global coordinates
                xmesh = xmesh - 81.43998312;
                ymesh = ymesh + 30.71494931;
                X = xmesh(:);
                Y = ymesh(:);
                %Rotate to fit clobal coordinates
                % create a matrix of these points, which will be useful in future calculations
                v = [X.';Y.'];
                % choose a point which will be the center of rotation
                x_center = -81.43998312;
                y_center = 30.71494931;
                % create a matrix which will be used later in calculations
                center = repmat([x_center; y_center], 1, length(X));
                % define a 60 degree counter-clockwise rotation matrix
                theta = -26;       % pi/3 radians = 60 degrees
                R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
                % do the rotation...
                s = v - center;     % shift points in the plane so that the center of rotation is at the origin
                so = R*s;           % apply the rotation about the origin
                vo = so + center;   % shift again so the origin goes back to the desired center of rotation
                % this can be done in one line as:
                % vo = R*(v - center) + center
                % pick out the vectors of rotated x- and y-data
                x_rotated = vo(1,:);
                X = x_rotated.';
                y_rotated = vo(2,:);
                Y = y_rotated.';
                Z = zmeshnew2(:);
                P = [X, Y, Z];
                fid = fopen(fullfile(wrkdir, strcat(filename, '.xyz')), 'w');
                fprintf(fid, [repmat('%.5f\t', 1, size(P,2)-1) '%.5f\n'], P.');
                fclose(fid);
            elseif write2file == 1
                zmeshnew2 = zmeshnew2.*(-1);
                X = xmesh(:);
                Y = ymesh(:);
                Z = zmeshnew2(:);
                P = [X, Y, Z];
                fid = fopen(fullfile(wrkdir, strcat(filename, '.xyz')), 'w');
                fprintf(fid, [repmat('%.5f\t', 1, size(P,2)-1) '%.5f\n'], P.');
                fclose(fid);
            end
        end
    else
        zmeshnew1 = base_data.zmesh;
        rmeshnew1 = base_data.rmesh;
    end
end
%% Optional Additional Plotting 
%{
if makeplot == 1
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
    set(hFig, 'visible','on');
    set(hFig, 'color', 'None');

    %han(4) = subplot(2,3,[1 2 3]);
    hold on; grid on; box on
    surf(ymesh,xmesh,zmeshnew2+0.003,'FaceAlpha',0.9, 'EdgeColor', 'none')
    % [cmap,climits] = demcmap(zmeshnew2);
    % shading interp
    c1 = colorbar;
    colormap (jet)
    caxis([-20 10])
    ylabel(c1, 'depth [m]','FontSize',12)
    view(-45,50)
    ylim([60000 80000])
    zlim([-40 5])

    xlabel('y [m]','FontSize',12)
    ylabel('x [m]','FontSize',12)
    zlabel('z [m]','FontSize',12)
    title("Generated Landscape")

    %han(3) = subplot(2,3,4);
    %hold on; grid on; box on
    %plot(xq,zq,'LineWidth',1,'Color','k','LineStyle','-')
    %xlabel('x [m]','FontSize',12)
    %ylabel('y [m]','FontSize',12)
    %text( 10 ,9,'base profile(cs)','FontSize',12)
    %ylim([-100 10])
    %xlim([0 max(xq)])

    %han(2) = subplot(2,3,5);
    %hold on; grid on; box on
    %plot(xq,zmeshnew2(:,1),'LineWidth',1,'Color','k','LineStyle','-')
    %plot(xbarrier,zbarrier,'LineWidth',1,'Color','r','LineStyle','-')
    %plot(xmarsh,zmarsh,'LineWidth',1,'Color','g','LineStyle','-')
    %xlabel('x [m]','FontSize',12)
    %ylabel('y [m]','FontSize',12)
    %text( 10 ,9,'super imposed profile(cs)','FontSize',12)
    %ylim([-20 10])
    %xlim([0 max(xq)])

    %han(1) = subplot(2,3,6);
    %contourf(ymesh,xmesh,rmeshnew2)
    %c2 = colorbar ;
    %colormap(han(1),parula)
    %ylabel(c2, 'Normal roughness','FontSize',12)
    %xlabel('x [m]','FontSize',12)
    %ylabel('y [m]','FontSize',12)
    %text( 100 ,600,'roughness','FontSize',12,'Color','w')

    %pos = NaN(4,4);

    %for q1 = 1:4
     %   pos(q1,:) = get(han(q1),'Position');
    %end

    %for q1 = 1:3
    %    pos(q1,2)   = pos(q1,2) *.6;
    %    pos(q1,4)   = pos(q1,4) *.7;
    %end

    %pos(4,1)   = pos(4,1)*.7;
    %pos(4,2)   = pos(3,2) +pos(3,4)+0.1;
    %pos(4,3)   = pos(4,3)*1.3;
    %pos(4,4)   = pos(4,4)*1.7;


    %for q1 = 1:4
    %    set(han(q1),'Position',pos(q1,:));
    %end
end

if saveplot == 1
    if makeplot == 0
        error("No plot to save, please change makeplot to 1")
    else
        ff1 = ('Fig_dem');
        fname=fullfile(wrkdir, ff1);
        print('-dpng','-r300', fname);
    end
end
%}
%% Writing Output to File - Optional
%{
if write2file == 1 && SMS == 1
    zmeshnew2 = zmeshnew2.*(-1);
    zmeshnew2 = flipud(zmeshnew2);
    xmesh = xmesh .* (1.4876056/len_y);
    ymesh = ymesh .* (1.4876056/len_y);
    xmesh = xmesh - 81.43998312;
    ymesh = ymesh + 30.71494931;
    X = xmesh(:);
    Y = ymesh(:);
    v = [X.';Y.'];
    % choose a point which will be the center of rotation
    x_center = -81.43998312;
    y_center = 30.71494931;
    % create a matrix which will be used later in calculations
    center = repmat([x_center; y_center], 1, length(X));
    % define a 60 degree counter-clockwise rotation matrix
    theta = -26;       % pi/3 radians = 60 degrees
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    % do the rotation...
    s = v - center;     % shift points in the plane so that the center of rotation is at the origin
    so = R*s;           % apply the rotation about the origin
    vo = so + center;   % shift again so the origin goes back to the desired center of rotation
    % this can be done in one line as:
    % vo = R*(v - center) + center
    % pick out the vectors of rotated x- and y-data
    x_rotated = vo(1,:);
    X = x_rotated.';
    y_rotated = vo(2,:);
    Y = y_rotated.';
    Z = zmeshnew2(:);
    P = [X, Y, Z];
    fid = fopen(fullfile(wrkdir, filename + '.xyz'), 'w');
    fprintf(fid, [repmat('%.5f\t', 1, size(P,2)-1) '%.5f\n'], P.');
    fclose(fid);
elseif write2file == 1
    zmeshnew2 = zmeshnew2.*(-1);
    X = xmesh(:);
    Y = ymesh(:);
    Z = zmeshnew2(:);
    P = [X, Y, Z];
    fid = fopen(fullfile(wrkdir, filename + '.xyz'), 'w');
    fprintf(fid, [repmat('%.5f\t', 1, size(P,2)-1) '%.5f\n'], P.');
    fclose(fid);
end
%}
