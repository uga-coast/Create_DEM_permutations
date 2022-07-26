function [xbarrier, zbarrier, zmeshnew1, rmeshnew1] = Func_add_barrier_island (ptype, ...
    bay_depth,back_bay,dist2edge,bay_check,bay_Depth,inlet_depth,xres,xbarr_mid,w_req, ...
    h_req,xbarrier_loc,barrier_pati,aux_slopeb,chn_bwidth,chn_slope, ...
    r_barr,r_norm,zthresh,base_data,datadir,BI_data)

zmesh   = base_data.zmesh;
rmesh   = base_data.rmesh;
xq      = base_data.xq ;
yq      = base_data.yq ;
zq      = base_data.zq;
yres    = base_data.yres;
len_y   = base_data.len_y;

% import BI profile data
load(fullfile(datadir,BI_data));

x1_BI         = Barri.xOrg; %BIdat.Barri_xOrg;
z1_BI         = Barri.zOrg; %BIdat.Barri_zOrg;
w_bi          = Barri.w_bi; % original width of the BI width at MSL
h_bi          = max(Barri.zOrg) + Barri.b_basez;%Barri.h_bi; %BIdat.h_bi; % original hiehg of the BI
h_bidx        = (h_bi-Barri.b_basez) == Barri.zOrg; %max BI height index


chn_bwidth    = floor(chn_bwidth/yres)*yres;

aspr_h        = h_req/h_bi; % ratio to change the original to suit the required
aspr_w        = w_req/w_bi; % ratio to change the original to suit the required
x2_BI         = x1_BI*aspr_w;
z2_BI         = z1_BI*aspr_h;

figure
hold on; grid on
%plot(xq,zq)
%plot(xbarr_mid,0,'o')
plot(x2_BI,z2_BI)

if ptype == 2 && back_bay == 1
    [~,MSL_idx]  = min(abs(zq));
    %xbarr_mid    = xq(MSL_idx) + (w_req/2); % X val of the mid point location of the BI
    xshft        = abs(x2_BI(Barri.hidx)- xbarr_mid);
    % shifted BI coordinates to put the BI at the correct location & height
    x3_BI        = x2_BI + xshft;
    z3_BI        = z2_BI + Barri.b_basez*aspr_h;
else
    [~,MSL_idx]  = min(abs(zq));
    %xbarr_mid    = xq(MSL_idx) - xbarrier_loc; % X val of the mid point location of the BI
    xshft        = abs(x2_BI(Barri.hidx) - xbarr_mid);
    % shifted BI coordinates to put the BI at the correct location & height
    x3_BI        = x2_BI + xshft;
    z3_BI        = z2_BI + Barri.b_basez*aspr_h;
end

%figure
%hold on; grid on
plot(xq,zq)
plot(xbarr_mid,0,'o')
plot(x3_BI,z3_BI)

if x3_BI(1) < xq(1)
    error("Barrier Island resides further offshore than extent of transect."...
        + " Please consider providing a longer transect, moving island closer to shoreline, or creating an idealized transect.")
end

% check if the profiles intersect
L1       = [xq ; zq];
L2       = [x3_BI; z3_BI];
int1     = Func_InterX(L1,L2);

if isempty(int1)
    % no intersection- needs auxilary slopes to extend until base prof for
    % front and back
    xt1 = linspace(0,x3_BI(1),ceil(x3_BI(1)/.5));
    zt1 = z3_BI(1) + aux_slopeb*(xt1-xt1(end));
    xt2 = linspace(x3_BI(end),xq(end),ceil((xq(end)-x3_BI(end))/.5));
    zt2 = z3_BI(end) - aux_slopeb*(xt2-xt2(1));
      
    L1       = [xq; zq];
    L2       = [xt1; zt1];
    int2     = Func_InterX(L1,L2);
    fint     = int2(1,end);
    if zt1 > zq(1)
        error("BI slope is unable to intersect with transect."...
            + " Consider lengthening transect, increasing BI slope, moving BI closer to shore, or using idealized transect.")
    end
    [~, lfx] = min(abs(xt1-fint));
    
    L2       = [xt2; zt2];
    int2     = Func_InterX(L1,L2);
    bint     = int2(1,1);
    [~, lbx] = min(abs(xt2-bint));
    
    xt_bi = [ xt1(lfx:end-1) x3_BI xt2(2:lbx) ];
    zt_bi = [ zt1(lfx:end-1) z3_BI zt2(2:lbx) ];
    
else
    % at least one profile intersection exsist
    dif  = int1(1,:)- xbarr_mid;
    idf  =  dif<0; % intersection locations at front
    idb  =  dif>0; % intersection locations at back
    
    if( sum(idf) ==0) % no intersections at the front
        % need aux slope for front
        xt1 = abs(linspace(0,x3_BI(1),ceil(abs(x3_BI(1)/0.5))));
        zt1 = z3_BI(1) +aux_slopeb*(xt1-xt1(end));
        
        L1       = [xq ; zq];
        L2       = [xt1; zt1];
        int3     = Func_InterX(L1,L2);
        fint     = int3(1,end);
        if isempty(int3)
            error("Front end of BI does not intersect with transect.")
        end
        [~, lfx] = min(abs(xt1-fint));
        
        xtf_bi = [xt1(lfx:end-1) x3_BI(1:Barri.hidx)];
        ztf_bi = [zt1(lfx:end-1) z3_BI(1:Barri.hidx)];
    else
        t_int    = int1(1,idf);
        fint     = max(t_int); % front intersection location
        [~, lfx] = min(abs(x3_BI-fint));
        xtf_bi   = x3_BI(lfx:Barri.hidx);
        ztf_bi   = z3_BI(lfx:Barri.hidx);
    end
    
    if( sum(idb) ==0) % no intersections at the back
        % need aux slope for back
        xt2 = linspace(x3_BI(end),xq(end),ceil((xq(end)-x3_BI(end))/0.5));
        zt2 = z3_BI(end) -aux_slopeb*(xt2-xt2(1));
        
        L1       = [xq ; zq];
        L2       = [xt2; zt2];
        int4     = Func_InterX(L1,L2);
        bint     = int4(1,1);
        [~, lbx] = min(abs(xt2-bint));
        
        xtb_bi   = [x3_BI(Barri.hidx+1:end) xt2(2:lbx)];
        ztb_bi   = [z3_BI(Barri.hidx+1:end) zt2(2:lbx)];
    else
        t_int = int1(1,idb);
        bint  = min(t_int); % back intersection location
        [~, lbx] = min(abs(x3_BI-bint));
        xtb_bi   = x3_BI(Barri.hidx+1:lbx);
        ztb_bi   = z3_BI(Barri.hidx+1:lbx);
    end
    
    xt_bi = [xtf_bi xtb_bi];
    zt_bi = [ztf_bi ztb_bi];
    
end
%%
x_BI = xt_bi;
z_BI = zt_bi;
figure
hold on; grid on
ax = gca;
ax.FontSize = 24;
xlabel("Cross-shore distance (m)")
ylabel("Elevation (m)")
plot(x_BI,z_BI,'LineWidth',3)
plot(xq,zq,'LineWidth',3)
legend("Superimposed Barrier Island", "Base Profile")
ylim([-10,10])
xlim([55000,80000])
set(gcf, "Color", "None")

%%

% interploate the BI and base profile to same grid
zq_BI   = interp1(x_BI,z_BI,xq,'spline',-99);
idz     = zq_BI == -99;
tidx    = 1:numel(xq);
BIidx   = tidx(~idz);
z_comb  = zq;
z_comb(~idz) = zq_BI(~idz);
z_comb(idz)  = zq(idz);

% create one long barrier island
zmesh_BI  = repmat(z_comb',1,length(yq));

% channel creation
OSdepth = 1;
OSdepth2 = 1;
ISdepth = 1;
ISdepth2 = 1;
count = 0;
inlet_check = (inlet_depth > bay_Depth); % Returns 1 if inlet depth is shallower than bay depth, 0 if inlet depth is deeper than bay depth
while OSdepth > inlet_depth
    OSchidx = BIidx(1) - count;
    count = count + 1;
    OSdepth = z_comb(BIidx(1) - count);
end
count = 0;
if inlet_check == 1 % If inlet depth is shallower than bay depth ISdepth(index) will continue counting until it gets to the inlet depth
    while ISdepth > inlet_depth
        ISchidx = BIidx(end) + count;
        count = count + 1;
        ISdepth = z_comb(BIidx(end) + count);
    end
else
    try
        count2 = count;
        while ISdepth > bay_Depth % If inlet depth is deeper than or equal to bay depth ISdepth(index) will continue counting until it gets to the bay depth
            ISchidx = BIidx(end) + count;
            count = count + 1;
            ISdepth = z_comb(BIidx(end) + count);
        end
    catch ME
        count = count2;
        diff = 1;
        while diff > 0
            ISchidx = BIidx(end) + count;
            count = count + 1;
            ISdepth = z_comb(BIidx(end) + count);
            diff = z_comb(ISchidx) - ISdepth;
        end
    end
                
end

chidx = [OSchidx:ISchidx];
zt_BI = z_comb(chidx);
zt_bs = repelem(inlet_depth, length(chidx));
delzt = zt_BI-zt_bs;
chn_bwidth1 = yres:yres:chn_bwidth-yres ;
dist2edge1 = yres:yres:dist2edge;
chidx2 = chidx;
zq_BI2 = zq(chidx);
zt_BI2 = zt_BI;
zt_bs2 = zt_bs;
delzt2 = delzt;

%This creates an index to create a cut from the BI height down to the bay
%depth along the edges of the grid.  This only has an effect when back_bay
%==1
if back_bay == 1
    count = 0;
    while OSdepth2 > bay_depth
        OSchidx2 = BIidx(1) - count;
        count = count + 1;
        OSdepth2 = z_comb(BIidx(1) - count);
    end
    count = 0;
    try
       count2 = count;
        while ISdepth2 > bay_Depth % If inlet depth is deeper than bay depth ISdepth(index) will continue counting until it gets to the bay depth
            ISchidx2 = BIidx(end) + count;
            count = count + 1;
            ISdepth2 = z_comb(BIidx(end) + count);
        end
    catch ME
        count = count2;
        diff = 1;
        while diff > 0
            ISchidx2 = BIidx(end) + count;
            count = count + 1;
            ISdepth2 = z_comb(BIidx(end) + count);
            diff = z_comb(ISchidx2) - ISdepth2;
        end
    end
    chidx2 = [OSchidx2:ISchidx2];
    zt_BI2 = z_comb(chidx2);
    zt_bs2 = repelem(bay_depth, length(chidx2));
    delzt2 = zt_BI2-zt_bs2;
    zq_BI2 = repelem(bay_depth,length(chidx2));
end

ch_Lbank = delzt./chn_slope; % left bank y length
ch_Rbank = delzt./chn_slope; % right bank y length

ch_Lbank2 = delzt2./chn_slope; % left bank y length along BI edge
ch_Rbank2 = delzt2./chn_slope; % right bank y length along BI edge

% at each cross-section
ch = cell(numel(zt_BI), 2);
ch2 = cell(numel(zt_BI2),2);
ch3 = ch2;
%Regular inlet cut
for i1 = 1:numel(zt_BI) %Loops from oceanside of the inlet inwards, element by element
    ch_y     = [ 0 ch_Lbank(i1) (chn_bwidth1 + ch_Lbank(i1)) (ch_Rbank(i1) + chn_bwidth1(end) + ch_Lbank(i1)) ]; % y coord from left most end
    ch_z     = [zt_BI(i1) repmat(zt_bs(i1),1,chn_bwidth/yres) zt_BI(i1)]; % depth at each loc
    yq1      = 0:yres:ch_y(end); % interpolation grid to suit the base mesh resolution
    ch_zq    = interp1(ch_y,ch_z,yq1); % interploated depth of chanel at each grid loc
    ch{i1,1} = yq1;
    ch{i1,2} = ch_zq;
end
%Cut at left-most end of grid
for i1 = 1:numel(zt_BI2) % Loops from oceanside of the inlet inwards, element by elemen
    ch_y     = [0 dist2edge1 (ch_Rbank2(i1) + dist2edge) ]; % y coord from left most end
    ch_z     = [zq_BI2(i1) repmat(zq_BI2(i1),1,dist2edge/yres) zt_BI2(i1)]; % depth at each loc
    yq1      = 0:yres:ch_y(end); % interpolation grid to suit the base mesh resolution
    ch_zq    = interp1(ch_y,ch_z,yq1); % interploated depth of chanel at each grid loc
    ch2{i1,1} = yq1;
    ch2{i1,2} = ch_zq;
end
%Cut at right-most end of grid
for i1 = 1:numel(zt_BI2) % Loops from oceanside of the inlet inwards, element by elemen
    ch_y     = [0 ch_Lbank2(i1) (ch_Lbank2(i1) + dist2edge1) ]; % y coord from left most end
    ch_z     = [zt_BI2(i1) repmat(zq_BI2(i1),1,dist2edge/yres) zq_BI2(i1)]; % depth at each loc
    yq1      = 0:yres:ch_y(end); % interpolation grid to suit the base mesh resolution
    ch_zq    = interp1(ch_y,ch_z,yq1); % interploated depth of chanel at each grid loc
    ch3{i1,1} = yq1;
    ch3{i1,2} = ch_zq;
end

% adding channel to partition locations
for j1 = 0:barrier_pati % loop over the number of channels
    %First partition along edge of mesh
    if j1 == 0
        y_chloc = 0;
        [~, yidx] = min(abs(yq -y_chloc));
        for i1 = 1:numel(zt_BI2)
            rowval = chidx2(i1);
            ycell  = length(ch2{i1,1}); % number of cells in the channel at this "rowval"
            ymid   = floor(ycell/2);
            colval = 1:ycell;
            %zmesh_BI(rowval,colval) = ch{i1,2}((end-ymid+1):end);
            zmesh_BI(rowval,colval) = ch2{i1,2};%((end-ymid+1):end);
        end
    %Last partition along edge of mesh    
    elseif j1 == barrier_pati
        y_chloc = len_y;
        [~, yidx] = min(abs(yq -y_chloc));
        for i1 = 1:numel(zt_BI2)
            rowval = chidx2(i1);
            ycell  = length(ch3{i1,1}); % number of cells in the channel at this "rowval"
            ymid   = floor(ycell/2);
            colval = yidx-ycell+1:yidx;
            zmesh_BI(rowval,colval) = ch3{i1,2};%(1:ymid);
        end
    %All other partitions    
    else     
        y_chloc   =  round(((len_y -(2*dist2edge))/barrier_pati)*j1 + dist2edge);
        [~, yidx] = min(abs(yq -y_chloc));
        for i1 = 1:numel(zt_BI)
            rowval = chidx(i1);
            ycell  = length (ch{i1,1}); % number of cells in the channel at this "rowval"
            ymid   = floor(ycell/2);
            colval = yidx-ymid:yidx+ycell-ymid-1;
            zmesh_BI(rowval,colval) = ch{i1,2};
        end
    end
end

% superposition-barrier island - roughness
barrmask = repmat(~idz',1,size(zmesh,2));
idr = zmesh_BI > zthresh ;
rmeshnew1                = rmesh;
rmeshnew1(idr&barrmask)  = r_barr;
zmeshnew1                = zmesh_BI;
zbarrier                 = zq_BI(BIidx);
xbarrier                 = xq(BIidx);
end
