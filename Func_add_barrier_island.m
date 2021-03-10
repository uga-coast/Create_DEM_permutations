function [xbarrier, zbarrier, zmeshnew1, rmeshnew1] = Func_add_barrier_island (w_req, h_req, xbarrier_loc , barrier_pati, aux_slopeb , chn_bwidth, chn_slope, r_barr, r_norm, zthresh, base_data, datadir,BI_data)

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
h_bi          = Barri.h_bi; %BIdat.h_bi; % original hiehg of the BI


chn_bwidth    = floor(chn_bwidth/yres)*yres;

aspr_h        = h_req/h_bi; % ratio to change the original to suit the required
aspr_w        = w_req/w_bi; % ratio to change the original to suit the required
x2_BI         = x1_BI*aspr_w;
z2_BI         = z1_BI*aspr_h;

[~,MSL_idx]  = min(abs(zq));
xbarr_mid    = xq(MSL_idx) - xbarrier_loc; % X val of the mid point location of the BI
xshft        = x2_BI(Barri.hidx)- xbarr_mid;
% shifted BI coordinates to put the BI at the correct location & height
x3_BI        = x2_BI - xshft;
z3_BI        = z2_BI + Barri.b_basez*aspr_h;

figure
hold on; grid on
plot(xq,zq)
plot(xbarr_mid,0,'o')
plot(x3_BI,z3_BI)

% check if the profiles intersect
L1       = [xq ; zq];
L2       = [x3_BI; z3_BI];
int1     = Func_InterX(L1,L2);

if (isempty(int1))
    % no intersection- needs auxilary slopes to extend until base prof for
    % front and back
    xt1 = linspace(0,x3_BI(1),ceil(x3_BI(1)/0.5));
    zt1 = z3_BI(1) +aux_slopeb*(xt1-xt1(end));
    xt2 = linspace(x3_BI(end),xq(end),ceil((xq(end)-x3_BI(end))/0.5));
    zt2 = z3_BI(end) -aux_slopeb*(xt2-xt2(1));
      
    L1       = [xq; zq];
    L2       = [xt1; zt1];
    int2     = Func_InterX(L1,L2);
    fint     = int2(1,end);
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
        xt1 = linspace(0,x3_BI(1),ceil(x3_BI(1)/0.5));
        zt1 = z3_BI(1) +aux_slopeb*(xt1-xt1(end));
        
        L1       = [xq ; zq];
        L2       = [xt1; zt1];
        int3     = Func_InterX(L1,L2);
        fint     = int3(1,end);
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

x_BI = xt_bi;
z_BI = zt_bi;
figure
hold on; grid on
plot(x_BI,z_BI)
plot(xq,zq)

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
zt_BI = zq_BI(~idz);
zt_bs = zq(~idz);
delzt = zt_BI-zt_bs;
chn_bwidth1 = yres:yres:chn_bwidth-yres ;

% at each cross-section
for i1 = 1:numel(zt_BI)
    ch_Lbank = delzt(i1)/chn_slope; % left bank y length
    ch_Rbank = delzt(i1)/chn_slope; % right bank y length
    ch_y     = cumsum([ 0 ch_Lbank chn_bwidth1 ch_Rbank ]); %  y coord from left most end
    ch_z     = [zt_BI(i1) repmat(zt_bs(i1),1,chn_bwidth/yres) zt_BI(i1)]; % depth at each loc
    yq1      = 0:yres:ch_y(end); % interpolation grid to suit the base mesh resolution
    ch_zq    = interp1(ch_y,ch_z,yq1); % interploated depth of chanel at each grid loc
    ch{i1,1} = yq1;
    ch{i1,2} = ch_zq;
end

% adding channel to partition locations
for j1 = 1:barrier_pati-1 % loop over the number of channels
    y_chloc   =  round(len_y/barrier_pati)*j1;
    [~, yidx] = min(abs(yq -y_chloc));
    
    for i1 = 1:numel(zt_BI)
        rowval = BIidx(i1);
        ycell  = length (ch{i1,1}); % number of cells in the channel at this "rowval"
        ymid   = floor(ycell/2);
        colval = yidx-ymid:yidx+ycell-ymid-1;
        zmesh_BI(rowval,colval) = ch{i1,2};
        
    end
    
end

% superposition-barrier island- roughness
barrmask = repmat(~idz',1,size(zmesh,2));
idr = zmesh_BI > zthresh ;
rmeshnew1                = rmesh;
rmeshnew1(idr&barrmask)  = r_barr;
zmeshnew1                = zmesh_BI;
zbarrier                 = zq_BI(BIidx);
xbarrier                 = xq(BIidx);
end