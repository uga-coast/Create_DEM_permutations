function [base_data] = Func_setup_base_profile(len_y, xres, yres, r_norm, r_base, zthresh, datadir)

duck_base   = 0; 
panama_base = 1; 
if (duck_base ==1)
    load(fullfile(datadir,'profiles_Duck'));
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
    
elseif (panama_base ==1)
    load(fullfile(datadir,'panama_base_beach.mat'));
    zq  = z_vals; 
    idx = xq<4500;
    xq  = xq(idx);
    zq  = zq(idx);
end

yq              = linspace(0,len_y,ceil(len_y/yres)+1);
[ymesh,xmesh]   = meshgrid(yq,xq);
zmesh           = repmat(zq',1,length(yq));


rmesh           = repmat(r_norm,numel(xq),numel(yq));
% roughness mask for base profile
base_mask       = zmesh > zthresh;
rmesh(base_mask)=  r_base;

base_data.xmesh = xmesh;
base_data.ymesh = ymesh;
base_data.zmesh = zmesh;
base_data.rmesh = rmesh;
base_data.xq    = xq;
base_data.yq    = yq;
base_data.zq    = zq;
end