function [zmeshnew2, rmeshnew2, xmarsh, zmarsh] = Func_add_marsh (f_slopem, b_slopem, width, marsh_h, marsh_loc ,r_marsh, r_norm, zmeshnew1, rmeshnew1, base_data)

%zmesh   = base_data.zmesh;
% rmesh   = base_data.rmesh;
xq      = base_data.xq ;
yq      = base_data.yq ;
zqint   = base_data.zq;
zq      = zmeshnew1(:,1);
zq      = reshape(zq,1,length(zq));
zqmall  = -999*ones(length(zq),1); % to be filled


[~,MSL_idx]  = min(abs(zqint)); % location of MSL
xmars_mid    = xq(MSL_idx) - marsh_loc; % mid location of marsh
xmars_frt    = xmars_mid-width/2;  % Front of marsh
xmars_bac    = xmars_mid+width/2;  % back of marsh


[~,mfrt_idx] = min(abs(xq-xmars_frt)); % index of front of marsh
[~,mbac_idx] = min(abs(xq-xmars_bac)); % index of back of marsh

% construct front marsh slope
xtempf       = xq(1:mfrt_idx);
ztempf       = zq(1:mfrt_idx);
zm_frt       = marsh_h - (-xtempf +xq(mfrt_idx))*f_slopem;

% find the intersection point of the slope with profile
L1       = [xtempf ; ztempf];
L2       = [xtempf; zm_frt];
int1     = Func_InterX(L1,L2);

bint      = int1(1,end); % last intersection 
[~, lfx]  = min(abs(xq-bint)); % intersection index in global coord
[~, lfxt] = min(abs(xtempf-bint));  % intersection index in temp coord

% construct back marsh slope
xtempf       = xq(mbac_idx:end);
ztempf       = zq(mbac_idx:end);
zm_bac       = marsh_h - (xtempf -xq(mbac_idx))*b_slopem;

% find the intersection point of the back slope with profile

if(mbac_idx > MSL_idx) % marsh back is behind 0 MSL line. back width to 0 MSL
    
    mbac_idx = MSL_idx;
    lbx      = MSL_idx;
else
    L1       = [xtempf ; ztempf];
    L2       = [xtempf; zm_bac];
    int2     = Func_InterX(L1,L2);
    
    bint     = int2(1,1); % first intersection
    [~, lbx] = min(abs(xq-bint));
    [~, lbxt]= min(abs(xtempf-bint));  % intersection index in temp coord
    zqmall(mbac_idx:lbx) = zm_bac(1:lbxt);
end

zqmall(lfx:mfrt_idx)          = zm_frt(lfxt:end);
zqmall(mfrt_idx+1:mbac_idx)   = marsh_h;


zqmall2 = repmat(zqmall,1,length(yq));

idmar             = zqmall2 ~= -999; 
zmeshnew2         = zeros(size(zmeshnew1,1),size(zmeshnew1,2));
zmeshnew2(idmar)  = zqmall2(idmar);
zmeshnew2(~idmar) = zmeshnew1(~idmar);


xmarsh = xq(lfx:lbx);
zmarsh = zqmall(lfx:lbx);

%superposition-marsh- roughness
for j1 = 1:numel(yq)
    for i1 =1:numel(xq)
        idx = xq(i1)== xmarsh;
        if any(idx)
            rmeshnew2(i1,j1) = r_marsh;
        else
            rmeshnew2(i1,j1) = rmeshnew1(i1,j1);
        end
    end
end



end