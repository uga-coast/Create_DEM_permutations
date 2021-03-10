function [zmeshnew2, rmeshnew2, xmarsh, zmarsh, znew2] = Func_add_marsh (f_slopem, b_slopem, width, marsh_h, marsh_loc ,r_marsh, r_norm, zmeshnew1, rmeshnew1, base_data)

%zmesh   = base_data.zmesh;
% rmesh   = base_data.rmesh;
xq      = base_data.xq ;
yq      = base_data.yq ;
zq      = base_data.zq;


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
    


end