function [base_data, xbarr_mid, bay_check] = Func_setup_base_profile(ptype,barrier_island,mainland,hmax_ml,...
    mlslope,back_bay,len_x,DOC,DOCSlope,OSSlope,bay_Depth,len_Bay,startDepth,...
    w_req,len_y, xres, yres, r_norm, r_base, zthresh, datadir, bay_override, ...
    DOC_len, halfBay, deepbay_len, OSoverride)
xbarr_mid = NaN;
%Duck Base Profile
if (ptype == 0)
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

%Tyndall Base Profile    
elseif (ptype == 1)
    load(fullfile(datadir,'panama_base_beach.mat'));
    zq      = z_vals; 
    idx     = xq < len_x;
    xq1     = xq(idx);
    zq1     = zq(idx);
    len_x   = xq(end);
    xqint   = linspace(xq1(1),xq1(end),ceil(len_x/xres)+1);    
    zq      = interp1(xq1,zq1,xqint,'spline');    
    xq      = xqint;
    
    
%Generating Idealized Base profile    
elseif (ptype == 2)
    xq = linspace(0,len_x,round(len_x/xres)+1);
    zq = zeros(1,round(len_x/xres)+1);
    i = 0;
    l = 1;
    %Creates mainland section of profile
    if mainland == 1
        zi = hmax_ml;
        i = i + 1;
        while zi > startDepth
            zq(i) = hmax_ml - (i-1)*mlslope*xres;
            zi  = zq(i);
            i = i + 1;
        end
        i = i-1;
    end
    %For barrier island profile with back bay
    if (back_bay == 1) %&& (bay_override == 0)
        len = round(len_Bay/xres)+1;
        bay_check = (bay_Depth > DOC);
        j = 1;
        %If depth of bay is less than depth of closure
        if bay_check == 1
            if mainland == 0
                while (zi - DOCSlope*xres) > bay_Depth
                    zq(i+j) = startDepth - (j-1)*DOCSlope*xres;
                    zi = zq(i+j);
                    zdown(l) = zi;
                    l = l+1;
                    j = j+1;
                end
            else
                while (zi - DOCSlope*xres) > bay_Depth 
                    zq(i+j) = zq(i+j-1) - j*DOCSlope*xres;
                    zi = zq(i+j);
                    if zi <= 0
                        zdown(l) = zi;
                        l = l+1;
                    end
                    j = j+1;
                end
            end
            l = l-1;
            idown = i + floor((abs(bay_Depth)/DOCSlope)/xres) + 1;
            if barrier_island == 1
                iup = len - (idown);
                for k = idown:(iup)
                    zq(k) = bay_Depth;
                    j = j+1;
                end
                zup = fliplr(zdown);
                zq((i+j):(i+j+l-1)) = zup;
                j = j+l;
                %while zi + OSSlope*xres >= DOC && (zi + DOCSlope*xres) <= 0
                 %   zq(i+j) = zq((i+j)-1) + DOCSlope*xres;
                 %   zi = zq(i+j);
                 %   j = j + 1;
                %end
                for k = 1:(ceil(w_req/xres)+1)
                    zq(i+j) = 0;
                    zi = zq(i+j);
                    j = j+1;
                    if k == round(ceil(w_req/xres)*0.5)
                        xbarr_mid = len_x - (i+j)*xres + xres;
                    end
                end
            else
                for k = idown:len
                    zq(k) = bay_Depth;
                    j = j+1;
                end
            end
            while zi > DOC
                zq(i+j) = zq((i+j)-1) - DOCSlope*xres;
                zi = zq(i+j);
                j = j + 1;
            end
            while zi <= DOC && (i+j) <= length(zq)
                zq(i+j) = zq((i+j)-1) - OSSlope*xres;
                zi = zq(i+j);
                j = j + 1;
            end
            zq = fliplr(zq);
        elseif bay_check == 0 %If bay depth is greater than depth of closure
            iDownDOC = round((abs(DOC)/DOCSlope)/xres)+1;
            iDownBay = round(abs(bay_Depth - DOC)/OSSlope/xres)+1;
            iupDOC = len - iDownDOC + 1;
            iupBay = len - iDownBay - iDownDOC;
            zi = startDepth;
            if mainland == 0
                while zi > DOC && (l -1) < ceil((len_Bay/yres)/2)
                    zq(i+j) = startDepth - (j-1)*DOCSlope*xres;
                    zi = zq(i+j);
                    zdown(l) = zi;
                    l = l+1;
                    j = j + 1;
                end
            else
                while zi > DOC && (l -1) < ceil((len_Bay/yres)/2)
                    zq(i+j) = zq(i+j-1) - DOCSlope*xres;
                    zi = zq(i+j);
                    if zi <= 0
                        zdown(l) = zi;
                        l = l+1;
                    end
                    j = j + 1;
                end
            end
            %iupDOC = len - (j - 1);
            while zi <= DOC && zi - OSSlope*xres >= bay_Depth %&& (l -1) < ceil((len_Bay/yres)/2)
                zq(i+j) = zq((i+j)-1) - OSSlope*xres;
                zi = zq(i+j);
                zdown(l) = zi;
                l = l+1;
                j = j + 1;
            end
            l = l -1;
            %iupBay = len - (j - 1);
            if barrier_island == 1
                for k = (j):(iupBay)
                    zq(k+i) = bay_Depth;
                    j = j+1;
                end
                %o = o-1;
                %zq(i+j:i+j+m+n) = fliplr(zq(i+j-o-m-n:i+j-o)); 
                zup = fliplr(zdown);
                zq((i+j):(i+j+l-1)) = zup;
                j = j+l;
                while (zi + OSSlope*xres) <= DOC
                    zq(i+j) = zq((i+j)-1) + OSSlope*xres;
                    zi = zq(i+j);
                    j = j + 1;
                end 
                while zi + OSSlope*xres >= DOC && (zi + DOCSlope*xres) <= 0 %(len - (iupDOC-1)):(len - iupBay-1)
                    zq(i+j) = zq((i+j)-1) + DOCSlope*xres; %zq(len - i) = zq(iupDOC) - (i-(len - (iupDOC)))*OSSlope*xres;
                    zi = zq(i+j);
                    j = j + 1;
                end
                for k = 1:(ceil(w_req/xres+1))
                    zq(i+j) = 0;
                    zi = zq(i+j);
                    if k == floor(ceil(w_req/xres+1)*0.5)
                        xbarr_mid = len_x - (i+j)*xres;
                    end
                    j = j+1;
                end
            else
                for k = (j):len
                    zq(k+i) = bay_Depth;
                    j = j+1;
                end
            end
            while zi - DOCSlope*xres > DOC
                zq(i+j) = zq((i+j)-1) - DOCSlope*xres;
                zi = zq(i+j);
                j = j + 1;
            end
            while (i+j) <= length(zq)
                zq(i+j) = zq((i+j)-1) - OSSlope*xres;
                zi = zq(i+j);
                j = j + 1;
            end
            zq = fliplr(zq);
        else
            error("Something went wrong. Check DOC and bay_Depth.")
        end 
        
%% Experimental        
    %For profile without barrier island/back_bay  
    %{
    elseif (back_bay == 1) && (bay_override == 1)
        len = round(len_Bay/xres)+1;
        bay_check = (bay_Depth > DOC);
        j = 1;
        %If depth of bay is less than depth of closure
        if bay_check == 1 || (bay_check == 0 && OSoverride == 1)
            pos = xres;
            if mainland == 0
                while pos < halfBay
                    zq(i+j) = startDepth - (j-1)*DOCSlope*xres;
                    pos = pos + xres;
                    j = j+1;
                end
            else
                while pos < halfBay 
                    zq(i+j) = zq(i+j-1) - j*DOCSlope*xres;
                    pos = pos + xres;
                    j = j+1;
                end
            end
            if barrier_island == 1
                while pos <= len_Bay
                    zq(i+j) = zq(i+j-1) + j*DOCSlope*xres;
                    pos = pos + xres;
                    j = j+1;
                end
                for k = 1:(ceil(w_req/xres)+1)
                    zq(i+j) = 0;
                    zi = zq(i+j);
                    j = j+1;
                    if k == round(ceil(w_req/xres)*0.5)
                        xbarr_mid = len_x - (i+j)*xres + xres;
                    end
                end
            end
            while zi > DOC
                zq(i+j) = zq((i+j)-1) - DOCSlope*xres;
                zi = zq(i+j);
                j = j + 1;
            end
            while zi <= DOC && (i+j) <= length(zq)
                zq(i+j) = zq((i+j)-1) - OSSlope*xres;
                zi = zq(i+j);
                j = j + 1;
            end
            zq = fliplr(zq);
        elseif bay_check == 1  %If bay depth is greater than depth of closure
            iDownDOC = round((abs(DOC)/DOCSlope)/xres)+1;
            iDownBay = round(abs(bay_Depth - DOC)/OSSlope/xres)+1;
            iupDOC = len - iDownDOC + 1;
            iupBay = len - iDownBay - iDownDOC;
            zi = startDepth;
            if mainland == 0
                while zi > DOC && (l -1) < ceil((len_Bay/yres)/2)
                    zq(i+j) = startDepth - (j-1)*DOCSlope*xres;
                    zi = zq(i+j);
                    zdown(l) = zi;
                    l = l+1;
                    j = j + 1;
                end
            else
                while zi > DOC && (l -1) < ceil((len_Bay/yres)/2)
                    zq(i+j) = zq(i+j-1) - DOCSlope*xres;
                    zi = zq(i+j);
                    if zi <= 0
                        zdown(l) = zi;
                        l = l+1;
                    end
                    j = j + 1;
                end
            end
            %iupDOC = len - (j - 1);
            while zi <= DOC && zi - OSSlope*xres >= bay_Depth && (l -1) < ceil((halfBay/yres)/2)
                zq(i+j) = zq((i+j)-1) - OSSlope*xres;
                zi = zq(i+j);
                zdown(l) = zi;
                l = l+1;
                j = j + 1;
            end
            l = l -1;
            %iupBay = len - (j - 1);
            if barrier_island == 1
                zup = fliplr(zdown);
                zq((i+j):(i+j+l-1)) = zup;
                j = j+l;
                for k = 1:(ceil(w_req/xres+1))
                    zq(i+j) = 0;
                    zi = zq(i+j);
                    if k == floor(ceil(w_req/xres+1)*0.5)
                        xbarr_mid = len_x - (i+j)*xres;
                    end
                    j = j+1;
                end
            else
                for k = (j):len
                    zq(k+i) = bay_Depth;
                    j = j+1;
                end
            end
            while zi - DOCSlope*xres > DOC
                zq(i+j) = zq((i+j)-1) - DOCSlope*xres;
                zi = zq(i+j);
                j = j + 1;
            end
            while (i+j) <= length(zq)
                zq(i+j) = zq((i+j)-1) - OSSlope*xres;
                zi = zq(i+j);
                j = j + 1;
            end
            zq = fliplr(zq);
        else
            error("Something went wrong. Check DOC and bay_Depth.")
        end 
        %}
    elseif(back_bay == 0)
        %zi = 0;
        bay_check = 0;
        zi = startDepth;
        for j = round(xq(1:(end-i))/xres) + 1
            if mainland == 0
                if zi >= (DOC + .1)
                    zq(i+j) = startDepth - (j-1)*DOCSlope*xres;
                    zi = zq(i+j);  
                else
                    zq(i+j) = DOC - (j-1)*OSSlope*xres;
                    zi = zq(i+j);
                end
            else
                if zi >= (DOC + .1)
                    zq(i+j) = zq(i+j-1) - j*DOCSlope*xres;
                    zi = zq(i+j);  
                else
                    zq(i+j) = DOC - (j-1)*OSSlope*xres;
                    zi = zq(i+j);
                end
            end
        end
        zq = fliplr(zq);
    else
        error("Please input 0 or 1 for back_bay variable")
    end
end

yq              = linspace(0,len_y,ceil(len_y/yres)+1);
[ymesh,xmesh]   = meshgrid(yq,xq);
zmesh           = repmat(zq',1,length(yq));

if ptype == 0 || ptype == 1
    bay_check = NaN;
end

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
