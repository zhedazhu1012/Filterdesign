function [results] = dataDriven_filterDesign_cosineBasis(css,cmf, refl,illum, filtInit, iteNo, methods)

    [p, K] = size(illum);
    [~, M] = size(refl);    
    
    filtBasis = methods{5};
    retfilter = filtInit; 
    xyzw = zeros(K, 3);
    ALLXYZ = zeros(M, 3, K);
    CSset = zeros(M, p, K);
    Lref = zeros(M, 3, K);
    ccMatrix = repmat(eye(3,3),1,1,K); 
    err = zeros(K, M);
    vecXYZ = zeros(M*3*K,1); 
    for k = 1:K
        CS = refl'*diag(illum(:, k));
        CSset(:,:,k) = CS;        
        RGBftemp= CS*retfilter*css;   
        xyzw(k, :) =  illum(:, k)'* cmf;
        xyztemp = CS*cmf;
        ALLXYZ(:,:,k) = xyztemp;              
        vecXYZ(1+M*3*(k-1):M*3*k,1) = xyztemp(:);
        Lref(:,:,k)=xyz2lab(ALLXYZ(:,:,k),xyzw(k, :));        
                
        curM =  inv(RGBftemp'*RGBftemp)*RGBftemp'*xyztemp;
        ccMatrix(:,:,k) = curM;      
        XYZ_est = CS*retfilter*css*curM;     
        lab1 = xyz2lab(XYZ_est, xyzw(k, :));
        err(k, :) = sqrt((squeeze(Lref(:,:,k)) - lab1).^2*ones(3,1));
    end
  
    meanXYZDE= zeros(iteNo,1);
    meanDE=zeros(iteNo,1);
    meanDEmin= 1000;    
    RGBS = zeros(M,3,p,K);
     for j= 1:p
            D0 = zeros(p);
            D0(j,j) = 1;
            for k = 1:K                 
                 RGBS(:,:,j,k)=squeeze(CSset(:,:,k))*D0*css;
            end
     end 
    
    for i = 1:iteNo          
         V = zeros(K*M*3,p);
         for j= 1:p
            V0 = zeros(K*M*3,1);
            for k = 1:K 
                rgb = squeeze(RGBS(:,:,j,k))*squeeze(ccMatrix(:,:,k));
                V0(1+M*3*(k-1):M*3*k,1) = rgb(:);
            end
            V(:,j) = V0; 
         end 
        V=V*filtBasis;   
        H = V'*V;
        f = -V'*vecXYZ;
        A = [-filtBasis; filtBasis];
        b = [-methods{4}(1)*ones(p,1); methods{4}(2)*ones(p,1)];
        [filtCoef, ~, ~, ~]= quadprog(H, f,A,b);         
        retfilter = diag(filtBasis*filtCoef);        
       
        err = zeros(K,M);
        errXYZ = zeros(K,M);
        ALLXYZest=zeros(M,3,K);  
        for k = 1:K
            RGBftemp = squeeze(CSset(:, :, k))*retfilter*css;
            curM =  inv(RGBftemp'*RGBftemp)*RGBftemp'*squeeze(ALLXYZ(:,:,k));           
            ccMatrix(:,:,k) = curM;      
            XYZ_est = squeeze(CSset(:,:,k))*retfilter*css*curM;     
            ALLXYZest(:,:,k)=XYZ_est;
            lab1 = xyz2lab(XYZ_est, xyzw(k, :));
            errXYZ(k,:) = sqrt((XYZ_est - squeeze(ALLXYZ(:,:,k))).^2*ones(3,1));
            err(k, :) = sqrt((squeeze(Lref(:,:,k)) - lab1).^2*ones(3,1));
        end                               
        meanDE(i) =mean(err(:));
        meanXYZDE(i)=mean(errXYZ(:));
        if (meanDE(i) < meanDEmin(1))
            idx = i;
            meanDEmin = [mean(err(:)), median(err(:)),prctile(err(:), 90), prctile(err(:), 95), prctile(err(:), 99), max(err(:))]; 
            filter=diag(retfilter); 
        end
    end    
  
    results.meanErrors = meanDEmin;
    results.filter=filter;
    results.iteNoIdx = idx;
end

