function [filters, VVs, costFcn] = NewtonMethod4VoraValue(f0, refBs, testBs, maxIter, lam)
         
    N=length(f0);
    fnew = f0;     
    VVs = ones(maxIter+1,1);
    VVs(1) = voraValue(refBs, diag(fnew)*testBs);
    costFcn = ones(maxIter+1,1);
    costFcn(1) = 1-VVs(1) + 0.5*lam*(fnew'*fnew);
    filters = [f0';zeros(maxIter,N)];
    
    step = 10;       alpha = 0.1;     beta = 0.5;   
    [gdf, ntf]= computeDerivatives(fnew, refBs, testBs,lam); 
    x_nt = -inv(ntf)*gdf; % newton step  
    
    iteNo = 0;
    figure; plot(400:10:700, fnew, 'r*-'); hold on;
    while(iteNo <= maxIter)        
        iteNo = iteNo +1;               
        cost = 1- voraValue(refBs, diag(fnew)*testBs)+ 0.5*lam*(fnew'*fnew);
        ftemp = fnew + step*x_nt;
        costtemp = 1- voraValue(refBs, diag(ftemp)*testBs)+0.5*lam*(ftemp'*ftemp);
        while(costtemp > (cost+alpha*step*(gdf'*x_nt))) % backtracking line search for stepsize
            step = step*beta
            ftemp = fnew + step*x_nt;
            costtemp = 1-voraValue(refBs, diag(ftemp)*testBs)+0.5*lam*(ftemp'*ftemp);
        end

        fnew = fnew + step*x_nt;
        VVs(iteNo+1) = voraValue(refBs, diag(fnew)*testBs);
        costFcn(iteNo+1) = 1- VVs(iteNo+1) + 0.5*lam*(fnew'*fnew);
        filters(iteNo,:) = fnew;      
    
        [gdf, ntf]= computeDerivatives(fnew, refBs, testBs,lam); 
        
        x_nt = - inv(ntf)*gdf; % newton step
        plot(400:10:700, fnew./max(fnew)); hold on;
    
    end
    plot(400:10:700, fnew./max(fnew),'g*--'); hold on;
    
    voraValues = 1 - costFcn;
    filters = filters(1:iteNo,:);
    
    figure; 
    subplot 121; plot(VVs,'g*'); title('Vora-Values');
    subplot 122; plot(costFcn,'r'); title('costFcn');


end

function [gdf, ntf]= computeDerivatives(fnew,u, uc, lam) 
    N = length(fnew);
    Px = u*inv(u'*u)*u';
    Pfq = diag(fnew)*uc*inv(uc'*diag(fnew.*fnew)*uc)*uc'* diag(fnew);
    invF = diag(1./fnew);
    In = eye(N);
    
    term1 = Px*Pfq*invF;
    term2 = Pfq*term1;
    gdf = -2*diag(term1 - term2)/3 + lam*fnew; 

    A1= invF*Pfq;
    A2= invF*Pfq*Px*(2*Pfq - 2*In);
    B1= invF*Pfq*invF;
    B2= Px*(In - Pfq);
    C1= invF*Pfq*Px*(2*Pfq - In);
    C2= invF*Pfq;
    D1= -invF*Pfq*Px*Pfq*invF;
    D2= In;
   
    sumTrace = A1'.*A2+B1'.*B2+C1'.*C2+D1'.*D2;
    ntf = -2/3*sumTrace + lam*In;
end

