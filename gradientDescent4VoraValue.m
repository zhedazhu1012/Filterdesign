function [filters, voraValues] = gradientDescent4VoraValue(f0, refBs, testBs, maxIter )
    % implement by Yuteng Zhu @UEA
    
    N=length(f0);
    fnew = f0;
    costFcn = ones(maxIter+1,1);
    costFcn(1) = 1-voraValue(refBs, diag(fnew)*testBs);
    filters = zeros(maxIter,N);

    gdf = computeGDFcnNew(fnew,refBs, testBs);  

    iteNo = 0;
    step = 10;   alpha = 0.1;     beta = 0.5; % for line search

    figure
    plot(400:10:700, fnew); hold on;
    while( iteNo < maxIter && max(abs(gdf))> 1e-10)   
        iteNo = iteNo+1;
        cost = voraValue(refBs, diag(fnew)*testBs);
        ftemp = fnew + step*gdf;
        costtemp = voraValue(refBs, diag(ftemp)*testBs);
        while(costtemp < (cost+alpha*step*(gdf'*gdf))) % backtracking line search for stepsize
            step = step*beta;
            ftemp = fnew + step*gdf;
            costtemp = voraValue(refBs, diag(ftemp)*testBs);
        end
        fnew = fnew + step*gdf;
        costFcn(iteNo+1) = 1-voraValue(refBs, diag(fnew)*testBs);
        filters(iteNo,:) = fnew;
        gdf = computeGDFcnNew(fnew,refBs, testBs);  

        plot(400:10:700, fnew./max(fnew)); hold on;
    end
    voraValues = 1- costFcn;
    filters = filters(1:iteNo,:);
end

function gdf = computeGDFcnNew(fnew,u, uc)
    
    Pfuc = diag(fnew)*uc*inv(uc'*diag(fnew.*fnew)*uc)*uc'*diag(fnew);
    Pu = u*inv(u'*u)*u';
    term1 = Pu*Pfuc*diag(1./fnew);
    term2 = Pfuc*term1;
    gdf = 2*diag(term1 - term2);
end

