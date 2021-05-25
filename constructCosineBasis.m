
% corrected on 10 Sept 2019 by YT
% as previously it did not deal with the sqrt(2) at the two ends
% and there is only 'wrongly' implemented version 1

function filterBasis = constructCosineBasis(N, basisNo, version)
    if (nargin <2)
        basisNo=N;
        version='2';
    elseif (nargin <3)
        version='2';
    end
    
    filterBasis = zeros(N,basisNo);
    switch version
        case '1'
            filterBasis = dct1(N, basisNo);
        case '2'
            filterBasis = dct2(N, basisNo);
        case '3'
            filterBasis = dct3(N, basisNo);
        case '4'
            filterBasis = dct4(N, basisNo);
        otherwise
            error('Please select the DCT type');
    end

end

function filterBasis = dct1(N, basisNo)
    %% DCT-1
    filterBasis = zeros(N,basisNo);
    for k = 1:basisNo
        for j=1:N
            delta_j1=0; 
            delta_jn=0;
            delta_k1=0;
            delta_kn=0;
            if(j==1)
                delta_j1=1;
            elseif (j==N)
                delta_jn=1;
            end
            if(k==1)
                delta_k1=1;
            elseif (k==N)
                delta_kn=1;
            end
            filterBasis(j, k) = sqrt(2/(N-1))*cos((k-1)*(j-1)*pi/(N-1))/sqrt(1+delta_j1+delta_jn)/sqrt(1+delta_k1+delta_kn);
        end
    end
end

function filterBasis = dct2(N, basisNo)
    %% DCT-2
    filterBasis = zeros(N,basisNo);
    for k = 1:basisNo
        for j=1:N
            filterBasis(j, k) = sqrt(2/N)*cos((k-1)*(2*j-1)*pi/2/N);
            if(k==1)
                filterBasis(j, k) = filterBasis(j, k)./sqrt(2);
            end
        end
    end
end

function filterBasis = dct3(N, basisNo)
    %% DCT-3
    filterBasis = zeros(N,basisNo);
    for k = 1:basisNo
        for j=1:N
            filterBasis(j, k) = sqrt(2/N)*cos((2*k-1)*(j-1)*pi/2/N);
            if(j==1)
                filterBasis(j, k) = filterBasis(j, k)./sqrt(2);
            end
        end
    end
end

function filterBasis = dct4(N, basisNo)
    %% DCT-4
    filterBasis = zeros(N,basisNo);
    for k = 1:basisNo
        for j=1:N
            filterBasis(j, k) = sqrt(2/N)*cos((2*k-1)*(2*j-1)*pi/4/N);
        end
    end
end