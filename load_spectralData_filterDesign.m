function [cmf, sensor, refSFU, illums, illumD65,illumA, wavel]=load_spectralData_filterDesign(cameraNo, dataSetting)

if nargin <2
    dataSetting =1;
end

camFolder = 'camera spectral data/';
files=dir(camFolder);
load('1931_CMF_380to780_5.mat'); 
load('SFU_illuminants_normalized_380to780by4nm.mat'); 
load('cie_illuminants_Smet.mat', 'D65', 'IllA');
load('SFU_reflectance_db.mat');

illum_SFU_interp = interp1(380:4:780, illum_SFU_normlized, 380:1:780, 'pchip');
refl = interp1(380:4:780,reflectance_SFU,380:1:780, 'pchip'); 

if (dataSetting ==1)
    % the default setting    
    load([camFolder files(cameraNo).name]);        
    wavel = 400:10:700;
    N = length(wavel);
    rgb_new = rgb(1:31,:);
    scaling = 1/max(rgb_new(:,2));
    sensor = scaling* rgb_new./repmat(sum(rgb_new), length(wavel),1); % ration correction: sum up to 1 for each channel
    cmf = CMF_2(5:2:65, 2:4); 
    
    ills_truc = illum_SFU_interp(21:10:321, :);
    ills_trucnorm = ills_truc./repmat(max(ills_truc),N,1);
    illums = ills_trucnorm;

    illumD65 = D65(21:10:321,2); 
    illumA = IllA(21:10:321,2);    
    refSFU = refl(21:10:321,:);
    
elseif (dataSetting ==2)
    % take smaller wavelength range as the two ends are small numbers    
    load([folder files(cameraNo).name]);  
    wavel = 420:10:700;
    N = length(wavel); % suppose the data is arranged as N*3
    cmf = CMF_2(9:2:65, 2:4); 
    cmf = cmf./repmat(sum(cmf), N,1);
    rgb_new = rgb(3:31,:);
    sensor = rgb_new./repmat(sum(rgb_new), N,1); 
    illums = illum_SFU_interp(41:10:321, :);
    illumD65 = D65(41:10:321,2);
    illumA = IllA(41:10:321,2);
    refSFU = refl(41:10:321,:);
elseif (dataSetting ==3)
    %% interpret data into smaller intervals    
    load([folder files(cameraNo).name]);        
    wavel=400:5:700;
    N = length(wavel); % suppose the data is arranged as N*3
    cmf = CMF_2(5:65, 2:4); 
    cmf = cmf./repmat(sum(cmf),N,1);
    rgb_new=interp1(400:10:700, rgb(1:31,:), 400:5:700, 'pchip');
    sensor = rgb_new./repmat(sum(rgb_new), N,1); % ration correction: sum up to 1 for each channel
    illums = illum_SFU_interp(21:5:321, :);
    illumD65 = D65(21:5:321,2);
    illumA = IllA(21:5:321,2);
    refSFU = refl(21:5:321,:);
else
   error('Plase enter an approariate setting'); 
end

end

