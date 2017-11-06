datafiles=dir('..\data\HSI_Havard\*.mat');
load('.\CRF_database\NIKON_D70.mat');
addpath('.\SupResPALM\include');

RMSE = zeros(1,size(datafiles,1));
SAM = zeros(1,size(datafiles,1));
time = zeros(1,size(datafiles,1));

for i = 1:size(datafiles,1)
    load(strcat('..\data\HSI_Havard\', datafiles(i).name));
    
    ref = ref(1:1024, 1:1024, :);
    ref = hyperConvert2d(ref);
    
    ref = ref ./ max(ref(:));
    [D_ref, RGB_ref] = hyperSynthetic(ref, crf, 32);
    D_ref = hyperConvert3d(D_ref);
    RGB_ref = hyperConvert3d(RGB_ref);
    
    tic;
    hyper = LR_latest4(D_ref, RGB_ref); %%%%%%
    time(i) = toc;
    hyper(hyper > 1) = 1;
    hyper(hyper < 0) = 0;

    hyper = round(hyper * 255);
    ref = round(ref * 255);

    RMSE(i) = hyperErrRMSE(ref, hyper);
    SAM(i) = hyperErrSam(ref, hyper);
    fprintf(datafiles(i).name);
    fprintf('  %f  %f  %f\n', RMSE(i), SAM(i), time(i));

end

fprintf('aver RMSE %f\n aver SAM %f\n, aver time %f\n', mean(RMSE), mean(SAM), mean(time));