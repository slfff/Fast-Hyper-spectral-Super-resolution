datafiles=dir('..\data\HSI_CAVE\*.mat');  % GroundTruth
load('.\CRF_database\NIKON_D700.mat');   % crf
addpath('.\SupResPALM\include');

RMSE = zeros(1,size(datafiles,1));
SAM = zeros(1,size(datafiles,1));
time = zeros(1,size(datafiles,1));

for i = 1:size(datafiles,1)
    load(strcat('..\data\HSI_CAVE\', datafiles(i).name));
    
    GroundTruth = hyperConvert2d(GroundTruth);
    GroundTruth = GroundTruth ./ max(GroundTruth(:));
    [D_GroundTruth, RGB_GroundTruth] = hyperSynthetic(GroundTruth, crf, 32);

    D_GroundTruth = hyperConvert3d(D_GroundTruth);
    RGB_GroundTruth = hyperConvert3d(RGB_GroundTruth);

    tic;
    hyper = LR_latest4(D_GroundTruth, RGB_GroundTruth); %%%%%%
    time(i) = toc;
    hyper(hyper > 1) = 1;
    hyper(hyper < 0) = 0;

    hyper = round(hyper * 255);
    GroundTruth = round(GroundTruth * 255);

    RMSE(i) = hyperErrRMSE(GroundTruth, hyper);
    SAM(i) = hyperErrSam(GroundTruth, hyper);
    fprintf(datafiles(i).name);
    fprintf('  %f  %f  %f\n', RMSE(i), SAM(i), time(i));

end

fprintf('aver RMSE %f\n aver SAM %f\n, aver time %f\n', mean(RMSE), mean(SAM), mean(time));