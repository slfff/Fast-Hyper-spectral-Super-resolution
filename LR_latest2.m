function hyper = LR_latest2(D_hyper, multi)

%%%%%%%%%%%% Lingfei Song 2017.10.31 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast hyperspectral super-resolution using linear regression
% Input: low resolution hyperspectral image 'D_hyper' (data cube);
%        high resolution multispectral image 'multi' (data cube);
% Output: high resolution hyperspectral image 'hyper' (data cube).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(multi, 1);     % CAVE: 512;  Havard: 1024;  ICVL: 1024

addpath('.\SupResPALM\include\');

multi = hyperConvert2d(multi);

D = hyperSpatialDown(sqrt(size(multi, 2)), sqrt(size(multi, 2)), 32);   % downsample operator
D_multi = multi * D;

D_multi = hyperConvert3d(D_multi);
multi = hyperConvert3d(multi);

hyper = zeros(N, N, 31, 'double');

for i = 1:N/64
    for j = 1:N/64
        d_hyper = D_hyper((i-1)*2+1:i*2, (j-1)*2+1:j*2, :);
        d_multi = D_multi((i-1)*2+1:i*2, (j-1)*2+1:j*2, :);
        d_hyper = (hyperConvert2d(d_hyper))';
        d_multi = (hyperConvert2d(d_multi))';
        mean_d_hyper = mean(d_hyper);
        mean_d_multi = mean(d_multi);
        d_hyper = d_hyper - repmat(mean_d_hyper, 4, 1);
        d_multi = d_multi - repmat(mean_d_multi, 4, 1);
        
        Lambda = (d_multi' * d_multi)^(-1) * d_multi' * d_hyper;
        
        tmp_multi = multi((i-1)*2*32+1:i*2*32, (j-1)*2*32+1:j*2*32, :);
        tmp_multi = (hyperConvert2d(tmp_multi))';
        tmp_multi = tmp_multi - repmat(mean(tmp_multi), 2*32*2*32, 1);
        tmp_hyper = tmp_multi * Lambda + repmat(mean_d_hyper, 2*32*2*32, 1);
        
        tmp_hyper = hyperConvert3d(tmp_hyper');
        hyper((i-1)*2*32+1:i*2*32, (j-1)*2*32+1:j*2*32, :) = tmp_hyper;
    end
end
end