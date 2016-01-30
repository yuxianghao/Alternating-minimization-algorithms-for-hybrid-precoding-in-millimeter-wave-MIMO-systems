clear,clc

addpath(pwd);
cd manopt;
addpath(genpath(pwd));
cd ..;

% load('Ns=3.mat');
load('Ns=3.mat');

Ns = 3;

NRF = 3;

SNR_dB = -35:5:5;
SNR = 10.^(SNR_dB./10);
realization = size(H,3);
smax = length(SNR);% enable the parallel

for reali = 1:realization
    [ FRF, FBB ] = MO_AltMin( Fopt(:,:,reali), NRF);
    FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
    [ WRF, WBB ] = MO_AltMin( Wopt(:,:,reali), NRF);

    for s = 1:smax
        R(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRF * WBB) * H(:,:,reali) * FRF * FBB * FBB' * FRF' * H(:,:,reali)' * WRF * WBB));
    end
end
plot(SNR_dB,sum(R,2)/realization,'k-p','LineWidth',1.5)
grid on
hold on