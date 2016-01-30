clear,clc

addpath(pwd);
cd cvx;
addpath(genpath(pwd));
cd ..;

load('Ns=3.mat');

Ns = 3;

NRF = 3;

SNR_dB = -35:5:5;
SNR = 10.^(SNR_dB./10);
realization = size(H,3);
smax = length(SNR);% enable the parallel

parfor reali = 1:realization
    [ FRF, FBB ] = SDR_AltMin( Fopt(:,:,reali), NRF);
    [ WRF, WBB ] = Receiver( Wopt(:,:,reali), NRF);
    for s = 1:smax
        R(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRF * WBB) * H(:,:,reali) * FRF * FBB * FBB' * FRF' * H(:,:,reali)' * WRF * WBB));
    end
end
plot(SNR_dB,sum(R,2)/realization,'Marker','diamond','LineWidth',1.5,'Color',[0.87058824300766 0.490196079015732 0]);
grid on
hold on