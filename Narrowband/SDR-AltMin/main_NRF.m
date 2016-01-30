clear,clc

addpath(pwd);
cd cvx;
addpath(genpath(pwd));
cd ..;

load('Ns=2.mat');
% load('Ns=8.mat');

Ns = 2;

NRF = [2,3,4,6,9,12,18];

SNR_dB = 0;
SNR = 10.^(SNR_dB./10);
realization = size(H,3);
smax = length(SNR);% enable the parallel

for r = 1:length(NRF)
    for reali = 1:realization
        [ FRF, FBB ] = SDR_AltMin( Fopt(:,:,reali), NRF(r) );
        [ WRF, WBB ] = Receiver( Wopt(:,:,reali), NRF(r) );
        R(r,reali) = log2(det(eye(Ns) + SNR/Ns * pinv(WRF * WBB) * H(:,:,reali) * FRF * FBB * FBB' * FRF' * H(:,:,reali)' * WRF * WBB));    
    end
end
plot(NRF,sum(R,2)/realization,'Marker','diamond','LineWidth',1.5,'Color',[0.87058824300766 0.490196079015732 0]);
grid on
hold on