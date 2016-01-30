clear,clc

addpath(pwd);
cd manopt;
addpath(genpath(pwd));
cd ..;

load('Ns=6.mat');
% load('Ns=8.mat');
Ns = 6;

NRF = 10;

SNR_dB = 0;
SNR = 10.^(SNR_dB./10);
realization = size(H,3);

R = zeros(1,length(NRF));

for r = 1:length(NRF)
    parfor reali = 1:realization
        [ FRF, FBB ] = MO_AltMin( Fopt(:,:,reali), NRF(r) );
        FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
        [ WRF, WBB ] = MO_AltMin( Wopt(:,:,reali), NRF(r) );
        R(r,reali) = log2(det(eye(Ns) + SNR/Ns * pinv(WRF * WBB) * H(:,:,reali) * FRF * FBB * FBB' * FRF' * H(:,:,reali)' * WRF * WBB));
    end
end
plot(NRF,sum(R,2)/realization,'m-p','LineWidth',1.5)
grid on
hold on
