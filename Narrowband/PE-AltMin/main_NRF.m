clear,clc

load('Ns=6.mat');
% load('Ns=8.mat');
Ns = 6;

NRF = 6:11;

SNR_dB = 0;
SNR = 10.^(SNR_dB./10);
realization = size(H,3);

for r = 1:length(NRF)
    parfor reali = 1:realization
        [ FRF, FBB ] = PE_AltMin( Fopt(:,:,reali), NRF(r));
        FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
        [ WRF, WBB ] = PE_AltMin( Wopt(:,:,reali), NRF(r));
        R(r,reali) = log2(det(eye(Ns) + SNR/Ns * pinv(WRF * WBB) * H(:,:,reali) * FRF * FBB * FBB' * FRF' * H(:,:,reali)' * WRF * WBB));
    end
end
plot(NRF,sum(R,2)/realization,'Marker','>','LineWidth',1.5,'Color',[0 0.447058826684952 0.74117648601532]);
grid on
hold on