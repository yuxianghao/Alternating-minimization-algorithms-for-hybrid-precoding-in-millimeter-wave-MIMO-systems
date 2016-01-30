clear,clc

load('Ns=3.mat');
% load('Ns=8, 2000.mat');
Ns = 3;

NRF = 3;

SNR_dB = -35:5:5;
SNR = 10.^(SNR_dB./10);
realization = size(H,3);
smax = length(SNR);% enable the parallel

for reali = 1:realization
    [ FRF, FBB ] = OMP( Fopt(:,:,reali), NRF, At(:,:,reali) );
    FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
    [ WRF, WBB ] = OMP( Wopt(:,:,reali), NRF, Ar(:,:,reali) );
    for s = 1:smax
        R(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRF * WBB) * H(:,:,reali) * FRF * FBB * FBB' * FRF' * H(:,:,reali)' * WRF * WBB));
        R_o(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(Wopt(:,:,reali)) * H(:,:,reali) * Fopt(:,:,reali) * Fopt(:,:,reali)' * H(:,:,reali)' * Wopt(:,:,reali)));
    end
end
figure()
plot(SNR_dB,sum(R,2)/realization,'Marker','^','LineWidth',1.5,'Color',[0 0.498039215803146 0])
grid on
hold on
plot(SNR_dB,sum(R_o,2)/realization,'r-o','LineWidth',1.5)