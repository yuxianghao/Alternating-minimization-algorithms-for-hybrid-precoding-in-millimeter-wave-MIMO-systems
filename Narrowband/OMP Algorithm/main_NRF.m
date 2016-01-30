clear,clc

% load('Ns=2.mat');
load('Ns=6.mat');
Ns = 6;

NRF = 6:11;

SNR_dB = 0;
SNR = 10.^(SNR_dB./10);
realization = size(H,3);

for r = 1:length(NRF)
    for reali = 1:realization
        [ FRF, FBB ] = OMP( Fopt(:,:,reali), NRF(r), At(:,:,reali) );
        FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
        [ WRF, WBB ] = OMP( Wopt(:,:,reali), NRF(r), Ar(:,:,reali) );
        R(r,reali) = log2(det(eye(Ns) + SNR/Ns * pinv(WRF * WBB) * H(:,:,reali) * FRF * FBB * FBB' * FRF' * H(:,:,reali)' * WRF * WBB));
        R_o(r,reali) = log2(det(eye(Ns) + SNR/Ns * pinv(Wopt(:,:,reali)) * H(:,:,reali) * Fopt(:,:,reali) * Fopt(:,:,reali)' * H(:,:,reali)' * Wopt(:,:,reali)));
    end
end
plot(NRF,sum(R_o,2)/realization,'r-o','LineWidth',1.5)
grid on
hold on
plot(NRF,sum(R,2)/realization,'Marker','^','LineWidth',1.5,'Color',[0 0.498039215803146 0])