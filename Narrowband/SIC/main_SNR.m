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

for reali = 1:realization
    [ FRF, FBB ] = SDR_AltMin( Fopt(:,:,reali), NRF);
    for s = 1:smax
        R(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(Wopt(:,:,reali)) * H(:,:,reali) * FRF * FBB * FBB' * FRF' * H(:,:,reali)' * Wopt(:,:,reali)));
        [ FRFS, FBBS ] = SIC( Fopt(:,:,reali), H(:,:,reali), NRF, SNR(s) );
        RS(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(Wopt(:,:,reali)) * H(:,:,reali) * FRFS * FBBS * FBBS' * FRFS' * H(:,:,reali)' * Wopt(:,:,reali)));
        Ro(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(Wopt(:,:,reali)) * H(:,:,reali) * Fopt(:,:,reali) * Fopt(:,:,reali)' * H(:,:,reali)' * Wopt(:,:,reali)));
    end
end
plot(SNR_dB,sum(R,2)/realization,'Marker','diamond','LineWidth',1.5,'Color',[0.87058824300766 0.490196079015732 0]);
grid on
hold on
plot(SNR_dB,sum(RS,2)/realization,'c-+','LineWidth',1.5);
plot(SNR_dB,sum(Ro,2)/realization,'r-o','LineWidth',1.5);