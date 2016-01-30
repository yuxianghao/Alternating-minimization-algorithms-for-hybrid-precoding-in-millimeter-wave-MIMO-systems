clear,clc

data = load('channel.mat');

Ns = 1;
SNR_dB = -35:5:5;
SNR = 10.^(SNR_dB./10);
realization = size(data.alpha,2);
smax = length(SNR);% enable the parallel

for reali = 1:realization
    [cc,num] = max(abs(data.alpha(:,reali)));
    FRF = data.At(:,num,reali);
    WRF = data.Ar(:,num,reali);
    for s = 1:smax
        R(s,reali) = log2(det( eye(Ns) + SNR(s)/Ns * pinv(WRF) * data.H(:,:,reali) * FRF * FRF' * data.H(:,:,reali)' * WRF) );
    end
end

plot(SNR_dB,sum(R,2)/realization,'b-*','LineWidth',1.5)
grid on
hold on