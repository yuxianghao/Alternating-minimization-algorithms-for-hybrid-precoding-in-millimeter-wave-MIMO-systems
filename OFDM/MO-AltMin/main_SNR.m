% This code is specialized for MO-AltMin, using parallel in MO-AltMin.m to speed up
clear,clc

addpath(pwd);
cd manopt;
addpath(genpath(pwd));
cd ..;

Ns = 4; % # of streams
NRF = 5;

Nc = 5; % # of clusters
Nray = 10; % # of rays in each cluster

Nt = 144; % # of transmit antennas
Nr = 36; % # of receive antennas

angle_sigma = 10/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx

gamma = sqrt((Nt*Nr)/(Nc*Nray)); %normalization factor
sigma = 1; %according to the normalization condition of the H

realization = 1000;
count = 0;
K = 128;

SNR_dB = -15:5:10;
SNR = 10.^(SNR_dB./10);
smax = length(SNR);

RM = zeros(smax,realization);

for reali = 1:realization
    reali
    
    for c = 1:Nc
        AoD_m = unifrnd(0,2*pi,1,2);
        AoA_m = unifrnd(0,2*pi,1,2);
        
        AoD(1,:) = laprnd(1,Nray,AoD_m(1),angle_sigma);
        AoD(2,:) = laprnd(1,Nray,AoD_m(2),angle_sigma);
        AoA(1,:) = laprnd(1,Nray,AoA_m(1),angle_sigma);
        AoA(2,:) = laprnd(1,Nray,AoA_m(2),angle_sigma);
        
        Ht(:,:,c) = zeros(Nr,Nt);
        for j = 1:Nray
            temp = (c-1)*Nray+j;
            At(:,temp) = array_response(AoD(1,j),AoD(2,j),Nt);
            Ar(:,temp) = array_response(AoA(1,j),AoA(2,j),Nr);
            alpha = normrnd(0,sqrt(sigma/2)) + 1i*normrnd(0,sqrt(sigma/2));
            Ht(:,:,c) = Ht(:,:,c) + alpha * Ar(:,temp) * At(:,temp)';
        end
    end
    
    for k = 1:K
        H(:,:,k) = zeros(Nr,Nt);
        for c = 1:Nc
            H(:,:,k) = H(:,:,k) + Ht(:,:,c) * exp(-1i*2*pi/K*(k-1)*(c-1));
        end
        H(:,:,k) = H(:,:,k) * gamma;
        if(rank(H(:,:,k))>=Ns)
            count = count + 1;
            
            [U,S,V] = svd(H(:,:,k));
            Fopt(:,:,k) = V([1:Nt],[1:Ns]);
            Wopt(:,:,k) = U([1:Nr],[1:Ns]);
        end
    end
        %% MO-AltMin
    [ FRFM, FBBM ] = MO_AltMin( Fopt, NRF );
    parfor k = 1:K
        FBBM(:,:,k) = sqrt(Ns) * FBBM(:,:,k) / norm(FRFM * FBBM(:,:,k),'fro');
    end
    [ WRFM, WBBM ] = MO_AltMin( Wopt, NRF );

        %% Calculate the spectral efficiency
    for k = 1:K
        for s = 1:smax
            RM(s,reali) = RM(s,reali) + log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRFM * WBBM(:,:,k)) * H(:,:,k) * FRFM * FBBM(:,:,k) * FBBM(:,:,k)' * FRFM' * H(:,:,k)' * WRFM * WBBM(:,:,k)))/K;
        end
    end
end

plot(SNR_dB,sum(RM,2)/realization,'k-p','LineWidth',1.5)
grid on
hold on
xlabel('SNR (dB)')
ylabel('Spectral Efficiency (bits/s/Hz)')
legend('Optimal Digital Precoder','AE-AltMin','OMP Algorithm','Location','NW')