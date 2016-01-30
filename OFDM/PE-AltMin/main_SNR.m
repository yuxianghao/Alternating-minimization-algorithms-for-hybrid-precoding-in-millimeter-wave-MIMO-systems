% In this code, 1000 realizations, evaluaation of AE-AltMin and OMP
clear,clc

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

R_o = zeros(smax,realization);
R = zeros(smax,realization);
RO = zeros(smax,realization);
RM = zeros(smax,realization);

for reali = 1:realization
    reali
    
    AoD = zeros(2,Nray);
    AoA = zeros(2,Nray);
    Ht = zeros(Nr,Nt,Nc);
    H = zeros(Nr,Nt,K);
    At = zeros(Nt,Nc*Nray);
    Ar = zeros(Nr,Nc*Nray);
    Fopt = zeros(Nt,Ns,K);
    Wopt = zeros(Nr,Ns,K);
    
    for c = 1:Nc
        AoD_m = unifrnd(0,2*pi,1,2);
        AoA_m = unifrnd(0,2*pi,1,2);
        
        AoD(1,:) = laprnd(1,Nray,AoD_m(1),angle_sigma);
        AoD(2,:) = laprnd(1,Nray,AoD_m(2),angle_sigma);
        AoA(1,:) = laprnd(1,Nray,AoA_m(1),angle_sigma);
        AoA(2,:) = laprnd(1,Nray,AoA_m(2),angle_sigma);
        
%         Ht(:,:,c) = zeros(Nr,Nt);
        for j = 1:Nray
            temp = (c-1)*Nray+j;
            At(:,temp) = array_response(AoD(1,j),AoD(2,j),Nt);
            Ar(:,temp) = array_response(AoA(1,j),AoA(2,j),Nr);
            alpha = normrnd(0,sqrt(sigma/2)) + 1i*normrnd(0,sqrt(sigma/2));
            Ht(:,:,c) = Ht(:,:,c) + alpha * Ar(:,temp) * At(:,temp)';
        end
    end
    
    for k = 1:K
%         H(:,:,k) = zeros(Nr,Nt);
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
        %% OMP Algorithm
    [ FRFO, FBBO ] = OMP( Fopt, NRF, At );
    for k = 1:K
        FBBO{k} = sqrt(Ns) * FBBO{k} / norm(FRFO * FBBO{k},'fro');
    end
    [ WRFO, WBBO ] = OMP( Wopt, NRF, Ar );
        %% PE-AltMin
    [ FRF, FBB ] = PE_AltMin( Fopt, NRF );
    for k = 1:K
        FBB(:,:,k) = sqrt(Ns) * FBB(:,:,k) / norm(FRF * FBB(:,:,k),'fro');
    end
    [ WRF, WBB ] = PE_AltMin( Wopt, NRF );
        %% Calculate the spectral efficiency
    for k = 1:K
        for s = 1:smax
            R_o(s,reali) = R_o(s,reali) + log2(det(eye(Ns) + SNR(s)/Ns * pinv(Wopt(:,:,k)) * H(:,:,k) * Fopt(:,:,k) * Fopt(:,:,k)' * H(:,:,k)' * Wopt(:,:,k)))/K;
            R(s,reali) = R(s,reali) + log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRF * WBB(:,:,k)) * H(:,:,k) * FRF * FBB(:,:,k) * FBB(:,:,k)' * FRF' * H(:,:,k)' * WRF * WBB(:,:,k)))/K;
            RO(s,reali) = RO(s,reali) + log2(det(eye(Ns) + SNR(s)/Ns * pinv(WRFO * WBBO{k}) * H(:,:,k) * FRFO * FBBO{k} * FBBO{k}' * FRFO' * H(:,:,k)' * WRFO * WBBO{k}))/K;
        end
    end
end

plot(SNR_dB,sum(R_o,2)/realization,'r-o','LineWidth',1.5)
grid on
hold on
plot(SNR_dB,sum(R,2)/realization,'Marker','>','LineWidth',1.5,'Color',[0 0.447058826684952 0.74117648601532])
plot(SNR_dB,sum(RO,2)/realization,'Marker','^','LineWidth',1.5,'Color',[0 0.498039215803146 0])
xlabel('SNR (dB)')
ylabel('Spectral Efficiency (bits/s/Hz)')
legend('Optimal Digital Precoder','AE-AltMin','OMP Algorithm','Location','NW')