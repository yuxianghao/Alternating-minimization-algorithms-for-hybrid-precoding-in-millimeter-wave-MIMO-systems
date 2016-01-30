function [ A, D ] = SIC( Fopt, H, N , SNR)
[Nt,Ns] = size(Fopt);
M = Nt/N;
T = eye(N);
A = [];
D = [];
P = [];
R = [eye(M),zeros(M,M*(N-1))];
G = R*H'*H*R';
for n = 1:N
    [U,S,V] = svd(G);
    v1 = V(:,1);
    
    A = blkdiag( A, exp(1i*angle(v1))/sqrt(M) );
    D = blkdiag( D, norm(v1,1)/sqrt(M) );
    
%     P = [P, [zeros(M*(n-1),1);norm(v1,1)/M*exp(1i*angle(v1));zeros(M*(N-n),1)] ];    
%     T = eye(size(H,1)) + SNR/Ns * H*P*P'*H';
%     R = [zeros(M,M*(n-1)),eye(M),zeros(M,M*(N-n))];
%     G = R*H'*inv(T)*H*R';

    G = G - SNR/Ns * S(1)^2 * v1*v1' / (1+SNR/Ns*S(1));
end

