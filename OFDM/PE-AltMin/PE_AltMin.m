function [ FRF,FBB ] = PE_AltMin( Fopt, NRF )

[Nt, Ns, K] = size(Fopt);
mynorm = [];
FRF = exp( 1i * unifrnd (0,2*pi,Nt,NRF) );
% FBB = zeros(NRF, Ns, K);
while (isempty(mynorm) || abs( mynorm(1) - mynorm(2) ) > 1e-1)
    mynorm = [0,0];
    temp = zeros(Nt, NRF);
    for k = 1:K
        [U,S,V] = svd(Fopt(:,:,k)'*FRF);
        FBB(:,:,k) = V(:,[1:Ns])*U';
        mynorm(1) = mynorm(1) + norm(Fopt(:,:,k) * FBB(:,:,k)' - FRF,'fro')^2;
        temp = temp + Fopt(:,:,k) * FBB(:,:,k)';
    end

    FRF = exp(1i * angle(temp));
    for k = 1:K
        mynorm(2) = mynorm(2) + norm(Fopt(:,:,k) * FBB(:,:,k)' - FRF,'fro')^2;
    end
end
end