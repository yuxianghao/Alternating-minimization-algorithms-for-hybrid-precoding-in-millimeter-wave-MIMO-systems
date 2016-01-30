function [ FRF,FBB ] = PE_AltMin( Fopt, NRF )

[Nt, Ns] = size(Fopt);
mynorm = [];
FRF = exp( sqrt(-1) * unifrnd (0,2*pi,Nt,NRF) );
while (isempty(mynorm) || abs( mynorm(end) - mynorm(end-1) ) > 1e-3)
    [U,S,V] = svd(Fopt'*FRF);
    FBB = V(:,[1:Ns])*U';
    mynorm = [mynorm, norm(Fopt * FBB' - FRF,'fro')^2];
    FRF = exp(1i * angle(Fopt * FBB'));
    mynorm = [mynorm, norm(Fopt * FBB' - FRF,'fro')^2];
end
end