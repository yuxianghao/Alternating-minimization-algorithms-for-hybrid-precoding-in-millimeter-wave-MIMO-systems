function [ FRF,FBB ] = MO_AltMin( Fopt, NRF )

[Nt, Ns] = size(Fopt);
y = [];
FRF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
while(isempty(y) || abs(y(1)-y(2))>1e-3)
    FBB = pinv(FRF) * Fopt;
    y(1) = norm(Fopt - FRF * FBB,'fro')^2;
    [FRF, y(2)] = sig_manif(Fopt, FRF, FBB);
end

end