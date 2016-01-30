function [ FRF,FBB ] = MO_AltMin( Fopt, NRF )

[Nt, Ns, K] = size(Fopt);
y = [];
FRF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
while(isempty(y) || abs(y(1)-y(2))>1e-1)
    y = [0,0];
    for k = 1:K      
        FBB(:,:,k) = pinv(FRF) * Fopt(:,:,k); 
        y(1) = y(1) + norm(Fopt(:,:,k) - FRF * FBB(:,:,k),'fro')^2;
    end
    [FRF, y(2)] = sig_manif(Fopt, FRF, FBB);
end

end