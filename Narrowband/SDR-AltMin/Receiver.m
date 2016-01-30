function [FRF,FBB] = Receiver(Fopt,NRF)

% randomly generate FRF
[Nt,Ns] = size(Fopt);
FRF = [];
for i = 1:NRF
    FRF = blkdiag(FRF, exp(sqrt(-1) * unifrnd (0,2*pi,[Nt/NRF,1])));
end
FRF = 1/sqrt(Nt)*FRF;

y = [];
while(isempty(y) || abs(y(1)-y(2))>1e-3)
    % fix FRF, optimize FBB
    FBB = pinv(FRF) * Fopt;
    
    y(1) = norm(Fopt-FRF*FBB,'fro')^2;
    
    % fix FBB, optimize FRF
    for i = 1:Nt
        m = ceil(i*NRF/Nt);
        FRF(i,m) = 1/sqrt(Nt) * exp( sqrt(-1) * angle( Fopt(i,:)*FBB(m,:)' ) );
    end
    
    y(2) = norm(Fopt-FRF*FBB,'fro')^2;
end

end