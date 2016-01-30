function [ FRF, FBB ] = OMP( Fopt, NRF, At )
K = size(Fopt,3);
FRF = [];
Fres = Fopt;
for i = 1:NRF
    temp = 0;
    for k = 1:K
        PU(:,:,k) = At' * Fres(:,:,k);
        temp = temp + sum( abs(PU(:,:,k)).^2, 2 );
    end
    [aa,bb] = max(temp);
    FRF = [FRF , At(:,bb)];
    for k = 1:K
        FBB{k} = pinv(FRF) * Fopt(:,:,k);
        Fres(:,:,k) = (Fopt(:,:,k) - FRF * FBB{k}) / norm(Fopt(:,:,k) - FRF * FBB{k},'fro');
    end
end

end