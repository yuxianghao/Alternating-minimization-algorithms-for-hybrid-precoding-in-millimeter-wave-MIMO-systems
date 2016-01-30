function [y, cost] = sig_manif(Fopt, FRF, FBB)
[Nt, NRF] = size(FRF);
K = size(FBB,3);

manifold = complexcirclefactory(Nt*NRF);
problem.M = manifold;

parfor k = 1:K
    temp = Fopt(:,:,k);
    A = kron(FBB(:,:,k).', eye(Nt));
    C1(:,:,k) = temp(:)'*A;
    C2(:,k) = A'*temp(:);
    C3(:,:,k) = A'*A;
    C4(k) = norm(temp,'fro')^2;
end
B1 = sum(C1,3);
B2 = sum(C2,2);
B3 = sum(C3,3);
B4 = sum(C4);

problem.cost = @(x) -B1*x - x'*B2 + trace(B3*x*x') + B4;
problem.egrad = @(x) -2*B2 + 2*B3*x;

% checkgradient(problem);
warning('off', 'manopt:getHessian:approx');

[x,cost,info,options] = conjugategradient(problem,FRF(:));
% [x,cost,info,options] = trustregions(problem, FRF(:));
y = reshape(x,Nt,NRF);

end


% problem.cost = @(x) mycost(Fopt, FBB, x);
% function g = mycost(Fopt, FBB, x)
%     g = 0;
%     for k = 1:K
%         temp = Fopt(:,:,k);
%         A = kron(FBB(:,:,k).', eye(Nt));
%         g = g + ( temp(:) -  A*x )' * ( temp(:) - A*x );
%     end
% end
%
% problem.egrad = @(x) mygrad(Fopt, FBB, x);
% function g = mygrad(Fopt, FBB, x)
%     g = 0;
%     for k = 1:K
%         temp = Fopt(:,:,k);
%         A = kron(FBB(:,:,k).', eye(Nt));
%         g = g -2*A'*temp(:) + 2*A'*A*x;
%     end
% end

% problem.costgrad = @(x) mycostgrad(Fopt, FBB, x);
% function [h, g] = mycostgrad(Fopt, FBB, x)
%     h = 0;
%     g = 0;
%     for k = 1:K
%         temp = Fopt(:,:,k);
%         A = kron(FBB(:,:,k).', eye(Nt));
%         g = g -2*A'*temp(:) + 2*A'*A*x;
%         h = h + ( temp(:) -  A*x )' * ( temp(:) - A*x );
%     end
% end