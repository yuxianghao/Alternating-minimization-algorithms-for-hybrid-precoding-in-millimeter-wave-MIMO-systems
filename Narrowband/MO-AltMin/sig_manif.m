function [y, cost] = sig_manif(Fopt, FRF, FBB)
[Nt, NRF] = size(FRF);

manifold = complexcirclefactory(Nt*NRF);
problem.M = manifold;

% problem.cost  = @(x) norm( Fopt - reshape(x,Nt,NRF) * FBB,'fro')^2;
% problem.egrad = @(x) -2 * kron(conj(FBB), eye(Nt)) * (Fopt(:) - kron(FBB.', eye(Nt)) * x);
f = Fopt(:);
A = kron(FBB.', eye(Nt));

problem.cost  = @(x) (f-A*x)'*(f-A*x);
problem.egrad = @(x) -2*A'*(f-A*x);

% checkgradient(problem);
warning('off', 'manopt:getHessian:approx');

[x,cost,info,options] = conjugategradient(problem,FRF(:));
% [x,cost,info,options] = trustregions(problem, FRF(:));
% info.iter
y = reshape(x,Nt,NRF);

end