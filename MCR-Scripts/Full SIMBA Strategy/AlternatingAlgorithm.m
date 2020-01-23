function [Sigma, H, V, cost] = AlternatingAlgorithm(mu, iota, Sigma0, sizes, mu_total)

%{
A naive alternating algorithm to optimize the startegy and permutation for given
multiplication and isogeny evaluation costs.

--------------------------------------------------------------------------------

Inputs:
-------
      mu: Vector of multiplication costs for the set of "active" primes.
    iota: Vector of isogeny-evaluation costs for the set of "active" primes.
  Sigma0: Initial choice of permutation matrix. Typically Sigma0 = flip(eye(n)).
   sizes: SIMBA partition
mu_total: Vector of multiplcation costs for the set of all primes.

Outputs:
--------

Sigma: The best permutation found.
    H: The horizontal strategy matrix of the best strategy found.
    V: The vertical strategy matrix of the best strategy found.
 Cost: The cost of the best strategy and permutation found.
%}


n = length(mu);                        %Number of primes
m = length(sizes);                     %Size of SIMBA partition


T1 = [eye(n-1), zeros(n-1,1)];         %Some matrices required to define
T2 = [zeros(n-1,1), eye(n-1)];         %the objective function

indices = [0, cumsum(sizes)'];         %Indices of the the start points of the
                                       %SIMBA substrategies


cost = n*sum(mu) + n*sum(iota);        %An upper bound on the cost of any
                                       %strategy and permutation

old_cost = cost + 1;                   %Something bigger than the current cost
Sigma = Sigma0;                        %Initial permutation

while(cost < old_cost)                 %We search until our cost stops improving
  old_cost=cost;                       %Store previous cost

  for i=1:m
    %We find the optimal single strategy for each subset of primes
    muhat = Sigma*mu;                  %Subset of "active" multiplication costs
    iotahat = Sigma*iota;              %Subset of "active" isogeny costs
    [Hs{i}, Vs{i}] = OptimalStrategy(muhat(indices(i)+1:indices(i+1)), iotahat(indices(i)+1:indices(i+1)));
  end

  H = Hs{1};                           %Initial value
  V = Vs{1};                           %Initial value

  for i=2:m
    %We need to SIMBA-join the strategy matrices
    [H,V] = SIMBAJoin(Hs{i}, Vs{i}, H, V);
  end
   
  Sigma = OptimalPermutationMatrix(mu, iota, H, V); 
  cost = vec(T2'*V*ones(n-1,1)*iota' + T1'*H'*ones(n-1,1)*mu')'*vec(Sigma); 

end

cost = cost+m*sum(mu_total)-sum(mu);   %Accounts for additional multiplications
				       %required for SIMBA (to eliminite
				       %"inactive" primes)

