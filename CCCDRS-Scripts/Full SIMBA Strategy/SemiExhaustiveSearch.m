function [BestSigma, BestH, BestV, BestCost, BestSizes] = SemiExhaustiveSearch(mu,iota,m,lb,ub,mutot)

%{
Performs a (somewhat) exhaustive search for an optimal SIMBA strategy for the
two-point method. It is somewhat exhaustive in the sense that it is exhaustive
on the set of SIMBA strategies which have only multiplication-based SIMBA
substrategies of size between lb and ub (inclusive).

--------------------------------------------------------------------------------

Input:
------

   mu: Vector of multiplication costs for "active" primes.
 iota: Vector of isogeny evaluation costs for  "active" primes.
    m: Number of SIMBA substrategies in each strategy.
   lb: Lower bound on the size of a SIMBA substrategy.
   ub: Upper bound on the size of a SIMBA subststrategy.
mutot: Vector of multiplication costs for all primes.

Output:
-------

BestSigma: Matrix representation of best permutation found.
    BestH: Horiztontal strategy matrix of best strategy found.
    BestV: Vertical strategy matrix of best strategy found.
 BestCost: Cost of best strategy and permutation found.
BestSizes: SIMBA partition of the best strategy found.

%}

n = length(mu);                %Number of active primes

P = Partitions(n,m,lb,ub);     %Relevant partitions of n

T1 = [eye(n-1), zeros(n-1,1)]; %Required to build the
T2 = [zeros(n-1,1), eye(n-1)]; %coefficient vector

%We initialize with a naive multiplication-based strategy, just to have a
%starting point.

BestH = tril(ones(n-1,n-1));
BestV = [ones(n-1,1) zeros(n-1,n-2)];
BestSigma = OptimalPermutationMatrix(mu,iota,BestH,BestV);
BestCost = vec(T2'*BestV*ones(n-1,1)*2*iota' + T1'*BestH'*ones(n-1,1)*mu')'*vec(BestSigma) + 2*sum(mutot) - 2*sum(mu);
BestSizes = n;

for i=1:size(P)(1)     %We go through every admissible partition.
  sizes = P(i,:);      %The ith row of P is our SIMBA partition.

  %We construct the strategy matrices.
  H = tril(ones(sizes(1)-1, sizes(1)-1));
  V = [ones(sizes(1)-1,1) zeros(sizes(1)-1,sizes(1)-2)];
  
  for j=2:m
    [H,V] = SIMBAJoin( tril(ones(sizes(j)-1, sizes(j)-1)), [ones(sizes(j)-1,1) zeros(sizes(j)-1,sizes(j)-2)], H, V);
  end
  
  %Compute the optimal permutation matrix across all SIMBA substrategies.
  Sigma = OptimalPermutationMatrix(mu,iota,H,V);   

  %Compute the cost of this permutation and strategy.
  cost = vec(T2'*V*ones(n-1,1)*2*iota' + T1'*H'*ones(n-1,1)*mu')'*vec(Sigma) + 2*m*sum(mutot) - 2*sum(mu);

  if cost < BestCost
    %We have found a better permutation and strategy.
    BestSigma = Sigma;
    BestH = H;
    BestV = V;
    BestCost = cost;
    BestSizes = sizes;
  end
end
