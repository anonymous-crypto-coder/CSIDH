function [BestSigma, BestH, BestV, BestCost, BestSplit] = StochasticSearch(mu,iota,Sigma0,trials,mu_total,q)

%{
Searches for a SIMBA partition, strategy, and permutation which permforms well
for the given parameters.

--------------------------------------------------------------------------------

Inputs:
-------

      mu: Vector of multiplication costs for the set of "active" primes.
    iota: Vector of isogeny-evaluation costs for the set of "active" primes.
  Sigma0: Initial choice of permutation matrix. Typically Sigma0 = flip(eye(n)).
  trials: Number of iterations of the alternating algorithm.
mu_total: Vector of multiplcation costs for the set of all primes.
       q: An upper bound on the maximum partition size.

Outputs:
--------

BestSigma: The best permutation found.
    BestH: The horizontal strategy matrix of the best strategy found.
    BestV: The vertical strategy matrix of the best strategy found.
 BestCost: The cost of the best strategy and permutation found.
BestSplit: The SIMBA partition of the best SIMBA strategy found.
%}

n = length(mu);

%Start with an equitable SIMBA partition into 3 parts (this choice is
%essentially arbitrary; we could have the initial SIMBA partition be input to
%the function instead, but it doesn't matter).
BestSplit = floor(n/3)*[1,1,1]';
i = 1;
while sum(BestSplit) < n
  BestSplit(i) = BestSplit(i) + 1;
  i = i+1;
end

%Apply the alternating algorithm starting with the permutation Sigma0 using the
%SIMBA partition BestSplit. Save the result.
[BestSigma, BestH, BestV, BestCost] = AlternatingAlgorithm(mu,iota,Sigma0,BestSplit,mu_total);

%We find and store all possible partitions
P = {};
P{1} = [n];
for i=2:q
  P{i} = Partitions(n,i,floor(n/(i+2)),ceil(n/i) + 15);
  %These upper and lower bounds on the entries of the partitions are somewhat
  %arbitrary; however, it dramatically reduces the number of partitions we need
  %to search through.
end

%Now we try to improve on the result
for i=1:trials
  idx1 = randi(q);                 %Pick a random partition size
  idx2 = randi(length(P{idx1}));   %Pick a random partition
  s = P{idx1}(idx2,:)';

  %Search for a new strategy and permutation for the new SIMBA partition.
  [Sigma, H, V, cost] = AlternatingAlgorithm(mu,iota,BestSigma,s,mu_total);

  if(cost < BestCost)
    %We have found a better SIMBA partition, strategy, and permutation; save it.
    BestSigma = Sigma;
    BestH = H;
    BestV = V;
    BestCost = cost;
    BestSplit = s;
  end
end
