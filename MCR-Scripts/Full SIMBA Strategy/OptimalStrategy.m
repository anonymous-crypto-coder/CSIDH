function [H,V] = OptimalStrategy(mu, iota)

%{
Returns the optimal full strategy for a given (possibly permuted) cost model.

--------------------------------------------------------------------------------

Input:
------

  mu: Vector of multiplication costs.
iota: Vector of isogeny evaluation costs.

Output:
-------

H: Horizontal strategy matrix for the optimal strategy
V: Vertical strategy matrix for the optimal strategy
%}

n = size(mu)(1);                           %Number of primes

[~, Indices] = CostsAndIndices(mu, iota);  %Optimal substrategy sizes

if n==1;
  %There's only one strategy on one prime
  H = [];
  V = [];
elseif n==2;
  %There's only one strategy on two primes
  H = [1];
  V = [1];
else
  split_index = Indices(n,n);              %The optimal index to divide into
                                           %left and right substrategies

  %Find the optimal left and right substrategoes
  [HL,VL] = OptimalStrategy(mu((n-split_index+1):n), iota((n-split_index+1):n));
  [HR,VR] = OptimalStrategy(mu(1:(n-split_index)), iota(1:(n-split_index)));

  %Join the subsstrategies into one strategy.
  [H,V] = StrategyJoin(HL, VL, HR, VR);
end;
