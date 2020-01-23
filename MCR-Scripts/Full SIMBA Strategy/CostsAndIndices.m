function [Costs, Indices] = CostsAndIndices(mu, iota)
%{
Computes the matrix of costs of full strategies for all contiguous and
compatible subsets of given vectors of multiplication and isogeny evaluation
costs.

--------------------------------------------------------------------------------

Input:
------

  mu: Vector of multiplication costs.
iota: Vector of isogeny evaluation costs.

Output:
-------

  Costs: Costs(j,k) is the optimal cost of a strategy on primes k-j+1, ..., k.
Indices: Indices(j,k) is the smallest integer between k-j+k and k which, if we
         split the list at this index, there are corresponding left and right
	 substrategies whose join achieves cost Costs(j,k).
%}

n = size(mu)(1);

Costs = zeros(n,n);
%Costs(j,k) will be the optimal cost of a strategy on primes k-j+1, ..., k

Indices = zeros(n,n);
%Indices(j,k) will tell us where to split the list of primes k-j+1, ..., k to achieve Costs(j,k).

for j=2:n
  for k=j:n
    cost_list = zeros(j-1,1);                          %Contains the costs
    for l=1:(j-1)
      %The primes in play here are k-j+1, ... k. For mu we use the first j-l. For iota we use the last l.
      cost_list(l) = Costs(l,k) + Costs(j-l,k-l) + sum(mu((k-j+1):(k-l))) + sum(iota((k-l+1):(k)));
    end
    Costs(j,k) = min(cost_list);                       %The optimal cost
    Indices(j,k) = find(cost_list == Costs(j,k))(1);   %The corresponding index
  end
end
