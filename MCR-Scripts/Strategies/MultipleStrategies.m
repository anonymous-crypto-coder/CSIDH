function [perms, Hs, Vs, costs, sizeses] = MultipleStrategies(primes,mu,iota,b)

%{
Finds a complete set of strategies for a given set of primes, cost model, and
bound vector.

--------------------------------------------------------------------------------

Input:
------

primes: Vector of primes.
    mu: Vector of multiplication costs.
  iota: Vector of isogeny-evaluation costs.
     b: Bound vector defining a keyspace.

Output:
-------

  perms: List of permutations found.
     Hs: List of horizontal strategy matrices of the strategies found.
     Vs: List of vertical strategy matrices of the strategies found.
  costs: List of costs of (strategy,permutation) pairs in the cost model.
sizeses: List of SIMBA partitions of the strategies.

%}

n = length(mu);  %Number of primes

nrange = 1:n;    %Needed for slicing later

m = 1;             
M = max(b);      %Total number of strategies required

%Intializing lists
A = {};          %The lists of indices of "necessarily incomplete" primes
l=[];            %The lengths of the A{i}
perms = {};      %The best permutations found for each set of primes
Sigmas = {};     %The permutation matrix representations of the perms{i}
Hs = {};         %The horizontal strategy matrices of the best strategies found
Vs = {};         %The horizontal strategy matrices of the best strategies found
primeses = {};   %The lists of primes indexed by the A{i}
mus = {};        %The sub-vectors of mu indexed by the A{i}
iotas = {};      %The sub-vectors of iota indexed by the A{i}
costs = [];      %The costs of the best strategies found for each step
sizeses = {};    %The best SIMBA partitions found of the l(i)

%Initializing the list entires
for i = m:M
  A{i} = nrange(b >= i); 
  primeses{i} = primes(b >= i);
  mus{i} = mu(b >= i);
  iotas{i} = iota(b >= i);
  l(i) = length(A{i});
  Sigmas{i} = flip(eye(l(i)));
  perms{i} = flip(eye(l(i)))*(1:l(i))';
  costs(i) = 0;
  sizeses{i} = l(i);
end

%The actual search for the strategies and permutations
for i = m:M
  [Sigmas{i}, Hs{i}, Vs{i}, costs(i), sizeses{i}] = StochasticSearch(mus{i}, iotas{i}, Sigmas{i},100,mu,5);
  perms{i} = (Sigmas{i}*(1:l(i))')';
end

EverythingPrint(b, primes, primeses, perms, Vs, sizeses);
