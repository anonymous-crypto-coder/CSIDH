function Batch = MultipleStrategies(primes,mu,iota,b)

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

m = min(b);             
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
  S={};
  H={};
  V={};
  C=[];
  Sz={};
  for j = 1:4
    [S{j},H{j},V{j},C(j),Sz{j}] = SemiExhaustiveSearch(mus{i},iotas{i},j,floor(l(i)/(j+2)), ceil(l(i)/j) + 15, mu);
  end
  [~,jstar] = min(C);
  Sigmas{i} = S{jstar};
  Hs{i} = H{jstar};
  Vs{i} = V{jstar};
  costs(i) = C(jstar);
  sizeses{i}= Sz{jstar}';
  perms{i} = (Sigmas{i}*(1:l(i))')';
end

%We need a complete collection of SIMBA strategies but we don't want to do
%irrelevant work, so if there is more than one round where all primes are
%active, we search once for the last round where all the primes are active (in
%the above loop) and use what we find in all earlier rounds.
if m != 1
  for i = 1:m-1
    Sigmas{i} = Sigmas{m};
    Hs{i} = Hs{m};
    Vs{i} = Vs{m};
    costs(i) = costs(m);
    sizeses{i}= sizeses{m};
    perms{i} = perms{m};
    primeses{i} = primeses{m};
  end
end

%Saving everything in the format we use for making header files.
%Honestly this is a mess. Trust me that it saves the batches correctly.
Batch = {};
k=0;
for i = 1:M
  permutedprimes = zeros(length(primeses{i}),1);
  for j = 1:length(primeses{i})
    permutedprimes(j) = find(primes == primeses{i}(perms{i})(j));
  end
  indices = [0; cumsum(sizeses{i})];
  for j=1:length(indices)-1
    k=k+1;
    Batch{k} = (permutedprimes(indices(j)+1 : indices(j+1))-1)';
  end
end
