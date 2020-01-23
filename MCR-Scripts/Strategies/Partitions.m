function P = Partitions(n,k,lb,ub)
%{
Returns the partitions of n into k parts where each part 
is between lb and ub. 

--------------------------------------------------------------------------------

Input:
------

 n: Number to be partitioned.
 k: Number of parts in each partition.
lb: Lower bound on entries of each partition.
ub: Upper bound on entries of each partition.

Output:
-------

P: A matrix whose rows are all partitions of n into k parts between lb and ub.
%}

%This part gets all compositions of n into k parts
m = nchoosek(n+k-1,k-1); 
x = [zeros(m,1), nchoosek((1:(n+k-1))',k-1), (n+k)*ones(m,1)];
P = diff(x,1,2)-1;

%We filter out for compositions that don't satisfy the bounds
P = P((1:length(P))(min(P') >= lb),:);
P = P((1:length(P))(max(P') <= ub),:);

%We sort (to get partitions) and filter out duplicates
P = unique(sort(P,2),'rows');
