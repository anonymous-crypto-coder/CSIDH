function Main()

%{
The main function used to find strategies and permutations. 
%}

Setup;     %Load the relevant data: n, mu, iota.
load 'BoundVector.mat'

%Search for a set of strategies and permutations
[perms, Hs, Vs, Costs, Splits] = MultipleStrategies(primes,mu,iota,b);

%Save the output
save -mat-binary 'FullParameterSet.mat'
