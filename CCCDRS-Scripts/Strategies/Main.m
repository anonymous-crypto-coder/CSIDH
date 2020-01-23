function Main()

%{
The main function used to find a strategy and permutation. 
%}

load 'Parameters.mat'
load 'BoundVector.mat'

%Search for a strategy and permutation
Batch = MultipleStrategies(Primes,mu,iota,b);

save -mat7-binary 'FullParameterSet.mat' b Batch %Save output to 'Parameters.mat'
