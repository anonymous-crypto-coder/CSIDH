function Main()

%{
The main function used to find a strategy and permutation. 
%}

Setup;     %Load the relevant data: n, mu, iota.

%Search for a strategy and permutation
[Sigma, H, V, Cost, Split] = StochasticSearch(mu,iota,flip(eye(n)),100,mu,5); 

save -mat7-binary 'Parameters.mat'      %Save output to 'Parameters.mat'
