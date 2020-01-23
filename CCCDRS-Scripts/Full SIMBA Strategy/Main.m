function Main()

%{
The main function used to find a strategy and permutation. 
%}

Setup;     %Load the relevant data: n, mu, iota.

%Search for a strategy and permutation
[Sigma, H, V, Cost, Split] = SemiExhaustiveSearch(mu,iota,1,1,n,mu);
for i=2:5
  [NewSigma, NewH, NewV, NewCost, NewSplit] = SemiExhaustiveSearch(mu,iota,i,floor(n/(i+2)),ceil(n/i)+15,mu);
  if NewCost < Cost
    Sigma = NewSigma;
    H = NewH;
    V = NewV;
    Cost = NewCost;
    Split = NewSplit;
  end
end

save -mat7-binary 'Parameters.mat'      %Save output to 'Parameters.mat'
