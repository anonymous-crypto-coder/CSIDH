function Main()

%{
Computes a bound vector for the parameters and cost model contained in
Parameters.mat. We use our heuristic mathematical program which we solve
approximately by iterative rounding. The result is written to BoundVector.mat.
%}

%Load the output of the first step; in particular, we load in the primes, cost
%model, best permutation and best strategy found.
params = load('Parameters.mat');

%We need to load the variables into a data structure at first, because some of
%of them share names with built-in functions. We then reassign them to make the
%code more concise.

mu = params.mu;
iota = params.iota;
Primes = params.Primes;
n=params.n;
Sigma=params.Sigma;
H=params.H;
V=params.V;

simba_M = length(params.Split);   %The number of SIMBA substrategies.

%kappa is a vector of costs of isogeny kernel point and codomain computation
kappa  = 38*Primes + 10*ceil(log(Primes)/log(2)-1) + 16*ceil(log(Primes)/log(2)) - 36;

Hhat = [H,zeros(n-1,1)];   %Required for defining the objective.
Vhat = [zeros(1,n-1);V];   %Required for defining the objective.

c = zeros(n,1);        %Initialize a vector used to define the objective.
d = zeros(n,1);        %Initialize a vector used to define the objective.

perm = Sigma*(1:n)';

%Define the vectors entry-by-entry
 for i = 1:n
  c(perm(i)) = sum(Hhat(:,i))*mu(perm(i)) + sum(Vhat(i,:))*2*iota(perm(i)) + 1*kappa(perm(i));
  d(perm(i)) = 2*mu(perm(i));
 end

cvx_begin
  variable b(n);         %Bound vector.
  variable lambda(n);    %Measures security entry-by-entry of b.
  variable m;            %Equal to max(b) in an optimal solution.
  variable prr;          %The expected number of extra rounds required.

  minimize( c'*b + d'*(m*ones(n,1) - b) + 2*(simba_M-1)*sum(mu) + 2*sum(mu)*prr )
  subject to
    ones(1,n)*lambda >= 255.997          %Sufficient security
    2.^lambda <= 2*b + ones(n,1)        
    prr*ones(n,1) >= (1./(Primes-1)).*b  
    m*ones(n,1) >= b

cvx_end

%Parameters used to fix the rounded entries of b for iterative rounding.
ub = 100*ones(n,1);    %Upper bound
lb = zeros(n,1);       %Lower boud
nrange = 1:n;          %Needed for slicing
rb = round(b);

%Iterative rounding
while(sum(lb==0)>1)
  unknown_variables = size(lb(lb==0))   %Number of unfixed entires of b
  range = nrange(lb == 0);
  rounding_indices = range(abs(b(range)-rb(range)) == min(abs(b(range)-rb(range))))
  rb = round(b);
  ub(rounding_indices) = rb(rounding_indices); %Fix the rounded entries of
  lb(rounding_indices) = rb(rounding_indices); %b by upper and lower
                                               %bounding them.
  cvx_begin
    variable b(n);         %Bound vector.
    variable lambda(n);    %Measures security entry-by-entry of b.
    variable m;            %Equal to max(b) in an optimal solution.
    variable prr;          %The expected number of extra rounds required.
  
    minimize( c'*b + d'*(m*ones(n,1) - b) + 2*(simba_M-1)*sum(mu) + 2*sum(mu)*prr )
    subject to
      ones(1,n)*lambda >= 255.997
      2.^lambda <= 2*b + ones(n,1)
      prr*ones(n,1) >= (1./(Primes-1)).*b
      b <= ub
      b >= lb
      m*ones(n,1) >= b
  
  cvx_end
end

idx = find(lb-ub);         %Find the last unrounded index
b(idx) = ceil(b(idx));     %Round it up
b = round(b);              %Cast as integer
bhat=b;

%We decrease entries of the bound vector until the protocol is minimally-secure
while sum(log(2*bhat+1)/log(2)) >= 256
  max_idx = max(find(b == max(b)));
  bhat(max_idx) = bhat(max_idx)-1;
  if sum(log(2*bhat+1)/log(2)) >= 256
    b=bhat;
  end
end

save('BoundVector.mat','b')    %Save the final bound vector to a file.
