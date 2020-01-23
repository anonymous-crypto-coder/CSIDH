function Main()

%{
Computes a bound vector for the parameters and cost model contained in
Parameters.mat. We use our heuristic mathematical program which we solve
approximately by iterative rounding. The result is written to BoundVector.mat.
%}

load 'Parameters.mat'      %Load the output of the first step; in particular, we
                           %load in the primes, cost model, best permutation and
		           %best strategy found.

simba_M = length(Split);   %The number of SIMBA substrategies.

%kappa is a vector of costs of isogeny kernel point and codomain computation
kappa  = 38*primes + 10*ceil(log(primes)/log(2)-1) + 16*ceil(log(primes)/log(2)) - 36;


Hhat = [H,zeros(n-1,1)];   %Required for defining the objective function.
Vhat = [zeros(1,n-1);V];   %Required for defining the objective function.

c = zeros(n,1);            %Initialize a vector used to define the objective.
d = zeros(n,1);            %Initialize a vector used to define the objective.

%Define the vectors entry-by-entry
 for i = 1:n
  c(perm(i)) = sum(Hhat(:,i))*mu(perm(i)) + sum(Vhat(i,:))*(2*iota(perm(i))) + 1*kappa(perm(i));
  d(perm(i)) = 2*mu(perm(i));
 end

cvx_begina
	variable b(n);         %Bound vector.
	variable lambda(n);    %Measures security entry-by-entry of b.
	variable m;            %Equal to max(b) in an optimal solution.
	variable prr;          %The expected number of extra rounds required.

	minimize( c'*b + d'*(m*ones(n,1) - b) + 2*(simba_M-1)*sum(mu) + 2 * sum(mu) * prr )
	subject to
		ones(1,n)*lambda >= 255.999          %Sufficient security
		2.^lambda <= 2*b + ones(n,1)        
		prr*ones(n,1) >= (1./(primes-1)).*b  
		m*ones(n,1) >= b

cvx_end

%Parameters used to fix the rounded entries of b for iterative rounding.
ub = 100*ones(n,1);    %Upper bound
lb = zeros(n,1);       %Lower boud
nrange = 1:n;
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
		variable b(n);
		variable lambda(n);
		variable m;
		variable prr;
	
		minimize( c'*b + d'*(m*ones(n,1) - b) + 2 *(simba_M-1) * sum(mu) + 2 * sum(mu) * prr )
		subject to
			ones(1,n)*lambda >= 255.999
			2.^lambda <= 2*b + ones(n,1)
			prr*ones(n,1) >= (1./(primes-1)).*b
			b <= ub
			b >= lb
			m*ones(n,1) >= b
	
	cvx_end
end

save('BoundVector.mat',b)  %Save the final bound vector to a file.
