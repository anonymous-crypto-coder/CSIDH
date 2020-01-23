function Sigma = OptimalPermutationMatrix(mu,iota,H,V)

%{
Finds the optimal permutation for a given cost model and strategy when the
two-point method is being used.

--------------------------------------------------------------------------------

Input:
------

  mu: Vector of multiplication costs.
iota: Vector of isogeny-evaluation costs.
   H: Horizontal strategy matrix.
   V: Vertical strategy matrix.

Output:
-------

Sigma: Permutation matrix encoding the optimal permutation for the given
       strategy and cost model.
%}

%Suppress LP solver output messages
%options = optimset('linprog');
%options.Display = 'off';

n=length(mu);                                              %Number of primes

T1 = [eye(n-1), zeros(n-1,1)];                             %Used to define the
T2 = [zeros(n-1,1), eye(n-1)];                             %coefficient vector

c=T2'*V*ones(n-1,1)*2*iota' + T1'*H'*ones(n-1,1)*mu';       %Coefficient matrix
c=c(:);	                                                   %Vectorize

A = zeros(2*n,n*n);    %Will hold the equality constraint matrix

for i=1:n
  A(i,1+n*(i-1):1:n*i)=1;
  A(n+i,i:n:n*(n-1)+i)=1;
end

b = ones(2*n,1);       %Equality constraint righthand side vector

%xopt = linprog(c,[],[],A,b,zeros(size(c)),ones(size(c)),options); %Solve linear program
xopt = glpk(c,A,b);

Sigma = reshape(xopt, [n,n]);      %Reshape the solution to be square
