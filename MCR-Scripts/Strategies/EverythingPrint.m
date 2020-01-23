function EverythingPrint(b, primes, primeses, perms, Vs, sizeses)

%{
Prints the vertical strategy matrices and permutations in a way that is
convenient for us, and into the files we like to use.

-------------------------------------------------------------------------------

Input:
------

  primes: The complete list of primes for the CSIDH parameter set.
primeses: The lists of primes used in the smaller strategies.
   perms: The permutations of the primeses{i}.
      Vs: The vertical strategy matrices.
 sizeses: The SIMBA partitions of the strategies.

Output:
-------

None.

%}

m = 1;
M = max(b);

file = fopen('Vs.txt','a');    %Open the matrix output file

%Print all the matrices
for i = m:M
  MatPrint(Vs{i},sizeses{i},file,i);
end

file = fopen('Perms.txt','a'); %Open the permutation output file

%Print all the permutations
for i = m:M
  PermPrint(primes, primeses{i}, perms{i}, sizeses{i}', file, i);
end
