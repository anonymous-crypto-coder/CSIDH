function permutedprimes = PermPrint(fullprimes, subsetprimes, subsetperm, sizes, file, j)

%{
Prints a list of permutations in a way that is convenient for postprocessing.

--------------------------------------------------------------------------------

Input:
------

  fullprimes: The complete list of primes for the CSIDH parameters.
subsetprimes: The subset of primes which subsetperm indexes
  subsetperm: A permutation of the {1, 2, ..., length(fullprimes)}
       sizes: The sizes of the SIMBA substrategies in the underlying strategy.
        file: The output file
           j: An integer indexing the underlying permutation (useful when many
	      permutations are relevant).

Output:
-------

permutedprimes: 

%}


permutedprimes = zeros(length(subsetprimes),1);
for i = 1:length(subsetprimes)
  permutedprimes(i) = find(fullprimes == subsetprimes(subsetperm)(i));
end

indices = [0; cumsum(sizes)'];

%Honestly this is a mess. Trust me that it prints the permutation correctly.
for i = 1:length(indices)-1
  fprintf(file, 'perm%1d%1d = [', j, i)
  for k = indices(i)+1 : indices(i+1)-1
    fprintf(file, '%2d,', permutedprimes(k)-1)
  end
  fprintf(file, '%2d];\n', permutedprimes(indices(i+1))-1)
end

fflush(file);  %Ensures that all the required text is printed
