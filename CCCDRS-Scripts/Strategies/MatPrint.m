function MatPrint(V,sizes,file,j)

%{
Prints a vertical strategy matrix to a file, in a way that is useful to us.

--------------------------------------------------------------------------------

Input:
------

    V: A vertical strategy matrix.
sizes: The sizes of the SIMBA substrategies in the underlying strategy.
 file: The output file.
    j: An integer indexing the underlying strategy (useful when many strategies
       are relevant).

Output:
-------

None.
%}

indices = [0, cumsum(sizes)']; %The indices where the SIMBA substrategies begin

%Honestly this is a mess. Trust me that it prints the matrix correctly.
for i = 1:length(indices)-1
  fprintf(file, '\nV%1d%1d = [\n', j, i)
  Vij = V(indices(i)+1:indices(i+1)-1,indices(i)+1:indices(i+1)-1);
  for k = 1:size(Vij)(1)-1
    fprintf(file, '[');
    for l = 1:size(Vij)(2)-1
      fprintf(file, '%1d, ',Vij(k,l))
    end
    fprintf(file, '%1d],\n', V(k,end))
  end
  fprintf(file, '[');
  for l = 1:size(Vij)(2)-1
    fprintf(file, '%1d, ',Vij(end,l))
  end
  fprintf(file, '%1d]\n];', Vij(end,end))
end

fflush(file);  %Ensures that all the required text is printed
