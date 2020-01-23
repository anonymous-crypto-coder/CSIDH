function SIMBA_WITHDUMMY_2()

file = fopen('simba_withdummy_2.h', 'w');

load 'FullParameterSet.mat'

n=length(b);             %Number of primes
B=b;                     %Fixes inconsistent variable naming
batches = length(Batch); %number of SIMBA substrategies across planned rounds

complement_Batch = {};

%Construct the complement batches
for i = 1:batches
  complement_Batch{i} = [];
  for j = 0:n-1
    if max(Batch{i} == j) == 0
      complement_Batch{i} = [complement_Batch{i}, j];
    end
  end
end

sizes = [];
complement_sizes = [];
for i = 1:batches
  sizes(i) = length(Batch{i});
  complement_sizes(i) = length(complement_Batch{i});
end

%Printing some constants
fprintf(file, '#ifndef _SIMBA_PARAMETERS_H_\n')
fprintf(file, '#define _SIMBA_PARAMETERS_H_\n')
fprintf(file, '\n')
fprintf(file, '// SIMBA-(NUMBER_OF_BATCHES)-MY\n')
fprintf(file, '#define NUMBER_OF_BATCHES %2d\n', batches) %Total number of batches (NOT batches per round)
fprintf(file, '#define MY %2d\n', 1) %Always 1
fprintf(file, '\n')
fprintf(file, '// (each entry corresponds to the number of degree-(l_i) to be required in the action: this the one given in Onuki et al. work)\n')

%Bound vector
fprintf(file, 'static int8_t B[] = {\n')
for i = 1:n-1
  fprintf(file, '%3d, ', B(i))
  if mod(i,8) == 0
    fprintf(file, '\n')
  end
end
fprintf(file, '%3d', B(n))
fprintf(file, '\n};\n\n')

%Batch primes
fprintf(file, '// (NUMBER_OF_BATCHES) different subsets (i.e., batches)\n')

%This code requires the batch entries to be listed in the order that the
%isogenies will be computed, rather than the order in which the primes
%will be multiplied in, so we reverse the order of each batch.
for i=1:batches
  fprintf(file, 'static uint8_t BATCH_%1d[] = {\n', i-1)
  for j=length(Batch{i}):-1:2
    fprintf(file, '%3d, ', Batch{i}(j))
    if(mod(j,8) == mod(length(Batch{i})+1,8))
      fprintf(file, '\n')
    end
  end
  fprintf(file, '%3d' ,Batch{i}(1))
  fprintf(file, '\n};\n\n')
end

%Batch sizes
fprintf(file, 'static uint8_t SIZE_OF_EACH_BATCH[NUMBER_OF_BATCHES] = {')
for i = 1:length(sizes)-1
  if(mod(i,8) == 1)
    fprintf(file, '\n');
  end
  fprintf(file, '%3d, ', sizes(i))
end
fprintf(file, '%3d', sizes(end))
fprintf(file, '\n};\n\n')

%Array of batches
fprintf(file, 'static uint8_t *BATCHES[NUMBER_OF_BATCHES] = {')
for i = 0:batches-2
  if(mod(i,4) == 0)
    fprintf(file, '\n')
  end
  fprintf(file, 'BATCH_%1d, ', i)
end
fprintf(file, 'BATCH_%1d', batches-1)
fprintf(file, '\n};\n\n')

%Last isogenies
fprintf(file, 'static uint8_t LAST_ISOGENY[NUMBER_OF_BATCHES] = {')
for i=1:batches-1
  if(mod(i,8) == 1)
    fprintf(file, '\n')
  end
  fprintf(file, '%2d, ', Batch{i}(1))
end
fprintf(file, '%2d', Batch{batches}(1))
fprintf(file, '\n};\n\n')

%Total number of isogenies
fprintf(file, 'static uint16_t NUMBER_OF_ISOGENIES = %3d;\n\n', sum(B))

%Batch complements/sizes
fprintf(file, '// The complement of each batch\n')

fflush(file);

%Batch complement sizes
fprintf(file, 'static uint8_t SIZE_OF_EACH_COMPLEMENT_BATCH[NUMBER_OF_BATCHES] = {\n')
for i=1:batches-1
  fprintf(file, '%2d, ', complement_sizes(i))
  if(mod(i,8) == 0)
    fprintf(file, '\n')
  end
end
fprintf(file, '%2d', complement_sizes(end))
fprintf(file, '};\n\n')

%Batch complements
fprintf(file, 'static uint8_t COMPLEMENT_OF_EACH_BATCH[NUMBER_OF_BATCHES][N] = {\n')
fprintf(file, '{\n') 
for i = 1:length(complement_Batch{1})-1
  fprintf(file, '%2d, ', complement_Batch{1}(i))
  if (mod(i,8) == 0)
    fprintf(file, '\n')
  end
end
fprintf(file, '%2d, \n', complement_Batch{1}(end))
for i = 1:sizes(1)-1
  if (mod(i,8) == 1)
    fprintf(file, '\n')
  end
  fprintf(file, '%2d, ', n)
end
fprintf(file, '%2d', n)
fprintf(file, '}'),

for j = 2:batches
  fprintf(file, ',\n\n{\n') 
  for i = 1:complement_sizes(j)-1
    fprintf(file, '%2d, ', complement_Batch{j}(i))
    if (mod(i,8) == 0)
      fprintf(file, '\n')
    end
  end
  fprintf(file, '%2d, \n', complement_Batch{j}(end))
  for i = 1:sizes(j)-1
    fprintf(file, '%2d, ', n)
    if (mod(i,8) == 0)
      fprintf(file, '\n')
    end
  end
fprintf(file, '%2d', n)
fprintf(file, '\n}'),
end

fflush(file);
fprintf(file, '\n};\n\n#endif\n')

fflush(file);
