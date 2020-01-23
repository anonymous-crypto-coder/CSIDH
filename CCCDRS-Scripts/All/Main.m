function Main()

%{
The main function used to find a strategy and permutation. 
%}

%Load CSIDH-512 parameters and cost model

n=74;      %Number of primes

Primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587]';     %The primes

mu = (3/4)*(112*ceil(log(Primes)/log(2))-56);      %Multiplication costs
iota = 20*Primes-4;                                %Isogeny evaluation costs

fprintf("Searching for Initial Strategy and Permutation...")

[Sigma, H, V, Cost, Split] = SemiExhaustiveSearch(mu,iota,1,1,n,mu);
for i=2:5
  [NewSigma, NewH, NewV, NewCost, NewSplit] = SemiExhaustiveSearch(mu,iota,i,max(2,floor(n/(i+2))),ceil(n/i)+15,mu);
  if NewCost < Cost
    Sigma = NewSigma;
    H = NewH;
    V = NewV;
    Cost = NewCost;
    Split = NewSplit;
  end
end

fprintf(" Found.\n")
%disp("Initial Strategy and Permutation Found.")
fprintf("Searching for Bound Vector...")

simba_M = length(Split);   %The number of SIMBA substrategies.

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

cvx_begin quiet
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

fprintf(" Fractional solution found.\n")
%disp("Optimal fractional solution found.")
fprintf("Begin iterative rounding...\n")

%Parameters used to fix the rounded entries of b for iterative rounding.
ub = 100*ones(n,1);    %Upper bound
lb = zeros(n,1);       %Lower bound
nrange = 1:n;          %Needed for slicing
rb = round(b);

%Iterative rounding
while(sum(lb==0)>1)
  unknown_variables = size(lb(lb==0));   %Number of unfixed entires of b
  range = nrange(lb == 0);
  rounding_indices = range(abs(b(range)-rb(range)) == min(abs(b(range)-rb(range))));
  if(max(unknown_variables) == 74)
      fprintf('Rounding entry %2d of b. %2d entries remain to be rounded.', rounding_indices, max(unknown_variables)-1);
  elseif(max(unknown_variables) == 2)
      fprintf(repmat('\b',1,56))
      fprintf('Rounding entry %2d of b. %2d entry remains to be rounded.', rounding_indices, max(unknown_variables)-1);
  else
      fprintf(repmat('\b',1,56))
      fprintf('Rounding entry %2d of b. %2d entries remain to be rounded.', rounding_indices, max(unknown_variables)-1);
  end
      
  rb = round(b);
  ub(rounding_indices) = rb(rounding_indices); %Fix the rounded entries of
  lb(rounding_indices) = rb(rounding_indices); %b by upper and lower
                                               %bounding them.
                                               
  cvx_begin quiet
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

disp("All entries of b have been rounded.")


idx = find(lb-ub);         %Find the last unrounded index
b(idx) = ceil(b(idx));     %Round it up
b = round(b);              %Cast as integer
bhat=b;

fprintf("Making some final adjustments...")

%We decrease entries of the bound vector until the protocol is minimally-secure
while sum(log(2*bhat+1)/log(2)) >= 256
  max_idx = max(find(b == max(b)));
  bhat(max_idx) = bhat(max_idx)-1;
  if sum(log(2*bhat+1)/log(2)) >= 256
    b=bhat;
  end
end

fprintf(" Bound vector finalized.\n")

fprintf("Searching for strategies and permutations for each submeasure...")
%Search for a strategy and permutation
Batch = MultipleStrategies(Primes,mu,iota,b);

fprintf(" Strategies and permutations found.\n")
fprintf("Writing header file to simba_withdummy_2.h...")

file = fopen('simba_withdummy_2.h', 'w');

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

sizes = zeros(batches,1);
complement_sizes = zeros(batches,1);
for i = 1:batches
  sizes(i) = length(Batch{i});
  complement_sizes(i) = length(complement_Batch{i});
end

%Printing some constants
fprintf(file, '#ifndef _SIMBA_PARAMETERS_H_\n');
fprintf(file, '#define _SIMBA_PARAMETERS_H_\n');
fprintf(file, '\n');
fprintf(file, '// SIMBA-(NUMBER_OF_BATCHES)-MY\n');
fprintf(file, '#define NUMBER_OF_BATCHES %2d\n', batches); %Total number of batches (NOT batches per round)
fprintf(file, '#define MY %2d\n', 1); %Always 1
fprintf(file, '\n');
fprintf(file, '// (each entry corresponds to the number of degree-(l_i) to be required in the action: this the one given in Onuki et al. work)\n');

%Bound vector
fprintf(file, 'static int8_t B[] = {\n');
for i = 1:n-1
  fprintf(file, '%3d, ', B(i));
  if mod(i,8) == 0
    fprintf(file, '\n');
  end
end
fprintf(file, '%3d', B(n));
fprintf(file, '\n};\n\n');

%Batch primes
fprintf(file, '// (NUMBER_OF_BATCHES) different subsets (i.e., batches)\n');

%This code requires the batch entries to be listed in the order that the
%isogenies will be computed, rather than the order in which the primes
%will be multiplied in, so we reverse the order of each batch.
for i=1:batches
  fprintf(file, 'static uint8_t BATCH_%1d[] = {\n', i-1);
  for j=length(Batch{i}):-1:2
    fprintf(file, '%3d, ', Batch{i}(j));
    if(mod(j,8) == mod(length(Batch{i})+1,8))
      fprintf(file, '\n');
    end
  end
  fprintf(file, '%3d' ,Batch{i}(1));
  fprintf(file, '\n};\n\n');
end

%Batch sizes
fprintf(file, 'static uint8_t SIZE_OF_EACH_BATCH[NUMBER_OF_BATCHES] = {');
for i = 1:length(sizes)-1
  if(mod(i,8) == 1)
    fprintf(file, '\n');
  end
  fprintf(file, '%3d, ', sizes(i));
end
fprintf(file, '%3d', sizes(end));
fprintf(file, '\n};\n\n');

%Array of batches
fprintf(file, 'static uint8_t *BATCHES[NUMBER_OF_BATCHES] = {');
for i = 0:batches-2
  if(mod(i,4) == 0)
    fprintf(file, '\n');
  end
  fprintf(file, 'BATCH_%1d, ', i);
end
fprintf(file, 'BATCH_%1d', batches-1);
fprintf(file, '\n};\n\n');

%Last isogenies
fprintf(file, 'static uint8_t LAST_ISOGENY[NUMBER_OF_BATCHES] = {');
for i=1:batches-1
  if(mod(i,8) == 1)
    fprintf(file, '\n');
  end
  fprintf(file, '%2d, ', Batch{i}(1));
end
fprintf(file, '%2d', Batch{batches}(1));
fprintf(file, '\n};\n\n');

%Total number of isogenies
fprintf(file, 'static uint16_t NUMBER_OF_ISOGENIES = %3d;\n\n', sum(B));

%Batch complements/sizes
fprintf(file, '// The complement of each batch\n');

%fflush(file);

%Batch complement sizes
fprintf(file, 'static uint8_t SIZE_OF_EACH_COMPLEMENT_BATCH[NUMBER_OF_BATCHES] = {\n');
for i=1:batches-1
  fprintf(file, '%2d, ', complement_sizes(i));
  if(mod(i,8) == 0)
    fprintf(file, '\n');
  end
end
fprintf(file, '%2d', complement_sizes(end));
fprintf(file, '};\n\n');

%Batch complements
fprintf(file, 'static uint8_t COMPLEMENT_OF_EACH_BATCH[NUMBER_OF_BATCHES][N] = {\n');
fprintf(file, '{\n'); 
for i = 1:length(complement_Batch{1})-1
  fprintf(file, '%2d, ', complement_Batch{1}(i));
  if (mod(i,8) == 0)
    fprintf(file, '\n');
  end
end
fprintf(file, '%2d, \n', complement_Batch{1}(end));
for i = 1:sizes(1)-1
  if (mod(i,8) == 1)
    fprintf(file, '\n');
  end
  fprintf(file, '%2d, ', n);
end
fprintf(file, '%2d', n);
fprintf(file, '}');

for j = 2:batches
  fprintf(file, ',\n\n{\n');
  for i = 1:complement_sizes(j)-1
    fprintf(file, '%2d, ', complement_Batch{j}(i));
    if (mod(i,8) == 0)
      fprintf(file, '\n');
    end
  end
  fprintf(file, '%2d, \n', complement_Batch{j}(end));
  for i = 1:sizes(j)-1
    fprintf(file, '%2d, ', n);
    if (mod(i,8) == 0)
      fprintf(file, '\n');
    end
  end
fprintf(file, '%2d', n);
fprintf(file, '\n}');
end

%fflush(file);
fprintf(file, '\n};\n\n#endif\n');

%fflush(file);

fprintf(" Header file written.\n")
fprintf("Saving bound vector and batch information to FullParameterSet.mat\n")
save -mat7-binary 'FullParameterSet.mat' b Batch %Save output to 'FullParameterSet.mat'
