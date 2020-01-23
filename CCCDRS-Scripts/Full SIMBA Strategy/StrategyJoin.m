function [H,V] = StrategyJoin(HL, VL, HR, VR)

%{
Joins a left substrategy to a right substrategy.

--------------------------------------------------------------------------------

Input:
------

HL: The horizontal strategy matrix of the left substrategy.
VL: The vertical strategy matrix of the left substrategy.
HR: The horizontal strategy matrix of the right substrategy.
VR: The vertical strategy matrix of the right substrategy.

Output:
-------

H: The horizontal strategy matrix of the joined strategy.
V: The vertical strategy matrix of the joined strategy.

%}

%This is self-explanatory.

k = 1 + size(HL)(1); 
l = 1 + size(HR)(1);

if k==1 && l==1
  H = [1];
  V = [1];
elseif k==1
  H = [ [HR, zeros(l-1,1)]; ones(1,l)];
  V = [ [VR, zeros(l-1,1)]; [1, zeros(1,l-1)]];
elseif l==1
  H = [ [zeros(k-1,1);1], [zeros(1,k-1); HL]];
  V = [ ones(k,1), [zeros(1,k-1); VL]];
else
  H = [ [HR, zeros(l-1,k)]; zeros(1,k+l-1); [ [zeros(k-2,l); ones(1,l)], HL]];
  V = [ [VR, zeros(l-1,k)]; [1, zeros(1,k+l-2)];[ones(k-1,1), zeros(k-1,l-1), VL] ];
end
