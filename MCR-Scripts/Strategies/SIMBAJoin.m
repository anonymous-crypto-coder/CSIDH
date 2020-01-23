function [H,V] = SIMBAJoin(HL,VL,HR,VR)

%{
Joins a left SIMBA substrategy to a right SIMBA substrategy.

--------------------------------------------------------------------------------

Input:
------

HL: The horizontal strategy matrix of the left SIMBA substrategy.
VL: The vertical strategy matrix of the left SIMBA substrategy.
HR: The horizontal strategy matrix of the right SIMBA substrategy.
VR: The vertical strategy matrix of the right SIMBA substrategy.

Output:
-------

H: The horizontal strategy matrix of the joined SIMBA strategy.
V: The vertical strategy matrix of the joined SIMBA strategy.

%}

%This is self-explanatory.

k=size(HL)(1);
l=size(HR)(1);

H = [[HR, zeros(l,k+1)];zeros(1,k+l+1);[zeros(k,l+1), HL]];
V = [[VR, zeros(l,k+1)];zeros(1,k+l+1);[zeros(k,l+1), VL]];
