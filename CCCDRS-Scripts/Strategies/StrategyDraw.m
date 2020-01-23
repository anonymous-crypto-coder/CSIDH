function StrategyDraw(H,V)

%{
Draws a strategy
--------------------------------------------------------------------------------

Input:
------

H: Horiztonal strategy matrix
V: Vertical strategy matrix

Output:
-------

None.
%}

%This code is self-explanatory.

n = size(H)(1)+1;

figure();

for i=1:(n-1)
  for j=1:(n-1)
    if(H(i,j) == 1)
      plot([j;j+1], [n-i;n-i],'Color','k');
      hold on;
    end
    if(V(i,j) == 1)
      plot([j;j], [n-i;n-i+1],'Color','k');
      hold on;
    end
  end
end

set(gca,'visible','off');
