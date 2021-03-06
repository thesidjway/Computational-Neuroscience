syms V n
figure
fun1=@(V,n) 10-120*((-1*(0.1*(V+35))/(exp(-(V+35)/10)-1))/((-1*(0.1*(V+35))/(exp(-(V+35)/10)-1))+(4*exp(-(V+60)/18))))^3*0.4*(V-55)-36*n^4*(V+72)-0.3*(V+49);
fun2=@(V,n) (-(0.01*(V+50))/(exp(-(V+50)/10)-1))*(1-n)-(0.125*exp(-(V+60)/80))*n; 


set(gca, 'fontsize', 14);
a1=ezplot(@(V,n) fun1(V,n), [-80 60 0 1])
hold on
a2=ezplot(@(V,n) fun2(V,n), [-80 60 -1 1])
set(a1,'color',[1 0 0])
title('Nullclines')

eqn1= 10-120*((-1*(0.1*(V+35))/(exp(-(V+35)/10)-1))/((-1*(0.1*(V+35))/(exp(-(V+35)/10)-1))+(4*exp(-(V+60)/18))))^3*0.4*(V-55)-36*n^4*(V+72)-0.3*(V+49)==0;
eqn2=(-(0.01*(V+50))/(exp(-(V+50)/10)-1))*(1-n)-(0.125*exp(-(V+60)/80))*n==0;

 [solV,soln]=solve(eqn1,eqn2);

%  ezplot(eqn1,[-80 60 -1 1]);
%   hold on
%   ezplot(eqn2,[-80 60 -1 1]);
