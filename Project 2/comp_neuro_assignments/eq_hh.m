syms V m h n 

% alpha_m=-1*(0.1*(V+35))/(exp(-(V+35)/10)-1);
% beta_m=4*exp(-(V+60)/18);
% 
% alpha_h=0.07*exp(-(V+60)/20);
% beta_h=1/(exp(-(V+30)/10)+1);
% 
% alpha_n=-(0.01*(V+50))/(exp(-(V+50)/10)-1);
% beta_n=0.125*exp(-(V+60)/80);
eq_pts=[];
for i=8:12
eqn1=i-120*m^3*h*(V-55)-36*n^4*(V+72)-0.3*(V+49)==0;
eqn2=(-1*(0.1*(V+35))/(exp(-(V+35)/10)-1))*(1-m)-(4*exp(-(V+60)/18))*m==0;
eqn3=(0.07*exp(-(V+60)/20))*(1-h)-(1/(exp(-(V+30)/10)+1))*h==0;
eqn4=(-(0.01*(V+50))/(exp(-(V+50)/10)-1))*(1-n)-(0.125*exp(-(V+60)/80))*n==0;
S=solve(eqn1,eqn2,eqn3,eqn4);
new=[S.V S.m S.h S.n];
eq_pts=[eq_pts; new];
end


