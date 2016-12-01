function [ variables ] = deriv_hh(t,states,I_ext)
global G_Na G_K G_l E_Na E_K E_l Q T Cm T_celcius 
V=states(1);
n=states(2);
h=1-n;

alpha_m=-(0.1*(V+35))/(exp(-(V+35)/10)-1);
beta_m=4*exp(-(V+60)/18);

m=alpha_m/(alpha_m+beta_m);


g_Na=G_Na*m^3*h;
I_Na=g_Na*(V-E_Na);

g_K=G_K*n^4;
I_K=g_K*(V-E_K);

g_l=G_l;
I_l=g_l*(V-E_l);

I_ion=I_Na+I_K+I_l;




alpha_n=-(0.01*(V+50))/(exp(-(V+50)/10)-1);
beta_n=0.125*exp(-(V+60)/80);

D_n=alpha_n*(1-n)-beta_n*n;
D_V=(I_ext-I_ion)/Cm;

variables=[D_V;D_n];

return
