function [ variables ] = deriv_hh_defect(t,states,I_ext,f_ni)
global G_Na G_K G_l E_Na E_K E_l Q T Cm T_celcius
V=states(1);
m=states(2);
h=states(3);
n=states(4);

g_Na_nm=G_Na*(1-f_ni)*m^3*h;
g_Na_m=G_Na*f_ni*m^3;
I_Na=g_Na_nm*(V-E_Na)+g_Na_m*(V-E_Na);

g_K=G_K*n^4;
I_K=g_K*(V-E_K);

g_l=G_l;
I_l=g_l*(V-E_l);

I_ion=I_Na+I_K+I_l;

phi=Q^((T_celcius-6.3)/10);

alpha_m=-(0.1*phi*(V+35))/(exp(-(V+35)/10)-1);
beta_m=4*phi*exp(-(V+60)/18);

alpha_h=0.07*phi*exp(-(V+60)/20);
beta_h=phi/(exp(-(V+30)/10)+1);

alpha_n=-(0.01*phi*(V+50))/(exp(-(V+50)/10)-1);
beta_n=0.125*phi*exp(-(V+60)/80);

D_m=alpha_m*(1-m)-beta_m*m;
D_h=alpha_h*(1-h)-beta_h*h;
D_n=alpha_n*(1-n)-beta_n*n;
D_V=(I_ext-I_ion)/Cm;

variables=[D_V;D_m;D_h;D_n];


