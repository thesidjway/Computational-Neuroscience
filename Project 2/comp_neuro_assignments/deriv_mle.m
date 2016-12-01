function [ variables ] = deriv_mle(t,states,I_ext)
V=states(1);
w=states(1);
global g_Ca g_K g_l V_Ca V_K V_l V1 V2 V3 V4 C

m_infi=(1+tanh((V-V1)/V2))/2;
I_Ca=g_Ca*m_infi*(V-V_Ca);

w_infi=(1+tanh((V-V3)/V4))/2;
tau_w=sech((V-V3)/(2*V4));
I_K=g_K*w*(V-V_K);

I_l=g_l(V-V_l);

I_ion=I_Ca+I_K+I_l;

D_V=(I_ext-I_ion)/C;
D_w=phi*(w_infi-w)/tau_w;

variables=[D_V; D_w];




