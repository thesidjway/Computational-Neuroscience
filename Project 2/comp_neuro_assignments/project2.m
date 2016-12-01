function project2
global c;
global gcal;
global vca;
global gk;
global vk;
global gl;
global vl;
global v1;
global v2;
global v3;
global v4;
global v5;
global v6;
global phi;
global Iext;
global Weq;
global Veq;

gcal=4.4; % mS/ cm^2
gk=8.0; % mS/ cm^2
gl=2; % mS/ cm^2
vca=120; % mV
vk=-84; % mV
vl=-60; % mV
phi=0.02; %ms^-1
v1=-1.2; % mV
v2=18; % mV
v3=2; % mV
v4=30; % mV
v5=2; % mV
v6=30; % mV
c=20; % uF/cm^2
Iext=0; % uA/cm^2

initial=[-60 0.015]; %initial conditions on v and w
tspan=[0 100]; %time period for simulation

func=@ml; %v func
[t,S]=ode15s(func, tspan, initial); % V differential 2n solver

syms Vnull;
syms wn;
null1 =(1/c)*(gcal*(0.5*(1+tanh((Vnull-v1)/v2)))*(vca-Vnull) + gk*(0.5*(1+tanh((Vnull-v3)/v4)))*(vk-Vnull) + gl*(vl-Vnull)+Iext);

Veq = solve(null1,Vnull)
Weq = 0.5*(1+tanh((Veq-v3)/v4))

jac = jacobian([(1/c)*(gcal*(0.5*(1+tanh((Vnull-v1)/v2)))*(vca-Vnull) + gk*wn*(vk-Vnull) + gl*(vl-Vnull)+Iext), phi*((0.5*(1+tanh((Vnull-v3)/v4)))-wn)/(1/cosh((Vnull-v3)/v4))],[Vnull,wn])

jaceq = subs(jac, {sym('Vnull'), sym('wn')}, {Veq, Weq})
eigenv = eig(jaceq);


figure(1)
myfun1 = @(Vn,wn) (1/c)*(gcal*(0.5*(1+tanh((Vn-v1)/v2)))*(vca-Vn) + gk*wn*(vk-Vn) + gl*(vl-Vn)+Iext);
myfun2 = @(Vn,wn) phi*((0.5*(1+tanh((Vn-v3)/v4)))-wn)/(1/cosh((Vn-v3)/v4));
set(gca, 'fontsize', 14);
a2=ezplot(@(Vn,wn) myfun2(Vn,wn), [-80 40 -0.2 1])
hold on
a1=ezplot(@(Vn,wn) myfun1(Vn,wn), [-80 60 -0.2 1])
set(a1,'color',[1 0 0])
title('Nullclines')

[tp3,Sp3]=ode15s(func, [0 300], [0  , 0.014]);
[tp4,Sp4]=ode15s(func, [0 300], [-12, 0.014]);
[tp5,Sp5]=ode15s(func, [0 300], [-24, 0.014]);

%plot trajectories
plot(Sp3(:,1), Sp3(:,2),'b',Sp4(:,1), Sp4(:,2),'b',Sp5(:,1), Sp5(:,2),'b');

phi=0.04;

[tp3,Sp3]=ode15s(func, [0 300], [0  , 0.014]);
[tp4,Sp4]=ode15s(func, [0 300], [-12, 0.014]);
[tp5,Sp5]=ode15s(func, [0 300], [-24, 0.014]);

%plot trajectories
plot(Sp3(:,1), Sp3(:,2),'m',Sp4(:,1), Sp4(:,2),'m',Sp5(:,1), Sp5(:,2),'m');
legend('W nullcline', 'V nullcline');
title('Nullclines and Trajectories. \phi=0.02 is blue and \phi=0.04 is magenta')


end

function dS = ml(t,S)

global c;
global gcal;
global vca;
global gk;
global vk;
global gl;
global vl;
global v1;
global v2;
global v3;
global v4;
global v5;
global v6;
global phi;
global Iext;

V=S(1);
w=S(2);

m_inf = (0.5*(1+tanh((V-v1)/v2)));
w_inf = (0.5*(1+tanh((V-v3)/v4)));
tau = 1/cosh((V-v3)/v4);

ddt_V = (1/c)*(-gcal*m_inf*(V - vca) - gk*w*(V - vk) - gl*(V - vl)+Iext);
ddt_w = phi*(w_inf-w)/(tau);

dS=[ddt_V; ddt_w];


end