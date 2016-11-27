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
OPTIONS = odeset('RelTol',1e-3,'AbsTol',1e-6);
func=@ml; %v func
[t,S]=ode15s(func, tspan, initial); % V differential eqn solver

%solving for V at equilibrium
syms Vnull;
syms wn;
null1 =(1/c)*(gcal*(0.5*(1+tanh((Vnull-v1)/v2)))*(vca-Vnull) + gk*(0.5*(1+tanh((Vnull-v3)/v4)))*(vk-Vnull) + gl*(vl-Vnull)+Iext);
Veq = solve(null1,Vnull);
Weq = 0.5*(1+tanh((Veq-v3)/v4));

%evaluating the jacobian and eigenvalues
jac = jacobian([(1/c)*(gcal*(0.5*(1+tanh((Vnull-v1)/v2)))*(vca-Vnull) + gk*wn*(vk-Vnull) + gl*(vl-Vnull)+Iext), phi*((0.5*(1+tanh((Vnull-v3)/v4)))-wn)*cosh((Vnull-v3)/v4)],[Vnull,wn]);
jaceq = subs(jac, {sym('Vnull'), sym('wn')}, {Veq, Weq})
eigenv = eig(jaceq);
eigen_stab(eigenv(1),eigenv(2))

figure(1)
%plotting the V and w nullclines
myfun1 = @(Vn,wn) (1/c)*(gcal*(0.5*(1+tanh((Vn-v1)/v2)))*(vca-Vn) + gk*wn*(vk-Vn) + gl*(vl-Vn)+Iext);
myfun2 = @(Vn,wn) phi*((0.5*(1+tanh((Vn-v3)/v4)))-wn)/(1/cosh((Vn-v3)/v4));
set(gca, 'fontsize', 14);
a2=ezplot(@(Vn,wn) myfun2(Vn,wn), [-80 40 -0.2 1]); %setting realistic limits on v,w
hold on
a1=ezplot(@(Vn,wn) myfun1(Vn,wn), [-80 60 -0.2 1]);
set(a1,'color',[1 0 0])
title('Nullclines')

[tp1,Sp1]=ode15s(func, [0 300], [-12, 0.014]); %generating action potential

%plot trajectories
plot(Sp1(:,1), Sp1(:,2));

phi=0.04;

[tp2,Sp2]=ode15s(func, [0 300], [-12, 0.014]); %generating another action potential

%plot trajectories
plot(Sp2(:,1), Sp2(:,2),'m');

legend('W nullcline', 'V nullcline','\phi=0.02','\phi=0.04');
title('Nullclines and Trajectories for Iext=0');
plot(Veq,Weq,'o');

figure(2)
phi=0.02
Iext=86.0; %New case with Iext 86uA/cm^2
myfun1 = @(Vn,wn) (1/c)*(gcal*(0.5*(1+tanh((Vn-v1)/v2)))*(vca-Vn) + gk*wn*(vk-Vn) + gl*(vl-Vn)+Iext);
myfun2 = @(Vn,wn) phi*((0.5*(1+tanh((Vn-v3)/v4)))-wn)/(1/cosh((Vn-v3)/v4));

set(gca, 'fontsize', 14);
%plotting nullclines again
a3=ezplot(@(Vn,wn) myfun2(Vn,wn), [-80 40 -0.2 1])
set(a3,'Linewidth', 3);
hold on
a4=ezplot(@(Vn,wn) myfun1(Vn,wn), [-80 60 -0.2 1])
set(a4,'color',[1 0 0],'Linewidth', 3);

title('Nullclines for Iext=86')

syms Vnull1;
syms wn1;
%evaluating new equilibrium point
null2 =(1/c)*(gcal*0.5*(1+tanh((Vnull1-v1)/v2))*(vca-Vnull1) + gk*0.5*(1+tanh((Vnull1-v3)/v4))*(vk-Vnull1) + gl*(vl-Vnull1)+Iext);

Veq1 = solve(null2,Vnull1)
Weq1 = 0.5*(1+tanh((Veq1-v3)/v4))

[tp3,Sp3]=ode15s(func, [0 500], [-60.8,0.0149]);
[tp4,Sp4]=ode15s(func, [0 500], [-27.9, 0.17]);
[tp5,Sp5]=ode15s(func, [0 500], [-27.9524,0.1195]); %eqlbrm point manually added
[tp6,Sp6]=ode15s(func, [0 -500], [-27.9, 0.17]);

p3=plot(Sp3(:,1), Sp3(:,2));
set(p3,'color',[0.5 0.2 0.2]);
hold on
p4=plot(Sp4(:,1), Sp4(:,2));
set(p4,'color',[0.8 0.8 0]);
hold on
p5=plot(Sp5(:,1), Sp5(:,2));
set(p5,'color',[0 1 1]);
hold on
p6=plot(Sp6(:,1), Sp6(:,2),'y');
set(p6,'color',[0.2 0.5 0.8]);
legend('W nullcline','V nullcline','start==eqlbrm for iext=0','start==[-27.9,0.17]','start==eqlbrm for iext=86','UPO for negative time, start==[-27.9,0.17]');
plot(Veq1,Weq1,'o');

%evaluating jacobian and eigenvalues
jac1 = jacobian([(1/c)*(gcal*(0.5*(1+tanh((Vnull1-v1)/v2)))*(vca-Vnull1) + gk*wn1*(vk-Vnull1) + gl*(vl-Vnull1)+Iext), phi*((0.5*(1+tanh((Vnull1-v3)/v4)))-wn1)*cosh((Vnull1-v3)/(v4))],[Vnull1,wn1]);
jaceq1 = subs(jac1, {sym('Vnull1'), sym('wn1')}, {Veq1, Weq1})
eigenv1 = eig(jaceq1)
eigen_stab(eigenv1(1),eigenv1(2));

%question no 10
gcal=4; % mS/ cm^2
gk=8.0; % mS/ cm^2
gl=2; % mS/ cm^2
vca=120; % mV
vk=-84; % mV
vl=-60; % mV
phi=0.667; %ms^-1
v1=-1.2; % mV
v2=18; % mV
v3=2; % mV
v4=12; % mV
v5=17.4; % mV
v6=12; % mV
c=20; % uF/cm^2
Iext=30; % uA/cm^2

syms Vnull2;
syms wn2;
%evaluating new equilibrium point
null3 =(1/c)*(gcal*(0.5*(1+tanh((Vnull2-v1)/v2)))*(vca-Vnull2) + gk*(0.5*(1+tanh((Vnull2-v3)/v4)))*(vk-Vnull2) + gl*(vl-Vnull2)+Iext);
Veq2_p1=vpasolve(null3,Vnull2,-10); %adding initial guesses to get 3 points
Veq2_p2=vpasolve(null3,Vnull2,-20)
Veq2_p3=vpasolve(null3,Vnull2,-40);
Weq2_p1=0.5*(1+tanh((Veq2_p1-v3)/v4));
Weq2_p2=0.5*(1+tanh((Veq2_p2-v3)/v4))
Weq2_p3=0.5*(1+tanh((Veq2_p3-v3)/v4));

%evaluating the jacobians and eigenvalues of each equilibrium point
jac2 = jacobian([(1/c)*(gcal*(0.5*(1+tanh((Vnull-v1)/v2)))*(vca-Vnull) + gk*wn*(vk-Vnull) + gl*(vl-Vnull)+Iext), phi*((0.5*(1+tanh((Vnull-v3)/v4)))-wn)*cosh((Vnull-v3)/v4)],[Vnull,wn]);
jaceq2 = subs(jac2, {sym('Vnull'), sym('wn')}, {Veq2_p1, Weq2_p1});
eigenv_1 = eig(jaceq2);
eigen_stab(eigenv_1(1),eigenv_1(2))
jaceq2 = subs(jac2, {sym('Vnull'), sym('wn')}, {Veq2_p2, Weq2_p2});
eigenv_2 = eig(jaceq2);
eigen_stab(eigenv_2(1),eigenv_2(2))
jaceq2 = subs(jac2, {sym('Vnull'), sym('wn')}, {Veq2_p3, Weq2_p3});
eigenv_3 = eig(jaceq2);
eigen_stab(eigenv_3(1),eigenv_3(2))



figure(3)
myfun1 = @(Vn,wn) (1/c)*(gcal*(0.5*(1+tanh((Vn-v1)/v2)))*(vca-Vn) + gk*wn*(vk-Vn) + gl*(vl-Vn)+Iext);
myfun2 = @(Vn,wn) phi*((0.5*(1+tanh((Vn-v3)/v4)))-wn)/(1/cosh((Vn-v3)/v4));
set(gca, 'fontsize', 14);

%plotting nullclines again
a5=ezplot(@(Vn,wn) myfun2(Vn,wn), [-80 40 -0.2 1]);
set(a5,'Linewidth',2,'color',[0.5 0 0.5]);
hold on
a6=ezplot(@(Vn,wn) myfun1(Vn,wn), [-80 60 -0.2 1]);
set(a6,'Linewidth',2,'color',[0.5 0.5 0])
title('Nullclines,Manifolds for Iext=30')


%evaluating manifolds
[tp7,Sp7]=ode15s(func, [0 1000], [-19.242, 0.0281]);
[tp8,Sp8]=ode15s(func, [0 -1000], [-19.242, 0.0281]);
[tp9,Sp9]=ode15s(func, [0 1000], [-19.244, 0.0283]);
[tp10,Sp10]=ode15s(func, [0 -1000], [-19.244, 0.0283]);
plot(Sp7(:,1), Sp7(:,2),'b');
hold on
plot(Sp8(:,1), Sp8(:,2),'m');
hold on
plot(Sp9(:,1), Sp9(:,2),'b');
hold on 
plot(Sp10(:,1), Sp10(:,2),'m');
legend('W nullcline','V nullcline','stable manifold','unstable manifold');
plot(Veq2_p1,Weq2_p1,'o')
plot(Veq2_p2,Weq2_p2,'o')
plot(Veq2_p3,Weq2_p3,'o')
hold off
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

function eigen_stab(p1,p2)

if isreal(p1)
    if p1*p2>0
        if p1<0
        display('The point')
        display(p1)
        display(',')
        display(p2)
        display('is stable')
        end
        if p1>0
        display('The point')
        display(p1)
        display(',')
        display(p2)
        display('is stable')
        end
    end
    if p1*p2<0
        display('The point')
        display(p1)
        display(',')
        display(p2)
        display('is a saddle node')
    end
else 
    if real(p1)>0
        display('The point')
        display(p1)
        display(',')
        display(p2)
        display('is unstable and diverging spiral')
    end
if real(p1)<0
    display('The point')
    display(p1)
    display(',')
    display(p2)
    display('is stable and converging spiral')
end
end
end
