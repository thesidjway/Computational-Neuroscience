
function[T,Y] = call_osc()
tspan =[0 300];
y1_0 = 1;
y2_0 = 0;
[T,Y] = ode15s(@osci,[0 60],[y1_0 y2_0]);
[T1,Y1] = ode15s(@osci1,[0 100],[y1_0 y2_0]);
[T2,Y2] = ode15s(@osci2,[0 300],[y1_0 y2_0]);

figure
subplot(3,1,1)
plot(T,Y(:,1))
title('u=0.1')
subplot(3,1,2)
plot(T1,Y1(:,1))
title('u=1')
subplot(3,1,3)
plot(T2,Y2(:,1))
title('u=100')

figure
subplot(3,1,1)
plot(Y(:,1),Y(:,2))
title('u=0.1')
subplot(3,1,2)
plot(Y1(:,1),Y1(:,2))
title('u=1')
subplot(3,1,3)
plot(Y2(:,1),Y2(:,2))
title('u=100')




end

function dydt =osci(t,y)
dydt =[y(2) ;  0.1*(1-y(1)^2)*y(2) - y(1)];
end


function dydt =osci1(t,y)
dydt =[y(2) ;  1*(1-y(1)^2)*y(2) - y(1)];
end

function dydt =osci2(t,y)
dydt =[y(2) ;  100*(1-y(1)^2)*y(2) - y(1)];
end
