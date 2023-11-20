%Problem 2(b)
solinit=bvpinit(linspace(0,1),[1,1,-2.2],1.2); %tau is in [0,1]
sol2=bvp4c(@odehw52,@bchw52,solinit);

x=sol2.y(1,:);
y=sol2.y(2,:);
theta = sol2.y(3,:);
tf= sol2.parameters;
time=tf*sol2.x;

figure
plot(x,y)
xlabel('x');
ylabel('y');
title('optimal trajectory');

function dydt = odehw52(t,y,tf)
dydt = tf*[cos(y(3))-0.2*(3*y(2)-y(2)^3)
           sin(y(3))
           0.6*(1-y(2)^2)*cos(y(3))^2];
end

function res = bchw52(ya,yb,tf)
res = [ya(1)-1
       yb(1)
       ya(2)-1 
       yb(2)];
end

%Problem 4(a)
syms theta t
eqn=[(15*cos(theta)+2)*t-20==-15,(15*sin(theta)-6)*t==35.5];
sol4a=solve(eqn);

tf=double(sol4a.t);
tf=tf(1); %choose positive solution
theta = double(sol4a.theta);
theta = theta(1);

time=linspace(0,tf);
x=(15*cos(theta)+2)*time-20;
y=(15*sin(theta)-6)*time;
figure
plot(x,y)
hold on

[X, Y] = meshgrid(-21:2:-14, -2:4:36);
u = 2*ones(size(X));
v = -6*ones(size(Y));

% Normalize vectors for better visualization
magnitude = sqrt(u.^2 + v.^2);
u = u ./ magnitude;
v = v ./ magnitude;

quiver(X, Y, u, v);
hold off
xlabel('x');
ylabel('y');
title('optimal trajectory');

%Problem 4(b)
syms tf theta
x=@(t) (15*cos(theta)+2)*t-20;
y=@(t) (15*sin(theta)-6)*t;
eqn=[theta==atan(1/(0.25+0.006*x(tf)^2)),25-0.25*x(tf) ...
    -0.002*x(tf)^3-y(tf)==0];
sol4b = solve(eqn);
theta4b = double(sol4b.theta);
tf4b = double(sol4b.tf);

x=@(t) (15*cos(theta4b)+2)*t-20;
y=@(t) (15*sin(theta4b)-6)*t;
xf4b = x(tf4b);
yf4b = y(tf4b);

f=@(x) 25-0.25*x-0.002*x.^3; 
xplot = linspace(-20,2);

[X, Y] = meshgrid(-21:2:0, -2:4:36);
u = 2*ones(size(X));
v = -6*ones(size(Y));
% Normalize vectors for better visualization
magnitude = sqrt(u.^2 + v.^2);
u = u ./ magnitude;
v = v ./ magnitude;
time=linspace(0,tf4b);

figure
plot(x(time),y(time))
xlabel('x');
ylabel('y');
title('optimal trajectory');
hold on
plot(xplot,f(xplot))
quiver(X, Y, u, v);
hold off
legend('trajectory','shoreline')

%Problem 4(c)
tinit = linspace(0,1);
solinit=bvpinit(tinit,@guess54,3); %tau is in [0,1]
sol4c=bvp4c(@odehw54,@bchw54,solinit);

x=sol4c.y(1,:);
y=sol4c.y(2,:);
theta = atan(sol4c.y(3,:));
tf= sol4c.parameters;
time=tf*sol4c.x;
f=@(x) 25-0.25*x-0.002*x.^3; 
xplot = linspace(-20,30);

[X, Y] = meshgrid(-18:2:26, -30:4:10);
u = -(Y-50);
v = 2*(X-15);
% Normalize vectors for better visualization
magnitude = sqrt(u.^2 + v.^2);
u = u ./ magnitude;
v = v ./ magnitude;

figure
plot(x,y)
xlabel('x');
ylabel('y');
title('optimal trajectory');
hold on
plot(xplot,f(xplot))
quiver(X, Y, u, v);
hold off
legend('trajectory','shoreline')

figure
plot(time,theta)
xlabel('time');
ylabel('theta');
title('optimal control');

function dydt = odehw54(t,y,tf)
dydt = tf*[15*cos(y(3))-y(2)+50 %x
           15*sin(y(3))+2*(y(1)-15) %y
           1+sin(y(3))^2]; %theta
end

function res = bchw54(ya,yb,tf)
res = [ya(1)+20 %x0
       ya(2) %y0
       25-0.25*yb(1)-0.002*yb(1)^3-yb(2) %endpoint constraint
       tan(yb(3))-1/(0.25+0.006*yb(1)^2)]; %transversality costates
end

function g=guess54(t)
xinit = @(t) -20+t;
yinit = @(t) -t;
g=[xinit(t) yinit(t) 1];
end

%Problem 5

syms p22(t) u0 tf
sol5=dsolve(diff(p22,t)==2*u0+p22^2,p22(tf)==0);