%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Problem 1a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = zeros(4,2); %store time

%%%%%%%%% initial point at (2,2) %%%%%%%%%
x0 = 2;
y0 = 2;
%time
t1 = 0; %optimal time
t2 = 0; %non-optimal time
%%%%%%%%%%%%%%%%%%%%%%%% optimal trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%blue circle centered at (-1,0)
th = linspace(0, 2*pi, 100);
R = norm([x0+1,y0]);  %radius
%circle
cx1 = R*cos(th) - 1;
cy1 = R*sin(th);
%find the first intersection of blue circle with switching curve
syms x y
%eyeball the switching curve it intersects with
eqns=[(x-3)^2+y^2==1,(x+1)^2+y^2==R^2];
sol = solve(eqns);
pt1 = [sol.x(1);sol.y(1)]
pt1 = double(pt1);
%first arc of optimal trajectory
th0 = atan(y0/(x0+1)); % initial angle for blue circle
th1 = atan(pt1(2)/(pt1(1)+1)); %first switch angle for blue circle
t1 = t1 + th0-th1; %time is the same as angle
th = linspace( th0, th1, 100);
%trajectory arc 1
trajx1 = R*cos(th) - 1;
trajy1 = R*sin(th);

%green circle centered at (1,0)
th = linspace(0, 2*pi, 100);
R = norm(pt1-[1;0]);  %radius
%circle
cx2 = R*cos(th) + 1;
cy2 = R*sin(th);
%find the next intersection of green circle with switching curve
syms x y
eqns=[(x+1)^2+y^2==1,(x-1)^2+y^2==R^2];
sol = solve(eqns);
pt2 = [sol.x(2);sol.y(2)]
pt2 = double(pt2);
%second arc of optimal trajectory
th0 = atan(pt1(2)/(pt1(1)-1)); % initial angle for green circle
th1 = acos((pt2(1)-1)/R); %first switch angle for green circle
t1 = t1 + th0+2*pi - th1; %time is the same as angle
th = linspace( th0+2*pi, th1, 100);
%trajectory arc 2
trajx2 = R*cos(th) + 1;
trajy2 = R*sin(th);

%semi-circle switching curves
th = linspace( 0, -pi, 100);
R = 1;  %radius
sx1 = R*cos(th) + 3;
sy1 = R*sin(th) ;
sx2 = R*cos(th) + 1;
sy2 = R*sin(th) ;
th = linspace( 0, pi, 100);
sx3 = R*cos(th) - 1;
sy3 = R*sin(th) ;
%third arc of optimal trajectory
th0 = asin(pt2(2)); % initial angle for semicircle
th1 = 0; %terminal angle for for semicircle
t1 = t1+ th0-th1;
th = linspace( th0,th1, 100);
trajx3 = R*cos(th) - 1;
trajy3 = R*sin(th);

figure;
plot(x0,y0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
hold on
plot(pt1(1),pt1(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(pt2(1),pt2(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(0,0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(sx1,sy1,Color=[1 .5 0],LineWidth=1.5); axis equal;
plot(sx2,sy2,Color=[1 .5 0],LineWidth=1.5); 
plot(sx3,sy3,Color=[1 .5 0],LineWidth=1.5); 
plot(cx1,cy1,Color='blue');
plot(trajx1,trajy1,Color='magenta',LineWidth=2);
plot(cx2,cy2,Color='green');
plot(trajx2,trajy2,Color='magenta',LineWidth=2);
plot(trajx3,trajy3,Color='magenta',LineWidth=2);
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
xlabel('wx1');
ylabel('wx2');
title('optimal trajectory wx1=2');
hold off

%%%%%%%%%%%%%%%%%%%%%%%% non-optimal trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%
%blue circle centered at (-1,0)
th = linspace(0, 2*pi, 100);
R = norm([x0+1,y0]);  %radius
cx1 = R*cos(th) - 1;
cy1 = R*sin(th);
%find the first intersection of blue circle with switching curve
syms x y
eqns=[y==0,(x+1)^2+y^2==R^2];
sol = solve(eqns);
pt1 = [sol.x(2);sol.y(2)]
pt1 = double(pt1);
%first arc of optimal trajectory
th0 = atan(y0/(x0+1)); % initial angle for blue circle
th1 = atan(pt1(2)/(pt1(1)+1)); %first switch angle for blue circle
t2 = t2 + th0-th1; %time is the same as angle
th = linspace( th0, th1, 100);
trajx1 = R*cos(th) - 1;
trajy1 = R*sin(th);

%green circle centered at (1,0)
th = linspace(0, 2*pi, 100);
R = norm(pt1-[1;0]);  %radius
cx2 = R*cos(th) + 1;
cy2 = R*sin(th);
%find the next intersection of green circle with switching curve
syms x y
eqns=[(x+1)^2+y^2==1,(x-1)^2+y^2==R^2];
sol = solve(eqns);
pt2 = [sol.x(2);sol.y(2)]
pt2 = double(pt2);
%second arc of non-optimal trajectory
th0 = atan(pt1(2)/(pt1(1)-1)); % initial angle for green circle
th1 = acos((pt2(1)-1)/R); %first switch angle for green circle
t2 = t2 + th0+2*pi - th1; %time is the same as angle
th = linspace( th0+2*pi, th1, 100);
trajx2 = R*cos(th) + 1;
trajy2 = R*sin(th);
%Gamma hat switching curves
th = linspace( 0, -pi, 100);
R = 1;  %radius
sx1 = linspace(-5,-2);
sx11 = linspace(2,5);
sy1 = 0*sx1;
sy11 = 0*sx11;
sx2 = R*cos(th) + 1;
sy2 = R*sin(th) ;
th = linspace( 0, pi, 100);
sx3 = R*cos(th) - 1;
sy3 = R*sin(th) ;
%third arc of optimal trajectory
th0 = acos(pt2(1)+1); % initial angle for semicircle
th1 = 0; %terminal angle for for semicircle
t2 = t2+ th0-th1;
th = linspace( th0,th1, 100);
trajx3 = R*cos(th) - 1;
trajy3 = R*sin(th);

figure;
plot(x0,y0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
hold on
plot(pt1(1),pt1(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(pt2(1),pt2(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(0,0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(sx1,sy1,Color=[1 .5 0],LineWidth=1.5); axis equal;
plot(sx11,sy11,Color=[1 .5 0],LineWidth=1.5);
plot(sx2,sy2,Color=[1 .5 0],LineWidth=1.5); 
plot(sx3,sy3,Color=[1 .5 0],LineWidth=1.5); 
plot(cx1,cy1,Color='blue');
plot(trajx1,trajy1,Color='magenta',LineWidth=2);
plot(cx2,cy2,Color='green');
plot(trajx2,trajy2,Color='magenta',LineWidth=2);
plot(trajx3,trajy3,Color='magenta',LineWidth=2);
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
xlabel('wx1');
ylabel('wx2');
title('non-optimal trajectory wx1=2');
hold off

fprintf('The optimal time is %f, the non-optimal time is %f.\n',t1,t2)

t(1,1) = t1;
t(1,2) = t2;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Problem 1b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% initial point at (3,2) %%%%%%%%%
x0 = 3;
y0 = 2;
%time
t1 = 0; %optimal time
t2 = 0; %non-optimal time
%blue circle centered at (-1,0)
th = linspace(0, 2*pi, 100);
R = norm([x0+1,y0]);  %radius
cx1 = R*cos(th) - 1;
cy1 = R*sin(th);
%find the first intersection of blue circle with switching curve
syms x y
eqns=[(x-3)^2+y^2==1,(x+1)^2+y^2==R^2];
sol = solve(eqns);
pt1 = [sol.x(1);sol.y(1)]
pt1 = double(pt1);
%first arc of optimal trajectory
th0 = atan(y0/(x0+1)); % initial angle for blue circle
th1 = atan(pt1(2)/(pt1(1)+1)); %first switch angle for blue circle
t1 = t1 + th0-th1; %time is the same as angle
th = linspace( th0, th1, 100);
trajx1 = R*cos(th) - 1;
trajy1 = R*sin(th);

%green circle centered at (1,0)
th = linspace(0, 2*pi, 100);
R = norm(pt1-[1;0]);  %radius
cx2 = R*cos(th) + 1;
cy2 = R*sin(th);
%find the next intersection of green circle with switching curve
syms x y
eqns=[(x+1)^2+y^2==1,(x-1)^2+y^2==R^2];
sol = solve(eqns);
pt2 = [sol.x(2);sol.y(2)]
pt2 = double(pt2);
%%%%%%second arc of optimal trajectory
th0 = atan(pt1(2)/(pt1(1)-1)); % initial angle for green circle
th1 = acos((pt2(1)-1)/R); %first switch angle for green circle
t1 = t1 + th0+2*pi - th1; %time is the same as angle
th = linspace( th0+2*pi, th1, 100);
trajx2 = R*cos(th) + 1;
trajy2 = R*sin(th);
%semi-circle switching curves
th = linspace( 0, -pi, 100);
R = 1;  %radius
sx1 = R*cos(th) + 3;
sy1 = R*sin(th) ;
sx2 = R*cos(th) + 1;
sy2 = R*sin(th) ;
th = linspace( 0, pi, 100);
sx3 = R*cos(th) - 1;
sy3 = R*sin(th) ;
%%%%third arc of optimal trajectory
th0 = acos(pt2(1)+1); % initial angle for semicircle
th1 = 0; %terminal angle for for semicircle
t1 = t1+ th0-th1;
th = linspace( th0,th1, 100);
trajx3 = R*cos(th) - 1;
trajy3 = R*sin(th);

figure;
plot(x0,y0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
hold on
plot(pt1(1),pt1(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(pt2(1),pt2(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(0,0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(sx1,sy1,Color=[1 .5 0],LineWidth=1.5); axis equal;
plot(sx2,sy2,Color=[1 .5 0],LineWidth=1.5); 
plot(sx3,sy3,Color=[1 .5 0],LineWidth=1.5); 
plot(cx1,cy1,Color='blue');
plot(trajx1,trajy1,Color='magenta',LineWidth=2);
plot(cx2,cy2,Color='green');
plot(trajx2,trajy2,Color='magenta',LineWidth=2);
plot(trajx3,trajy3,Color='magenta',LineWidth=2);
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
xlabel('wx1');
ylabel('wx2');
title('optimal trajectory wx1=3');
hold off

%%%%%%%%%%%%%%%%%%%%%%%% non-optimal trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%
%blue circle centered at (-1,0)
th = linspace(0, 2*pi, 100);
R = norm([x0+1,y0]);  %radius
cx1 = R*cos(th) - 1;
cy1 = R*sin(th);
%find the first intersection of blue circle with switching curve
syms x y
eqns=[y==0,(x+1)^2+y^2==R^2];
sol = solve(eqns);
pt1 = [sol.x(2);sol.y(2)]
pt1 = double(pt1);
%first arc of optimal trajectory
th0 = atan(y0/(x0+1)); % initial angle for blue circle
th1 = atan(pt1(2)/(pt1(1)+1)); %first switch angle for blue circle
t2 = t2 + th0-th1; %time is the same as angle
th = linspace( th0, th1, 100);
trajx1 = R*cos(th) - 1;
trajy1 = R*sin(th);

%green circle centered at (1,0)
th = linspace(0, 2*pi, 100);
R = norm(pt1-[1;0]);  %radius
cx2 = R*cos(th) + 1;
cy2 = R*sin(th);
%find the next intersection of green circle with switching curve
syms x y
eqns=[(x+1)^2+y^2==1,(x-1)^2+y^2==R^2];
sol = solve(eqns);
pt2 = [sol.x(2);sol.y(2)]
pt2 = double(pt2);
%second arc of non-optimal trajectory
th0 = atan(pt1(2)/(pt1(1)-1)); % initial angle for green circle
th1 = acos((pt2(1)-1)/R); %first switch angle for green circle
t2 = t2 + th0+2*pi - th1; %time is the same as angle
th = linspace( th0+2*pi, th1, 100);
trajx2 = R*cos(th) + 1;
trajy2 = R*sin(th);
%Gamma hat switching curves
th = linspace( 0, -pi, 100);
R = 1;  %radius
sx1 = linspace(-6,-2);
sx11 = linspace(2,6);
sy1 = 0*sx1;
sy11 = 0*sx11;
sy1 = 0*sx1;
sx2 = R*cos(th) + 1;
sy2 = R*sin(th) ;
th = linspace( 0, pi, 100);
sx3 = R*cos(th) - 1;
sy3 = R*sin(th) ;
%third arc of optimal trajectory
th0 = acos(pt2(1)+1); % initial angle for semicircle
th1 = 0; %terminal angle for for semicircle
t2 = t2+ th0-th1;
th = linspace( th0,th1, 100);
trajx3 = R*cos(th) - 1;
trajy3 = R*sin(th);

figure;
plot(x0,y0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
hold on
plot(pt1(1),pt1(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(pt2(1),pt2(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(0,0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(sx1,sy1,Color=[1 .5 0],LineWidth=1.5); axis equal;
plot(sx11,sy11,Color=[1 .5 0],LineWidth=1.5);
plot(sx2,sy2,Color=[1 .5 0],LineWidth=1.5); 
plot(sx3,sy3,Color=[1 .5 0],LineWidth=1.5); 
plot(cx1,cy1,Color='blue');
plot(trajx1,trajy1,Color='magenta',LineWidth=2);
plot(cx2,cy2,Color='green');
plot(trajx2,trajy2,Color='magenta',LineWidth=2);
plot(trajx3,trajy3,Color='magenta',LineWidth=2);
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
xlabel('wx1');
ylabel('wx2');
title('non-optimal trajectory wx1=3');
hold off


t(2,1) = t1;
t(2,2) = t2;

%%%%%%%%%%%%%%%%%%%%%%%% initial point at (4,2) %%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = 4;
y0 = 2;

%time
t1 = 0; %optimal time
t2 = 0; %non-optimal time
%blue circle centered at (-1,0)
th = linspace(0, 2*pi, 100);
R = norm([x0+1,y0]);  %radius
cx1 = R*cos(th) - 1;
cy1 = R*sin(th);
%find the first intersection of blue circle with switching curve
syms x y
eqns=[(x-5)^2+y^2==1,(x+1)^2+y^2==R^2];
sol = solve(eqns);
pt1 = [sol.x(1);sol.y(1)]
pt1 = double(pt1);
%first arc of optimal trajectory
th0 = atan(y0/(x0+1)); % initial angle for blue circle
th1 = atan(pt1(2)/(pt1(1)+1)); %first switch angle for blue circle
th0-th1
t1 = t1 + th0-th1; %time is the same as angle
th = linspace( th0, th1, 100);
trajx1 = R*cos(th) - 1;
trajy1 = R*sin(th);

%green circle centered at (1,0)
th = linspace(0, 2*pi, 100);
R = norm(pt1-[1;0]);  %radius
cx2 = R*cos(th) + 1;
cy2 = R*sin(th);
%find the next intersection of green circle with switching curve
syms x y
eqns=[(x+3)^2+y^2==1,(x-1)^2+y^2==R^2];
sol = solve(eqns);
pt2 = [sol.x(2);sol.y(2)]
pt2 = double(pt2);
%%%%%%second arc of optimal trajectory
th0 = atan(pt1(2)/(pt1(1)-1)); % initial angle for green circle
th1 = acos((pt2(1)-1)/R); %first switch angle for green circle
th0+2*pi - th1
t1 = t1 + th0+2*pi - th1; %time is the same as angle
th = linspace( th0+2*pi, th1, 100);
trajx2 = R*cos(th) + 1;
trajy2 = R*sin(th);

%blue circle centered at (-1,0)
th = linspace(0, 2*pi, 100);
R = norm(pt2+[1;0]);  %radius
cx3 = R*cos(th) - 1;
cy3 = R*sin(th);
%find the next intersection of blue circle with switching curve
syms x y
eqns=[(x-1)^2+y^2==1,(x+1)^2+y^2==R^2];
sol = solve(eqns);
pt3 = [sol.x(1);sol.y(1)]
pt3 = double(pt3);
%third arc of optimal trajectory
th0 = acos((pt2(1)+1)/R); % initial angle for blue circle
th1 = atan(pt3(2)/(pt3(1)+1)); %first switch angle for blue circle
th0-th1
t1 = t1 + th0-th1; %time is the same as angle
th = linspace( th0, th1, 100);
trajx3 = R*cos(th) - 1;
trajy3 = R*sin(th);

%semi-circle switching curves
th = linspace( 0, -pi, 100);
R = 1;  %radius
sx1 = R*cos(th) + 3;
sy1 = R*sin(th) ;
sx2 = R*cos(th) + 1;
sy2 = R*sin(th) ;
sx3 = R*cos(th) + 5;
sy3 = R*sin(th) ;
th = linspace( 0, pi, 100);
sx4 = R*cos(th) - 1;
sy4 = R*sin(th) ;
sx5 = R*cos(th) - 3;
sy5 = R*sin(th) ;
%%%%fourth arc of optimal trajectory
th0 = -acos(pt3(1)-1); % initial angle for semicircle
th1 = -pi; %terminal angle for for semicircle
abs(th1-th0)
t1 = t1+ abs(th1-th0);
th = linspace( th0,th1, 100);
trajx4 = R*cos(th) + 1;
trajy4 = R*sin(th);

figure;
plot(x0,y0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
hold on
plot(pt1(1),pt1(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(pt2(1),pt2(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(pt3(1),pt3(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(0,0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(sx1,sy1,Color=[1 .5 0],LineWidth=1.5); axis equal;
plot(sx2,sy2,Color=[1 .5 0],LineWidth=1.5); 
plot(sx3,sy3,Color=[1 .5 0],LineWidth=1.5); 
plot(sx4,sy4,Color=[1 .5 0],LineWidth=1.5); 
plot(sx5,sy5,Color=[1 .5 0],LineWidth=1.5); 
plot(cx1,cy1,Color='blue');
plot(trajx1,trajy1,Color='magenta',LineWidth=2);
plot(cx2,cy2,Color='green');
plot(trajx2,trajy2,Color='magenta',LineWidth=2);
plot(cx3,cy3,Color='blue');
plot(trajx3,trajy3,Color='magenta',LineWidth=2);
plot(trajx4,trajy4,Color='magenta',LineWidth=2);
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
xlabel('wx1');
ylabel('wx2');
title('optimal trajectory wx1=4');
hold off

%%%%%%%%%%%%%%%%%%%%%%%% non-optimal trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%
%blue circle centered at (-1,0)
th = linspace(0, 2*pi, 100);
R = norm([x0+1,y0]);  %radius
cx1 = R*cos(th) - 1;
cy1 = R*sin(th);
%find the first intersection of blue circle with switching curve
syms x y
eqns=[y==0,(x+1)^2+y^2==R^2];
sol = solve(eqns);
pt1 = [sol.x(2);sol.y(2)]
pt1 = double(pt1);
%first arc of optimal trajectory
th0 = atan(y0/(x0+1)); % initial angle for blue circle
th1 = atan(pt1(2)/(pt1(1)+1)); %first switch angle for blue circle
th0-th1
t2 = t2 + th0-th1; %time is the same as angle
th = linspace( th0, th1, 100);
trajx1 = R*cos(th) - 1;
trajy1 = R*sin(th);

%green circle centered at (1,0)
th = linspace(0, 2*pi, 100);
R = norm(pt1-[1;0]);  %radius
cx2 = R*cos(th) + 1;
cy2 = R*sin(th);
%find the next intersection of green circle with switching curve
syms x y
eqns=[y==0,(x-1)^2+y^2==R^2];
sol = solve(eqns);
pt2 = [sol.x(1);sol.y(1)]
pt2 = double(pt2);
%second arc of non-optimal trajectory
th0 = atan(pt1(2)/(pt1(1)-1)); % initial angle for green circle
th1 = acos((pt2(1)-1)/R); %first switch angle for green circle
th0+2*pi - th1
t2 = t2 + th0+2*pi - th1; %time is the same as angle
th = linspace( th0+2*pi, th1, 100);
trajx2 = R*cos(th) + 1;
trajy2 = R*sin(th);

%blue circle centered at (-1,0)
th = linspace(0, 2*pi, 100);
R = norm(pt2+[1;0]);  %radius
cx3 = R*cos(th) - 1;
cy3 = R*sin(th);
%find the next intersection of blue circle with switching curve
syms x y
eqns=[(x-1)^2+y^2==1,(x+1)^2+y^2==R^2];
sol = solve(eqns);
pt3 = [sol.x(1);sol.y(1)]
pt3 = double(pt3);
%third arc of non-optimal trajectory
th0 = acos((pt2(1)+1)/R); % initial angle for blue circle
th1 = atan(pt3(2)/(pt3(1)+1)); %first switch angle for blue circle
th0-th1
t2 = t2 + th0-th1; %time is the same as angle
th = linspace( th0, th1, 100);
trajx3 = R*cos(th) - 1;
trajy3 = R*sin(th);

%Gamma hat switching curves
th = linspace( 0, -pi, 100);
R = 1;  %radius
sx1 = linspace(-7,-2);
sx11 = linspace(2,7);
sy1 = 0*sx1;
sy11 = 0*sx11;
sy1 = 0*sx1;
sx2 = R*cos(th) + 1;
sy2 = R*sin(th) ;
th = linspace( 0, pi, 100);
sx3 = R*cos(th) - 1;
sy3 = R*sin(th) ;

%%%%fourth arc of non-optimal trajectory
th0 = -acos(pt3(1)-1); % initial angle for semicircle
th1 = -pi; %terminal angle for for semicircle
abs(th1-th0)
t2 = t2 + abs(th1-th0);
th = linspace( th0,th1, 100);
trajx4 = R*cos(th) + 1;
trajy4 = R*sin(th);

figure;
plot(x0,y0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
hold on
plot(pt1(1),pt1(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(pt2(1),pt2(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(pt3(1),pt3(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(0,0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(sx1,sy1,Color=[1 .5 0],LineWidth=1.5); axis equal;
plot(sx11,sy11,Color=[1 .5 0],LineWidth=1.5);
plot(sx2,sy2,Color=[1 .5 0],LineWidth=1.5); 
plot(sx3,sy3,Color=[1 .5 0],LineWidth=1.5); 
plot(cx1,cy1,Color='blue');
plot(trajx1,trajy1,Color='magenta',LineWidth=2);
plot(cx2,cy2,Color='green');
plot(cx3,cy3,Color='blue');
plot(trajx2,trajy2,Color='magenta',LineWidth=2);
plot(trajx3,trajy3,Color='magenta',LineWidth=2);
plot(trajx4,trajy4,Color='magenta',LineWidth=2);
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
xlabel('wx1');
ylabel('wx2');
title('non-optimal trajectory wx1=4');
hold off

t(3,1) = t1;
t(3,2) = t2;
%%%%%%%%%%%%%%%%%%%%%% initial point at (5,2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = 5;
y0 = 2;

%time
t1 = 0; %optimal time
t2 = 0; %non-optimal time
%blue circle centered at (-1,0)
th = linspace(0, 2*pi, 100);
R = norm([x0+1,y0]);  %radius
cx1 = R*cos(th) - 1;
cy1 = R*sin(th);
%find the first intersection of blue circle with switching curve
syms x y
eqns=[(x-5)^2+y^2==1,(x+1)^2+y^2==R^2];
sol = solve(eqns);
pt1 = [sol.x(1);sol.y(1)]
pt1 = double(pt1);
%first arc of optimal trajectory
th0 = atan(y0/(x0+1)); % initial angle for blue circle
th1 = atan(pt1(2)/(pt1(1)+1)); %first switch angle for blue circle
t1 = t1 + th0-th1; %time is the same as angle
th = linspace( th0, th1, 100);
trajx1 = R*cos(th) - 1;
trajy1 = R*sin(th);

%green circle centered at (1,0)
th = linspace(0, 2*pi, 100);
R = norm(pt1-[1;0]);  %radius
cx2 = R*cos(th) + 1;
cy2 = R*sin(th);
%find the next intersection of green circle with switching curve
syms x y
eqns=[(x+3)^2+y^2==1,(x-1)^2+y^2==R^2];
sol = solve(eqns);
pt2 = [sol.x(2);sol.y(2)]
pt2 = double(pt2);
%%%%%%second arc of optimal trajectory
th0 = atan(pt1(2)/(pt1(1)-1)); % initial angle for green circle
th1 = acos((pt2(1)-1)/R); %first switch angle for green circle
t1 = t1 + th0+2*pi - th1; %time is the same as angle
th = linspace( th0+2*pi, th1, 100);
trajx2 = R*cos(th) + 1;
trajy2 = R*sin(th);

%blue circle centered at (-1,0)
th = linspace(0, 2*pi, 100);
R = norm(pt2+[1;0]);  %radius
cx3 = R*cos(th) - 1;
cy3 = R*sin(th);
%find the next intersection of blue circle with switching curve
syms x y
eqns=[(x-1)^2+y^2==1,(x+1)^2+y^2==R^2];
sol = solve(eqns);
pt3 = [sol.x(1);sol.y(1)]
pt3 = double(pt3);
%third arc of optimal trajectory
th0 = acos((pt2(1)+1)/R); % initial angle for blue circle
th1 = atan(pt3(2)/(pt3(1)+1)); %first switch angle for blue circle
t1 = t1 + th0-th1; %time is the same as angle
th = linspace( th0, th1, 100);
trajx3 = R*cos(th) - 1;
trajy3 = R*sin(th);

%semi-circle switching curves
th = linspace( 0, -pi, 100);
R = 1;  %radius
sx1 = R*cos(th) + 3;
sy1 = R*sin(th) ;
sx2 = R*cos(th) + 1;
sy2 = R*sin(th) ;
sx3 = R*cos(th) + 5;
sy3 = R*sin(th) ;
th = linspace( 0, pi, 100);
sx4 = R*cos(th) - 1;
sy4 = R*sin(th) ;
sx5 = R*cos(th) - 3;
sy5 = R*sin(th) ;
%%%%fourth arc of optimal trajectory
th0 = -acos(pt3(1)-1); % initial angle for semicircle
th1 = -pi; %terminal angle for for semicircle
t1 = t1+ abs(th1-th0);
th = linspace( th0,th1, 100);
trajx4 = R*cos(th) + 1;
trajy4 = R*sin(th);

figure;
plot(x0,y0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
hold on
plot(pt1(1),pt1(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(pt2(1),pt2(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(pt3(1),pt3(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(0,0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(sx1,sy1,Color=[1 .5 0],LineWidth=1.5); axis equal;
plot(sx2,sy2,Color=[1 .5 0],LineWidth=1.5); 
plot(sx3,sy3,Color=[1 .5 0],LineWidth=1.5); 
plot(sx4,sy4,Color=[1 .5 0],LineWidth=1.5); 
plot(sx5,sy5,Color=[1 .5 0],LineWidth=1.5); 
plot(cx1,cy1,Color='blue');
plot(trajx1,trajy1,Color='magenta',LineWidth=2);
plot(cx2,cy2,Color='green');
plot(trajx2,trajy2,Color='magenta',LineWidth=2);
plot(cx3,cy3,Color='blue');
plot(trajx3,trajy3,Color='magenta',LineWidth=2);
plot(trajx4,trajy4,Color='magenta',LineWidth=2);
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
xlabel('wx1');
ylabel('wx2');
title('optimal trajectory wx1=5');
hold off

%%%%%%%%%%%%%%%%%%%%%%%% non-optimal trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%
%blue circle centered at (-1,0)
th = linspace(0, 2*pi, 100);
R = norm([x0+1,y0]);  %radius
cx1 = R*cos(th) - 1;
cy1 = R*sin(th);
%find the first intersection of blue circle with switching curve
syms x y
eqns=[y==0,(x+1)^2+y^2==R^2];
sol = solve(eqns);
pt1 = [sol.x(2);sol.y(2)]
pt1 = double(pt1);
%first arc of non-optimal trajectory
th0 = atan(y0/(x0+1)); % initial angle for blue circle
th1 = atan(pt1(2)/(pt1(1)+1)); %first switch angle for blue circle
t2 = t2 + th0-th1; %time is the same as angle
th = linspace( th0, th1, 100);
trajx1 = R*cos(th) - 1;
trajy1 = R*sin(th);

%green circle centered at (1,0)
th = linspace(0, 2*pi, 100);
R = norm(pt1-[1;0]);  %radius
cx2 = R*cos(th) + 1;
cy2 = R*sin(th);
%find the next intersection of green circle with switching curve
syms x y
eqns=[y==0,(x-1)^2+y^2==R^2];
sol = solve(eqns);
pt2 = [sol.x(1);sol.y(1)]
pt2 = double(pt2);
%second arc of non-optimal trajectory
th0 = atan(pt1(2)/(pt1(1)-1)); % initial angle for green circle
th1 = acos((pt2(1)-1)/R); %first switch angle for green circle
t2 = t2 + th0+2*pi - th1; %time is the same as angle
th = linspace( th0+2*pi, th1, 100);
trajx2 = R*cos(th) + 1;
trajy2 = R*sin(th);

%blue circle centered at (-1,0)
th = linspace(0, 2*pi, 100);
R = norm(pt2+[1;0]);  %radius
cx3 = R*cos(th) - 1;
cy3 = R*sin(th);
%find the next intersection of blue circle with switching curve
syms x y
eqns=[(x-1)^2+y^2==1,(x+1)^2+y^2==R^2];
sol = solve(eqns);
pt3 = [sol.x(1);sol.y(1)]
pt3 = double(pt3);
%third arc of non-optimal trajectory
th0 = acos((pt2(1)+1)/R); % initial angle for blue circle
th1 = atan(pt3(2)/(pt3(1)+1)); %first switch angle for blue circle
t2 = t2 + th0-th1; %time is the same as angle
th = linspace( th0, th1, 100);
trajx3 = R*cos(th) - 1;
trajy3 = R*sin(th);

%Gamma hat switching curves
th = linspace( 0, -pi, 100);
R = 1;  %radius
sx1 = linspace(-8,-2);
sx11 = linspace(2,8);
sy1 = 0*sx1;
sy11 = 0*sx11;
sy1 = 0*sx1;
sx2 = R*cos(th) + 1;
sy2 = R*sin(th) ;
th = linspace( 0, pi, 100);
sx3 = R*cos(th) - 1;
sy3 = R*sin(th) ;

%%%%fourth arc of non-optimal trajectory
th0 = -acos(pt3(1)-1); % initial angle for semicircle
th1 = -pi; %terminal angle for for semicircle
t2 = t2+ abs(th1-th0);
th = linspace( th0,th1, 100);
trajx4 = R*cos(th) + 1;
trajy4 = R*sin(th);

figure;
plot(x0,y0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
hold on
plot(pt1(1),pt1(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(pt2(1),pt2(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(pt3(1),pt3(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(0,0,'o',MarkerEdgeColor='r',MarkerFaceColor='r')
plot(sx1,sy1,Color=[1 .5 0],LineWidth=1.5); axis equal;
plot(sx11,sy11,Color=[1 .5 0],LineWidth=1.5);
plot(sx2,sy2,Color=[1 .5 0],LineWidth=1.5); 
plot(sx3,sy3,Color=[1 .5 0],LineWidth=1.5); 
plot(cx1,cy1,Color='blue');
plot(trajx1,trajy1,Color='magenta',LineWidth=2);
plot(cx2,cy2,Color='green');
plot(cx3,cy3,Color='blue');
plot(trajx2,trajy2,Color='magenta',LineWidth=2);
plot(trajx3,trajy3,Color='magenta',LineWidth=2);
plot(trajx4,trajy4,Color='magenta',LineWidth=2);
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
xlabel('wx1');
ylabel('wx2');
title('non-optimal trajectory wx1=5');
hold off

t(4,1) = t1;
t(4,2) = t2;

figure;
plot(t(:,2)./t(:,1)*100)
ylabel('percentage');
xlabel('wx1-1');
title('time ratio as wx1 increases');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Problem 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x y
eqns = [x==-1/2*(y+1/2)^2-1/2*(y+1/2),x==-2*y];
sol2 = solve(eqns);
sol2.x
sol2.y
eqns = [x==1/2*(y+1/2)^2-3/2*(y+1/2),x==-2*y];
sol2 = solve(eqns);
sol2.x
sol2.y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Problem 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = linspace(-5,5);
x0 = 0;
y0 = 1;
%up
Gamma1xf =@(t) 0.5*t.^2 + y0*t + x0;
Gamma1yf =@(t) t + y0;
%down
Gamma2xf =@(t) -0.5*t.^2 + y0*t + x0;
Gamma2yf =@(t) -t + y0;
Gamma1x = Gamma1xf(t);
Gamma1y = Gamma1yf(t);
Gamma2x = Gamma2xf(t);
Gamma2y = Gamma2yf(t);
figure;
plot(Gamma1x,Gamma1y,Color=[1 .5 0],LineWidth=1.5)
hold on
plot(Gamma2x,Gamma2y,Color=[1 .5 0],LineWidth=1.5)
for n=1:5
    k=2*n/3;
    plot(Gamma1x-k,Gamma1y,Color='cyan')
    plot(Gamma1x+k,Gamma1y,Color='cyan')
    plot(Gamma2x-k,Gamma2y,Color='cyan')
    plot(Gamma2x+k,Gamma2y,Color='cyan')
end
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
xlabel('x1');
ylabel('x2');
title('switching curves');
axis([-4,4,-4,4]);
hold off

%%%%%%%%%%% initial (0,0) %%%%%%%%%%%%
x0 = 0;
y0 = 0;
upx =@(t) 0.5*t.^2 + y0*t + x0;
upy =@(t) t + y0;
downx =@(t) -0.5*t.^2 + y0*t + x0;
downy =@(t) -t + y0;

%u=-1 goes down, intersects Gamma2 u=1 goes up
eqns = [x-(x0+1/2*y0^2)==-1/2*(y)^2,x+1/2==1/2*(y)^2];
sol31 = solve(eqns);
sol31.x
sol31.y
pt1 = [sol31.x(1);sol31.y(1)]; %pick the lowest x2 to go up
pt1 = double(pt1);
tf1 = -pt1(2);
t1 = linspace(0,tf1);

%u=1 goes up, intersects Gamma1 u=-1 goes down
eqns = [x-(x0-1/2*y0^2)==1/2*(y)^2,x-1/2==-1/2*(y)^2];
sol32 = solve(eqns);
sol32.x
sol32.y
pt2 = [sol32.x(1);sol32.y(2)]; %pick the highest x2 to go down
pt2 = double(pt2);
tf2 = pt2(2);
t2 = linspace(0,tf2);
figure;
plot(Gamma1x,Gamma1y,Color=[1 .5 0],LineWidth=1.5)
hold on
plot(Gamma2x,Gamma2y,Color=[1 .5 0],LineWidth=1.5)
plot(upx(t),upy(t),Color='cyan')
plot(downx(t),downy(t),Color='cyan')
plot(x0,y0,'o',MarkerEdgeColor='r',MarkerFaceColor='r',MarkerSize=6)
plot(pt1(1),pt1(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r',MarkerSize=3)
plot(pt2(1),pt2(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r',MarkerSize=3)
plot(downx(t1),downy(t1),Color='magenta',LineWidth=2)
plot(upx(t2),upy(t2),Color='green',LineWidth=2)
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
xlabel('x1');
ylabel('x2');
title('initial point (0,0)');
axis([-1,1,-1.1,1.1]);
hold off


%%%%%%%%%%% initial (1,0) %%%%%%%%%%%%
x0 = 1;
y0 = 0;
upx =@(t) 0.5*t.^2 + y0*t + x0;
upy =@(t) t + y0;
downx =@(t) -0.5*t.^2 + y0*t + x0;
downy =@(t) -t + y0;

%u=-1 goes down, intersects Gamma2 u=1 goes up
eqns = [x-(x0+1/2*y0^2)==-1/2*(y)^2,x+1/2==1/2*(y)^2];
sol31 = solve(eqns);
sol31.x
sol31.y
pt1 = [sol31.x(1);sol31.y(1)]; %pick the lowest x2 to go up
pt1 = double(pt1);
tf1 = -pt1(2);
t1 = linspace(0,tf1);

%u=1 goes up, intersects Gamma1 u=-1 goes down
eqns = [x-(x0-1/2*y0^2)==1/2*(y)^2,x-1/2==-1/2*(y)^2];
%no solution!
figure;
plot(Gamma1x,Gamma1y,Color=[1 .5 0],LineWidth=1.5)
hold on
plot(Gamma2x,Gamma2y,Color=[1 .5 0],LineWidth=1.5)
plot(upx(t),upy(t),Color='cyan')
plot(downx(t),downy(t),Color='cyan')
plot(x0,y0,'o',MarkerEdgeColor='r',MarkerFaceColor='r',MarkerSize=6)
plot(pt1(1),pt1(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r',MarkerSize=3)
plot(downx(t1),downy(t1),Color='magenta',LineWidth=2)
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
xlabel('x1');
ylabel('x2');
title('initial point (1,0)');
axis([-1,2,-2,2]);
hold off

%%%%%%%%%%% initial (-1,1) %%%%%%%%%%%%
x0 = -1;
y0 = 1;
upx =@(t) 0.5*t.^2 + y0*t + x0;
upy =@(t) t + y0;
downx =@(t) -0.5*t.^2 + y0*t + x0;
downy =@(t) -t + y0;

%u=-1 goes down, intersects Gamma2 u=1 goes up
eqns = [x-(x0+1/2*y0^2)==-1/2*(y)^2,x+1/2==1/2*(y)^2];
sol31 = solve(eqns);
sol31.x
sol31.y
pt1 = [sol31.x(1);sol31.y(1)]; %pick the lowest x2 to go up
pt1 = double(pt1);
tf1 = -pt1(2);
t1 = linspace(0,tf1);

%u=1 goes up, intersects Gamma1 u=-1 goes down
eqns = [x-(x0-1/2*y0^2)==1/2*(y)^2,x-1/2==-1/2*(y)^2];
sol32 = solve(eqns);
sol32.x
sol32.y
pt2 = [sol32.x(1);sol32.y(2)]; %pick the highest x2 to go down
pt2 = double(pt2);
tf2 = pt2(2);
t2 = linspace(0,tf2);
figure;
plot(Gamma1x,Gamma1y,Color=[1 .5 0],LineWidth=1.5)
hold on
plot(Gamma2x,Gamma2y,Color=[1 .5 0],LineWidth=1.5)
plot(upx(t),upy(t),Color='cyan')
plot(downx(t),downy(t),Color='cyan')
plot(x0,y0,'o',MarkerEdgeColor='r',MarkerFaceColor='r',MarkerSize=6)
plot(pt1(1),pt1(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r',MarkerSize=3)
plot(pt2(1),pt2(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r',MarkerSize=3)
plot(downx(t1),downy(t1),Color='magenta',LineWidth=2)
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
xlabel('x1');
ylabel('x2');
title('initial point (-1,1)');
axis([-2,1,-2,2]);
hold off

%%%%%%%%%%% initial (1,-1) %%%%%%%%%%%%
x0 = 1;
y0 = -1;

upx =@(t) 0.5*t.^2 + y0*t + x0;
upy =@(t) t + y0;
downx =@(t) -0.5*t.^2 + y0*t + x0;
downy =@(t) -t + y0;

%u=-1 goes down, intersects Gamma2 u=1 goes up
eqns = [x-(x0+1/2*y0^2)==-1/2*(y)^2,x+1/2==1/2*(y)^2];
sol31 = solve(eqns);
sol31.x
sol31.y
pt1 = [sol31.x(1);sol31.y(1)]; %pick the lowest x2 to go up
pt1 = double(pt1);
tf1 = -pt1(2);
t1 = linspace(0,tf1);

%u=1 goes up, intersects Gamma1 u=-1 goes down
eqns = [x-(x0-1/2*y0^2)==1/2*(y)^2,x-1/2==-1/2*(y)^2];
sol32 = solve(eqns);
sol32.x
sol32.y
pt2 = [sol32.x(1);sol32.y(1)]; %pick the highest x2 to go down
pt2 = double(pt2);
tf2 = pt2(2);
t2 = linspace(0,tf2);
figure;
plot(Gamma1x,Gamma1y,Color=[1 .5 0],LineWidth=1.5)
hold on
plot(Gamma2x,Gamma2y,Color=[1 .5 0],LineWidth=1.5)
plot(upx(t),upy(t),Color='cyan')
plot(downx(t),downy(t),Color='cyan')
plot(x0,y0,'o',MarkerEdgeColor='r',MarkerFaceColor='r',MarkerSize=6)
plot(pt1(1),pt1(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r',MarkerSize=3)
plot(pt2(1),pt2(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r',MarkerSize=3)
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
xlabel('x1');
ylabel('x2');
title('initial point (1,-1)');
axis([-1,2,-2,2]);
hold off

%%%%%%%%%%% initial (-1,0) %%%%%%%%%%%%
x0 = -1;
y0 = 0;
upx =@(t) 0.5*t.^2 + y0*t + x0;
upy =@(t) t + y0;
downx =@(t) -0.5*t.^2 + y0*t + x0;
downy =@(t) -t + y0;

%u=-1 goes down, intersects Gamma2 u=1 goes up
eqns = [x-(x0+1/2*y0^2)==-1/2*(y)^2,x+1/2==1/2*(y)^2];
sol31 = solve(eqns);
sol31.x
sol31.y
pt1 = [sol31.x(1);sol31.y(1)]; %pick the lowest x2 to go up
pt1 = double(pt1);
tf1 = -pt1(2);
t1 = linspace(0,tf1);

%u=1 goes up, intersects Gamma1 u=-1 goes down
eqns = [x-(x0-1/2*y0^2)==1/2*(y)^2,x-1/2==-1/2*(y)^2];
sol32 = solve(eqns);
sol32.x
sol32.y
pt2 = [sol32.x(1);sol32.y(2)]; %pick the highest x2 to go down
pt2 = double(pt2);
tf2 = pt2(2);
t2 = linspace(0,tf2);
figure;
plot(Gamma1x,Gamma1y,Color=[1 .5 0],LineWidth=1.5)
hold on
plot(Gamma2x,Gamma2y,Color=[1 .5 0],LineWidth=1.5)
plot(upx(t),upy(t),Color='cyan')
plot(downx(t),downy(t),Color='cyan')
plot(x0,y0,'o',MarkerEdgeColor='r',MarkerFaceColor='r',MarkerSize=6)
plot(pt1(1),pt1(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r',MarkerSize=3)
plot(pt2(1),pt2(2),'o',MarkerEdgeColor='r',MarkerFaceColor='r',MarkerSize=3)
plot(upx(t2),upy(t2),Color='green',LineWidth=2)
set(gca, 'XAxisLocation', 'origin')
set(gca, 'YAxisLocation', 'origin')
xlabel('x1');
ylabel('x2');
title('initial point (-1,0)');
axis([-2,1,-2,2]);
hold off
