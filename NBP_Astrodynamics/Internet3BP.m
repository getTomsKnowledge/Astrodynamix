clc 
clear all
close all
% Time step [s]
step = 0.001 ;  
%Body Mass 
m1 = 1; 
m2 = 1; 
m3 = 1; 
M = [m1 ; m2 ; m3] ; 
% Newton constant
G = 1;
% Time vector
t = 0:step:5; 
% Initial condiition were taken from:
%{
Numerical Search for Periodic Solutions in the Vicinity of the
Figure-Eight Orbit: Slaloming around Singularities on the Shape
Sphere
Milovan ˇSuvakov 
%C. Moore, Phys. Rev. Lett. 70, 3675 (1993).
%C. Sim ́o, Contemp Math 292 20928 (2002)
%}
% Initial conditions -- positions
x10 = -1;
y10 = 0;
z10 = 0 ;
 
x20 = -x10; 
y20 = -y10 ;  
z20 = 0 ; 
x30 = 0 ; 
y30 = 0 ; 
z30 = 0 ; 
% Initial condition -- velocity 
vx10 =  0.3471; 
vy10 =  0.5327;  
vz10 = 0 ; 
vx20 = vx10; 
vy20 = vy10;
vz20 = 0 ; 
vx30 = -2*vx10; 
vy30 = -2*vy10; 
vz30 = 0 ; 
% Preparing the vector to store the solution 
% There are 3 columns since I want to simulate three bodies 
x = zeros(numel(t),3) ; 
y = zeros(numel(t),3) ;
z = zeros(numel(t),3) ;  
vx = zeros(numel(t),3); 
vy = zeros(numel(t),3); 
vz = zeros(numel(t),3); 
% Assign the initial conditions
x(1,1) = x10 ;
x(1,2) = x20 ; 
x(1,3) = x30 ; 
y(1,1) = y10 ;
y(1,2) = y20 ; 
y(1,3) = y30 ; 
z(1,1) = z10 ;
z(1,2) = z20 ; 
z(1,3) = z30 ;
vx(1,1) = vx10  ;
vx(1,2) = vx20  ; 
vx(1,3) = vx30  ; 
vy(1,1) = vy10 ;
vy(1,2) = vy20 ; 
vy(1,3) = vy30 ; 
vz(1,1) = vz10 ;
vz(1,2) = vz20 ; 
vz(1,3) = vz30 ; 
%% Compute the Force
Fx = @(dx,dy,dz,mass) (-dx.*mass*G)./(sqrt (dx.^2 + dy.^2 + dz.^2) ).^3  ; 
Fy = @(dx,dy,dz,mass) (-dy.*mass*G) ./ ( sqrt (dx.^2 + dy.^2 + dz.^2 ) ).^3 ; 
Fz = @(dx,dy,dz,mass) (-dz.*mass*G) ./ ( sqrt (dx.^2 + dy.^2 + dz.^2 ) ).^3 ; 
%% Equation obtained by applying the simple Euler method 
% x(t+step) = x(t) + step * vx(t) 
% y(t+step) = y(t) + step * vy(t) 
% z(t+step) = z(t) + step * vz(t)
 
% vx(t+step) = vx(t) +step * fx(t) 
% vy(t+step) = vy(t) +step * fy(t) 
% vz(t+step) = vz(t) +step * fz(t)
% loop over each time step
for i = 1:length(t) 
    
	 % computing the postion at i+1
         x(i+1,:) = x(i,:) + step*vx(i,:) ;
         y(i+1,:) = y(i,:) + step*vy(i,:) ;
         z(i+1,:) = z(i,:) + step*vz(i,:) ;
    
	 % loop over the number of body 
         for j = 1:3
      
  
        dx = ( x(i,j)-x(i,:) )' ; 
        dy = ( y(i,j)-y(i,:) )' ; 
	    dz = ( z(i,j)-z(i,:) )' ; 
        
        Fxtot =  Fx (dx,dy,dz,M )  ; 
        Fytot =  Fy ( dx,dy,dz, M )  ; 
	    Fztot =  Fz ( dx,dy,dz, M )  ; 
                
        % change NaN into zeros
	% the NaN comes when computing the interaction of a body with itself
	% since (dx,dy,dz) = 0 
  
        Fxtot ( isnan(Fxtot)) = 0  ; 
        Fytot ( isnan(Fytot)) = 0 ; 
        Fztot ( isnan(Fztot)) = 0 ; 
 
	% total force exerted on the body i from all the other       
        Fxtot = sum(Fxtot)  ; 
        Fytot = sum(Fytot)  ;  
	    Fztot = sum(Fztot)  ; 
        
        
        
     
     % Computing the velocity at i+1
        vx(i+1,j)  = vx(i,j) + step * Fxtot;
        vy(i+1,j)  = vy(i,j) + step * Fytot;
        vz(i+1,j)  = vz(i,j) + step * Fztot;
        
	end
    
    
    
end
figure()
view([45 45 45 ]) 
hold on 
plot3(x(:,1),y(:,1),z(:,1))
plot3(x(1,1),y(1,1),z(1,1),'r*')
plot3( x(:,2) ,y(:,2),z(:,2)) 
plot3(x(2,2),y(2,2),z(2,2),'r*')
plot3( x(:,3) ,y(:,3),z(:,3))
plot3(x(3,3),y(3,3),z(3,3),'r*')
xlabel('x')
ylabel('y')
zlabel('z')
legend('B1','B1 start','B2','B2 start','B3','B3 start')
grid on
