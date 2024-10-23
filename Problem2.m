clear all

% Number of nodes
N = 21;
dt = 0.01; % second - Time step size
final_v = [];

%use for test number of nodes and timestep influence on the result
% for N = 1:2:51%
% dtArray = [0.01, 0.1, 1, 5, 10, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 35, 40 ,50];
% for dt = dtArray


ndof = N * 2; % number of degrees of freedom

RodLength = 0.1; %
deltaL = RodLength / (N-1);

% Radii of spheres
R = zeros(N,1); % Vector of size N - Radius of N nodes
R(:) = deltaL/10;
midNode = (N+1)/2;
R(midNode) = 0.025;

% Density
rho_metal = 7000; % kg/m^3
rho_f = 1000; % fluid
rho = rho_metal - rho_f;
r0 = 0.001; % meter - rod radius
Y = 1e9; % Young's modulus (Y instead of E for clarity)
g = 9.8; % m/s^2 - gravity
visc = 1000 / 100; % pa-s
totalTime = 50; % second - total simulation time

% Utility parameter
ne = N - 1; % number of edges
EI = Y * pi * r0^4 / 4; % Nm^2 - bending stiffness
EA = Y * pi * r0^2; % Newton

% Geometry - initial configuration
nodes = zeros(N,2);

for c=1:N % Loop over all the nodes
 nodes(c,1) = (c-1) * deltaL; % x coordinates
 nodes(c,2) = 0;
end

% Mass, M
M = zeros(ndof, ndof);
for k=1:N
 M(2*k-1, 2*k-1) = 4/3*pi*R(k)^3*rho_metal; % Mass for x_k
 M(2*k, 2*k) = M(2*k-1, 2*k-1); % Mass for y_k
end

% Viscous damping matrix, C
C = zeros(ndof,ndof);
for k=1:N
 C(2*k-1, 2*k-1) = 6 * pi * visc * R(k);
 C(2*k, 2*k) = C(2*k-1, 2*k-1);
end

% Weight vector, W
W = zeros(ndof, 1);
for k=1:N
 W(2*k-1) = 0; % weight along x is zero
 W(2*k) = -4/3*pi*R(k)^3*rho*g;
end

% Initial DOF
q0 = zeros(ndof, 1);
for c=1:N % loop over nodes
 q0( 2*c-1 ) = nodes(c,1); % x1, x2, x3
 q0( 2*c ) = nodes(c,2); % y1, y2, y3
end
u0 = zeros(ndof, 1); % old velocity (initial velocity)

% tolerance
tol = EI/RodLength^2 * 1e-3; % small enouch force that can be neglected

% Time marching scheme
Nsteps = round(totalTime/dt);

% Storage for y-velocity of the middle node
all_mid_v = zeros(Nsteps, 1);
all_mid_y = zeros(Nsteps, 1);


for c = 2:Nsteps
 
 fprintf('Time = %f\n', (c-1) * dt);
 
 
 % Guess
 q = q0; % New DOFs are initialized to be equal to old DOFs

 % Newton Raphson
 err = 10 * tol;
 while err > tol
 f = M / dt * ( (q-q0)/dt - u0 );
 J = M / dt^2;
 
 %
 % Elastic forces
 %
 
 % Linear spring
 for k=1:N-1
 xk = q(2*k-1);
 yk = q(2*k);
 xkp1 = q(2*k+1);
 ykp1 = q(2*k+2);
 l_k = deltaL;
 dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
 dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
 ind = [2*k-1, 2*k, 2*k+1, 2*k+2];
 f(ind) = f(ind) + dF;
 J(ind,ind) = J(ind,ind) + dJ;
 end
 
 % Bending spring
 for k=2:N-1 
 xkm1 = q(2*k-3);
 ykm1 = q(2*k-2);
 xk = q(2*k-1);
 yk = q(2*k);
 xkp1 = q(2*k+1);
 ykp1 = q(2*k+2);
 curvature0 = 0;
 l_k = deltaL;
 dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
 dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
 ind = [2*k-3, 2*k-2, 2*k-1, 2*k, 2*k+1, 2*k+2];
 f(ind) = f(ind) + dF;
 J(ind, ind) = J(ind, ind) + dJ;
 end
 
 % Viscous force
 f = f + C * (q-q0) / dt;
 J = J + C / dt;
 
 % Weight
 f = f - W;
 
 % Update
 q = q - J \ f;
 err = sum ( abs(f) );
 
 end
 
 % New velocity
 u = (q - q0) / dt;
 
 % Store some information
 all_mid_v(c) = u(2*midNode);
 all_mid_y(c) = q(2*midNode);

%  Plot
% if c==N
%  figure(1);
%  plot( q(1:2:end), q(2:2:end), 'ro-');
%  axis equal
%  xlabel('x [meter]');
%  ylabel('y [meter]');
% end

 
 % Update (new becomes old)
 q0 = q;
 u0 = u; 
end

% %use for test number of nodes and timestep influence on the result
% final_v = [final_v u(2*midNode)];
% 
% end

% Plot middle node downward velocity
% figure(2);
% timeArray = (1:Nsteps) * dt;
% plot(timeArray, all_mid_v, 'k-');
% xlabel('Time, t [sec]');
% ylabel('Velocity (vertical) of middle node, v [m/s]');
% 
% figure(3);
% plot(timeArray, all_mid_y, 'k-');
% xlabel('Time, t [sec]');
% ylabel('Position (vertical) of middle node, q [m/s]');

% plot for number of nodes or timestep vs final velocity
% NArray = 1:2:51;
% figure(4);
% plot(NArray, final_v, 'k-');
% xlabel('Node number, N');
% ylabel('Final velocity, q [m/s]');
% 
% figure(5);
% plot(dtArray, final_v, 'k-');
% xlabel('time step, dt[s]');
% ylabel('Final velocity, q [m/s]');
