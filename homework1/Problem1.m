clear all

method = ('im');
% method = ('ex'); %comment for changing the method

% Number of nodes
N = 3;
ndof = N * 2; % number of degrees of freedom

dti = 0.01; % second - Time step size
dte = 0.00001; %explicit

RodLength = 0.1; % l = 1 m
deltaL = RodLength / (N-1);

% Radii of spheres
R1 = 0.005; % meter
R2 = 0.025; % meter
R3 = 0.005; % meter

R1 = 0.005; % meter
R2 = 0.005; % meter
R3 = 0.005; % meter

% Density
rho_metal = 7000; % kg/m^3
rho_f = 1000; % fluid
rho = rho_metal - rho_f;
r0 = 0.001; % meter - rod radius
Y = 1e9; % Young's modulus (Y instead of E for clarity)
g = 9.8; % m/s^2 - gravity
visc = 1000; % pa-s
totalTime = 10.1; % second - total simulation time

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
M(1,1) = 4/3*pi*R1^3*rho_metal;
M(2,2) = 4/3*pi*R1^3*rho_metal;
M(3,3) = 4/3*pi*R2^3*rho_metal;
M(4,4) = 4/3*pi*R2^3*rho_metal;
M(5,5) = 4/3*pi*R3^3*rho_metal;
M(6,6) = 4/3*pi*R3^3*rho_metal;

% Viscous damping matrix, C
C = zeros(6,6);
C1 = 6 * pi * visc * R1;
C2 = 6 * pi * visc * R2;
C3 = 6 * pi * visc * R3;
C(1,1) = C1;
C(2,2) = C1;
C(3,3) = C2;
C(4,4) = C2;
C(5,5) = C3;
C(6,6) = C3;

% Weight vector, W
W = zeros(ndof, 1);
W(2) = -4/3*pi*R1^3*rho*g;
W(4) = -4/3*pi*R2^3*rho*g;
W(6) = -4/3*pi*R3^3*rho*g;

% Initial DOF
q0 = zeros(ndof, 1);

for c=1:N % loop over nodes
    q0( 2*c-1 ) = nodes(c,1); % x1, x2, x3
    q0( 2*c ) = nodes(c,2); % y1, y2, y3
end

q_init = q0;

u0 = zeros(ndof, 1); % old velocity (initial velocity)

% tolerance
tol = EI/RodLength^2 * 1e-3; % small enouch force that can be neglected

switch method

    case 'im'

% Time marching scheme
Nsteps = round(totalTime/dti);

% Storage for y-velocity of the middle node
all_mid_v = zeros(Nsteps, 1);
all_mid_y = zeros(Nsteps, 1);

for c = 2:Nsteps
 
     fprintf('Time = %f\n', (c-1) * dti);
     
     % Guess
     q = q0; % New DOFs are initialized to be equal to old DOFs
     % Newton Raphson
     err = 10 * tol;
     while err > tol
         f = M / dti * ( (q-q0)/dti - u0 );
         J = M / dti^2;
    
         % Linear spring between nodes 1 and 2
         xk = q(1);
         yk = q(2);
         xkp1 = q(3);
         ykp1 = q(4);
         l_k = deltaL;
         dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
         dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
         f(1:4) = f(1:4) + dF;
         J(1:4,1:4) = J(1:4,1:4) + dJ;
    
         % Linear spring between nodes 2 and 3
         xk = q(3);
         yk = q(4);
         xkp1 = q(5);
         ykp1 = q(6);
         l_k = deltaL;
         dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
         dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
         f(3:6) = f(3:6) + dF;
         J(3:6, 3:6) = J(3:6, 3:6) + dJ;
    
         % Bending spring at node 2
         xkm1 = q(1);
         ykm1 = q(2);
         xk = q(3);
         yk = q(4);
         xkp1 = q(5);
         ykp1 = q(6);
         curvature0 = 0;
         l_k = deltaL;
         dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
         dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
         f(1:6) = f(1:6) + dF;
         J(1:6, 1:6) = J(1:6, 1:6) + dJ;
    
             % Viscous force
         f = f + C * (q-q0) / dti;
         J = J + C / dti;
         
         % Weight
         f = f - W;
         
         % Update
         q = q - J \ f;
         
         err = sum ( abs(f) );
         
     end

    % New velocity
    u = (q - q0) / dti;
    
    % Store some information
    all_mid_v(c) = u(4);
    all_mid_y(c) = q(4);
     
     
     % Plot
%      figure(1);
%      plot( q(1:2:end), q(2:2:end), 'ro-');
%      axis equal
%      xlabel('x [meter]');
%      ylabel('y [meter]');
%      drawnow

     %plot the shape of the structure at t
     timenow = (c-1) * dti;
     if timenow == 0.01 || timenow == 0.05|| timenow == 0.1|| timenow == 1|| timenow == 10
        figure
        plot( q(1:2:end), q(2:2:end), 'o-');
        axis equal
        title(sprintf('strcuture at t = %.2f', timenow));
        xlabel('x [meter]');
        ylabel('y [meter]');
     end
     
     % Update (new becomes old)
     q0 = q;
     u0 = u; 
end


timeArray = (1:Nsteps) * dti;


% Plot middle node downward velocity
figure(2);
plot(timeArray, all_mid_v, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity (vertical) of middle node, v [m/s]');

% Plot middle node downward position
figure(3);
plot(timeArray, all_mid_y, 'k-');
xlabel('Time, t [sec]');
ylabel('Position (vertical) of middle node, v [m/s]');

%plot at t=0
figure(4)
plot( q_init(1:2:end), q_init(2:2:end), 'o-');
axis equal
title(sprintf('strcuture at t = 0'));
xlabel('x [meter]');
ylabel('y [meter]');


    case 'ex'

% Time marching scheme
Nsteps = round(totalTime/dte);

% Storage for y-velocity of the middle node
all_mid_v = zeros(Nsteps, 1);
all_mid_y = zeros(Nsteps, 1);

%calc for dt./M
invm = 1./M;
invm(isinf(invm)) = 0;


for c = 2:Nsteps
 
     fprintf('Time = %f\n', (c-1) * dte);

    f = zeros(ndof, 1);     

    % Linear spring between nodes 1 and 2
    xk = q0(1);
    yk = q0(2);
    xkp1 = q0(3);
    ykp1 = q0(4);
    l_k = deltaL;
    dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
    f(1:4) = f(1:4) + dF;


    % Linear spring between nodes 2 and 3
    xk = q0(3);
    yk = q0(4);
    xkp1 = q0(5);
    ykp1 = q0(6);
    l_k = deltaL;
    dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
    f(3:6) = f(3:6) + dF;
     
    % Bending spring at node 2
    xkm1 = q0(1);
    ykm1 = q0(2);
    xk = q0(3);
    yk = q0(4);
    xkp1 = q0(5);
    ykp1 = q0(6);
    curvature0 = 0;
    l_k = deltaL;
    dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
    f(1:6) = f(1:6) + dF;

    f_ext = C*u0 - W;

    %new position q
    q = q0 + u0*dte - (dte^2)*invm*(f + f_ext);

    % New velocity
     u = (q - q0) / dte;
     
     % Store some information
     all_mid_v(c) = u(4);
     all_mid_y(c) = q(4);

%      % Plot
%      figure(1);
%      plot( q(1:2:end), q(2:2:end), 'ro-');
%      axis equal
%      xlabel('x [meter]');
%      ylabel('y [meter]');
%      drawnow

     %plot the shape of the structure at t
     timenow = (c-1) * dte;
     if timenow == 0.01 || timenow == 0.05|| timenow == 0.1|| timenow == 1|| timenow == 10
        figure
        plot( q(1:2:end), q(2:2:end), 'o-');
        axis equal
        title(sprintf('strcuture at t = %.2f', timenow));
        xlabel('x [meter]');
        ylabel('y [meter]');
     end

     %update
     u0 = u;
     q0 = q;
end

timeArray = (1:Nsteps) * dti;

% Plot middle node downward velocity
figure(2);
plot(timeArray, all_mid_v, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity (vertical) of middle node, v [m/s]');

% Plot middle node downward position
figure(3);
plot(timeArray, all_mid_y, 'k-');
xlabel('Time, t [sec]');
ylabel('Position (vertical) of middle node, v [m/s]');

%plot at t=0
figure(4)
plot( q_init(1:2:end), q_init(2:2:end), 'o-');
axis equal
title(sprintf('strcuture at t = 0'));
xlabel('x [meter]');
ylabel('y [meter]');


end

