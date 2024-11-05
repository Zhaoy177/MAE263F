close all;
clear all;
clc

%% global variables
global Fg M dt
global kappaBar EI GJ voronoiLength
global EA refLen

%% Inputs
% (1) Define the number of DOF related quantities
nv = 50; % number of nodes or vertices
ne = nv - 1; % number of edges
ndof = 3*nv + ne; % number od DOG

% (2) Define time step related quantities
dt = 0.01; % time step size in second
totalTime = 5; % simulation time
Nsteps = round(totalTime / dt);

% (3) Physical parameters
% (3a) Geometry
RodLength = 0.2; % meter
natR = 0.02; % meter
r0 = 0.001; % meter

% (3b) Material parameters
Y = 10e6; % Young's modulus
nu = 0.5; % Poisson's ratio
G = Y / (2 * (1+nu)); % Shear modulus
rho = 1000; % density (kg/m^3)

% (3c) Geometry
g = [0;0;-9.81];

%% Stifness variables
EI = Y * pi * r0^4 / 4; % Bending stiffness
GJ = G * pi * r0^4/2; % Shearing stiffness
EA = Y * pi * r0^2; % Stretching stiffness

%% Tolerance
tol = EI / RodLength^2 * 1e-6;


%% Mass Matrix
totalM = pi*r0^2*RodLength*rho; % kg
dm = totalM / ne; % mass per edge
massVector = zeros(ndof, 1);
for c = 1:nv
 ind = [4*c-3; 4*c-2; 4*c-1]; % c-th node
 if c==1 || c==nv
 massVector(ind) = dm/2;
 else
 massVector(ind) = dm;
 end
end
for c=1:ne
 ind = 4*c;
 massVector(ind) = 1/2 * dm * r0^2;
end
M = diag(massVector);

%% Initial DOF vector
nodes = zeros(nv, 3);
dTheta = (RodLength / natR) * (1/ne);
for c=1:nv
 nodes(c,1) = natR * cos ( (c-1) * dTheta );
 nodes(c,2) = natR * sin ( (c-1) * dTheta );
 nodes(c,3) = 0;
end
q0 = zeros(ndof, 1); % Initial condition
for c=1:nv
 ind = [4*c-3; 4*c-2; 4*c-1]; % c-th node
 q0(ind) = nodes(c,:);
end
% All theta angles are zero, i.e., q0(4:4:end) = 0 but it's redundant
u = zeros(ndof, 1); % Velocity

%% Reference length for each edge (used in stretching force)
refLen = zeros(ne, 1);
for c=1:ne % loop over each edge
 dx = nodes(c+1,:) - nodes(c,:); % dx = x_{c+1} - x_c
 refLen(c) = norm(dx);
end

%% Voronoi length (length associated with each node; used for bending and twisting)
voronoiLength = zeros(nv, 1);
for c=1:nv
 if c==1
 voronoiLength(c) = 1/2 * refLen(c);
 elseif c==nv
 voronoiLength(c) = 1/2 * refLen(c-1);
 else
 voronoiLength(c) = 1/2 * refLen(c-1) + 1/2 * refLen(c);
 end
end

%% Reference frame (At t=0, initialize using space parallel transport)
a1 = zeros(ne, 3); % First reference director for all the edges
a2 = zeros(ne, 3); % Second reference director for all the edges
tangent = computeTangent( q0 ); % Tangent for all the edges; size is (ne,3)

% Compute a1 for the first edge
t0 = tangent(1, :); % tangent on the first edge
t1 = [0;0;-1]; % arbitrary
a1Tmp = cross(t0, t1); % This vector is perp. to t0
if abs(a1Tmp) < 1e-6
 t1 = [0;1;0]; % arbitrary
 a1Tmp = cross(t0, t1); % This vector is perp. to t0
end
a1(1,:) = a1Tmp / norm(a1Tmp);
a2(1,:) = cross(tangent(1,:), a1(1,:));
% Done with the first edge
for c=2:ne
 t0 = tangent(c-1,:); % previous tangent on c-1-th edge
 t1 = tangent(c,:); % current tangent on c-th edge
 a1_0 = a1(c-1,:); % previous a1 director on c-1-th edge
 a1_1 = parallel_transport( a1_0, t0, t1);
 a1(c,:) = a1_1 / norm( a1_1 ); 
 a2(c,:) = cross(t1, a1(c,:));
end

%% Material frame
theta = q0(4:4:end); % Vector of size ne
[m1, m2] = computeMaterialDirectors( a1, a2, theta);

%% Reference twist
refTwist = zeros(nv, 1);
% refTwist = computeRefTwist(a1, tangent, refTwist);

%% Natural curvature
kappaBar = getkappa( q0, m1, m2 ); % Natural curvature at each node

%% Gravity
Fg = zeros(ndof, 1); % Gravity vector
for c=1:nv % loop over the nodes
 ind = [4*c-3; 4*c-2; 4*c-1]; % c-th node
 Fg(ind) = massVector(ind) .* g; % vector of size 3
end

%% Fixed and free DOFs
fixedIndex = 1:7; % Clamped BC
freeIndex = 8:ndof;

%% Time stepping scheme
ctime = 0; % Current time
endZ = zeros(Nsteps, 1); % z-coordinate of the last node with time
for timeStep = 1:Nsteps
 fprintf('Current time = %f\n', ctime);
 [q, u, a1, a2] = objfun(q0, u, a1, a2, freeIndex, ...
 tol, refTwist);
 ctime = ctime + dt;
 
 
 % Update q0 (old position)
 q0 = q;
 
 % Store endZ
 endZ(timeStep) = q(end);
 
 % Plot
%  if mod(timeStep,100) == 0
%  theta = q(4:4:end);
%  [m1, m2] = computeMaterialDirectors( a1, a2, theta);
%  plotrod(q, a1, a2, m1, m2, ctime);
%  end
end

% Visualization
figure(2);
timeArray = (1:1:Nsteps) * dt;
plot( timeArray, endZ, 'ro-');
xlabel('Time, t [sec]');
ylabel('z-coordinate of last node, \delta_z [m]');


%% compute tangent
function tangent = computeTangent(q)
% tangent is a matrix of size ne x 3
% q is a vector of size 4*nv-1
nv = (length(q)+1)/4;
ne = nv - 1;
tangent = zeros(ne, 3);
% Each column of "tangent" matrix is a 3-element
% vector representing the tangent on each edge
for c = 1:ne
 xc = q(4*c-3: 4*c-1);
 xcp1 = q(4*c+1: 4*c+3);
 edge = xcp1 - xc;
 tangent(c,:) = edge / norm(edge);
end
end

%% parallel transport
function d = parallel_transport(u, t1, t2)
% This function parallel transports a vector u from
% tangent t1 to t2.
b = cross(t1,t2);
if norm(b) == 0
 d = u;
else
 b = b / norm(b);
 % Good practice for safety
 b = b - dot(b,t1)*t1;
 b = b / norm(b);
 b = b - dot(b,t2)*t2;
 b = b / norm(b);
 
 n1 = cross(t1, b);
 n2 = cross(t2, b);
 d = dot(u,t1) * t2 + dot(u,n1) * n2 + ...
 dot(u,b) * b;
end
end

%% compute material directions
function [m1,m2] = computeMaterialDirectors(a1, a2, theta)
% Inputs:
% a1 is a matrix of size ne x 3. This contains the first
% reference director on each edge.
% a2 is a matrix of size ne x 3.
% theta is a vector of size ne ( theta = q(4:4:end) )
%
% Outputs:
% m1 and m2 are matrices of size ne x 3. Each column of
% these matrices contain the material directors on each
% edge.
ne = length(theta);
m1 = zeros(ne,3);
m2 = zeros(ne,3);
for c=1:ne % loop over edges
 cs = cos(theta(c));
 ss = sin(theta(c));
 m1(c,:) = cs * a1(c,:) + ss * a2(c,:);
 m2(c,:) = - ss * a1(c,:) + cs * a2(c,:);
end
end

%% get kappa
function kappa = getkappa( q, m1, m2 )
% Inputs:
% q: DOF vector
% m1 is the first material director (size is ne x 3)
% m2 is the second material director (size is ne x 3)
% Output
% kappa is a matrix of size nv x 2. Each node has kappa_1 and kappa_2.
nv = ( length(q) + 1) / 4; % length(q) = 4*nv-1
ne = nv - 1;
kappa = zeros(nv, 2);
for c=2:ne % All internal nodes except first and last
 node0 = q(4*c-7:4*c-5);
 node1 = q(4*c-3:4*c-1);
 node2 = q(4*c+1:4*c+3);
 
 m1e = m1(c-1, :); % m1 vector of c-1-th edge
 m2e = m2(c-1, :); % m2 vector of c-1-th edge
 m1f = m1(c, :); % m1 vector of c-th edge
 m2f = m2(c, :); % m2 vector of c-th edge
 
 kappaLocal = computekappa(node0, node1, node2, m1e, m2e, m1f, m2f );
 
 kappa(c,1) = kappaLocal(1);
 kappa(c,2) = kappaLocal(2);
end
end

%% objfun
function [q, u, a1, a2] = objfun(q0, u, a1, a2, freeIndex, ...
 tol, refTwist)
 
% Global
global Fg M dt
% Guess 
q = q0;
iter = 1;
error = 10 * tol;
while error > tol
 
 % Compute reference frame
 [a1Iterate, a2Iterate] = computeTimeParallel(a1, q0, q);
 
 % Compute reference twist
 tangent = computeTangent(q);
 refTwist_iterate = computeRefTwist(a1Iterate, tangent, refTwist);
 
 % Compute material frame
 theta = q(4:4:end);
 [m1, m2] = computeMaterialDirectors( a1Iterate, a2Iterate, theta);
 
 % Compute force and Jacobian
 [Fb, Jb] = getFb(q, m1, m2); % Bending
 [Ft, Jt] = getFt(q, refTwist_iterate); % Twisting
 [Fs, Js] = getFs(q); % Stretching
 % Equations of motion
 f = M/dt* ((q-q0)/dt - u) - Fb - Ft - Fs - Fg;
 J = M/dt^2 - Jb - Jt - Js;
 
 f_free = f(freeIndex);
 J_free = J(freeIndex, freeIndex);
 
 % Newton's update
 dq_free = J_free \ f_free ;
 q(freeIndex) = q(freeIndex) - dq_free;
 
 % Error
 error = sum( abs( f_free ) );
 
 fprintf('Iter = %d, error=%f\n', iter, error);
 iter = iter + 1;
 
end
u = (q - q0) / dt;
a1 = a1Iterate;
a2 = a2Iterate;
end

%% compute time parallel
function [a1, a2] = computeTimeParallel(a1_old, q0, q)
% a1_old is a matrix of size nex3: It contains the
% time-parallel reference director (a1) at each of
% the ne edges at the old time step
% q0 is the DOF vector at the old time step
% q is the DOF vector at the new time step
% Outputs:
% a1 is the time-parallel reference director (a1) at
% each of the ne edges at the new time step. Size of
% a1 is ne x 3.
% a2 is the second reference director. It's size is
% ne x 3.
nv = (length(q)+1)/4;
ne = nv - 1;
tangent0 = computeTangent(q0); % Tangents at old step
tangent = computeTangent(q); % Tangents at new step
a1 = zeros(ne, 3);
a2 = zeros(ne, 3);
for c=1:ne % loop over edges 
 t0 = tangent0(c,:); % tangent on c-th edge at old step
 t = tangent(c,:); % tangent on c-th edge at new step
 
 a1_local_old = a1_old(c,:); % a1 director on c-th edge
 % at old step
 a1_local = parallel_transport( a1_local_old, t0, t);
 % a1_local is the first reference director on c-th
 % edge at new step
 
 % Just to be careful: enforce a1 and t are perp.
 a1_local = a1_local - dot(a1_local, t) * t;
 a1_local = a1_local / norm(a1_local);
 
 a1(c,:) = a1_local; % store
 a2(c,:) = cross(t, a1_local);
end
end

%% compute ref twist
function refTwist = computeRefTwist(a1, tangent, refTwist)
% a1 is a matrix of size ne x 3. Each column of a1
% contains the first time-parallel reference director
% on each edge.
%
% tangent is a matrix of size ne x 3. It contains all
% the tangent vectors on the edges.
%
% refTwist is a vector of size nv (one reference twist
% for each node). This (optionally) is provided as an
% input (guess solution).
[ne, ~] = size(a1); % number of edges
nv = ne + 1;
% refTwist = zeros(nv, 1);
for c=2:ne % all internal nodes
 
 u0 = a1(c-1,:); % a1 vector of previous edge
 u1 = a1(c, :); % a1 vector of current edge
 t0 = tangent(c-1,:); % tangent of previous edge
 t1 = tangent(c,:); % tangent of current edge
 ut = parallel_transport(u0, t0, t1);
 % ut and u1 are the same?
 
 % Method 1: Okay? But 2*pi issue?
 % refTwist(c) = signedAngle(ut, u1, t1);
 
 % Method 2
 ut = rotateAxisAngle( ut, t1, refTwist(c) );
 refTwist(c) = refTwist(c) + signedAngle(ut, u1, t1);
end
end

%% get Fb
function [Fb, Jb] = getFb(q, m1, m2)
global kappaBar EI voronoiLength
nv = (length(q)+1) / 4;
Fb = zeros( size(q) );
Jb = zeros( length(q), length(q) );
for c=2:nv-1 % Compute bending force at each internel node
 
 node0 = transpose(q(4*c-7:4*c-5));
 node1 = transpose(q(4*c-3:4*c-1));
 node2 = transpose(q(4*c+1:4*c+3));
 
 m1e = m1(c-1, :); % m1 vector of c-1-th edge
 m2e = m2(c-1, :); % m2 vector of c-1-th edge
 m1f = m1(c, :); % m1 vector of c-th edge
 m2f = m2(c, :); % m2 vector of c-th edge
 
 [dF, dJ] = ...
 gradEb_hessEb(node0, node1, node2, ...
 m1e, m2e, m1f, m2f, ...
 kappaBar(c,:), voronoiLength(c), EI);
 
 ind = 4*c-7:4*c+3; % 11 numbers
 Fb(ind) = Fb(ind) - dF;
 Jb(ind, ind) = Jb(ind, ind) - dJ;
 
end
end

%% get Ft

function [Ft, Jt] = getFt(q, refTwist)
global GJ voronoiLength
nv = (length(q)+1) / 4;
Ft = zeros( size(q) );
Jt = zeros( length(q), length(q) );
for c=2:nv-1 % Compute bending force at each internel node
 
 node0 = transpose(q(4*c-7:4*c-5));
 node1 = transpose(q(4*c-3:4*c-1));
 node2 = transpose(q(4*c+1:4*c+3));
 theta_e = q(4*c-4);
 theta_f = q(4*c);
 
 [dF, dJ] = ...
 gradEt_hessEt(node0, node1, node2, ...
 theta_e, theta_f, refTwist(c), ...
 voronoiLength(c), GJ);
 
 ind = 4*c-7:4*c+3; % 11 numbers
 Ft(ind) = Ft(ind) - dF;
 Jt(ind, ind) = Jt(ind, ind) - dJ;
 
end
end

%% get Fs
function [Fs, Js] = getFs(q)
global EA refLen
nv = (length(q)+1) / 4;
ne = nv - 1;
Fs = zeros( size(q) );
Js = zeros( length(q), length(q) );
for c=1:ne % Each edge
 
 node1 = transpose(q(4*c-3:4*c-1)); % c-th node
 node2 = transpose(q(4*c+1:4*c+3)); % c+1-th node
 
 [dF, dJ] = ...
 gradEs_hessEs(node1, node2, refLen(c), EA);
 
 ind = [4*c-3, 4*c-2, 4*c-1, 4*c+1, 4*c+2,4*c+3]; % 6 numbers
 Fs(ind) = Fs(ind) - dF;
 Js(ind, ind) = Js(ind, ind) - dJ;
 
end
end

%% plotrod
function plotrod(q, a1, a2, m1, m2, ctime)
nv = (length(q)+1)/4;
x1 = q(1:4:end);
x2 = q(2:4:end);
x3 = q(3:4:end);
L = sum(sqrt( (x1(2:end) - x1(1:end-1)).^2 + ...
 (x2(2:end) - x2(1:end-1)).^2 + ...
 (x3(2:end) - x3(1:end-1)).^2)) / 3;
a1 = 0.1*L * a1;
a2 = 0.1*L * a2;
m1 = 0.1*L * m1;
m2 = 0.1*L * m2;
h1 = figure(1);
% set(h1, 'visible', 'off');
clf()
plot3(x1,x2,x3, 'ko-');
hold on
plot3(x1(1),x2(1),x3(1), 'r^');
for c=1:nv-1
 xa = q(4*c-3:4*c-1);
 xb = q(4*c+1:4*c+3);
 xp = (xa+xb)/2;
 p1 = plot3( [xp(1), xp(1) + a1(c,1)], [xp(2), xp(2) + a1(c,2)], ...
 [xp(3), xp(3) + a1(c,3)], 'b--', 'LineWidth', 2);
 p2 = plot3( [xp(1), xp(1) + a2(c,1)], [xp(2), xp(2) + a2(c,2)], ...
 [xp(3), xp(3) + a2(c,3)], 'c--', 'LineWidth', 2);
 p3 = plot3( [xp(1), xp(1) + m1(c,1)], [xp(2), xp(2) + m1(c,2)], ...
 [xp(3), xp(3) + m1(c,3)], 'r-');
 p4 = plot3( [xp(1), xp(1) + m2(c,1)], [xp(2), xp(2) + m2(c,2)], ...
 [xp(3), xp(3) + m2(c,3)], 'g-');
end
hold off
legend([p1,p2,p3,p4], 'a_1','a_2','m_1','m_2');
title(num2str(ctime, 't=%f'));
axis equal
view(3);
xlabel('x');
ylabel('y');
zlabel('z');
drawnow
end

%% corssMat
function A = crossMat(a)
A = [0, -a(3), a(2); ...
    a(3), 0, -a(1); ...
    -a(2), a(1), 0];
end

%% compute kappa

function kappa = computekappa(node0, node1, node2, m1e, m2e, m1f, m2f )
%
% Inputs:
% node0: 1x3 vector - position of the node prior to the "turning" node
% node1: 1x3 vector - position of the "turning" node
% node2: 1x3 vector - position of the node after the "turning" node
%
% m1e: 1x3 vector - material director 1 of the edge prior to turning
% m2e: 1x3 vector - material director 2 of the edge prior to turning
% m1f: 1x3 vector - material director 1 of the edge after turning
% m2f: 1x3 vector - material director 2 of the edge after turning
%
% Outputs:
% kappa: 1x2 vector - curvature at the turning node

t0 = (node1-node0) / norm(node1-node0);
t1 = (node2-node1) / norm(node2-node1);
kb = 2.0 * cross(t0, t1) / (1.0 + dot(t0, t1));

kappa = zeros(1, 2);
kappa1 = 0.5 * dot( kb, m2e + m2f);
kappa2 = -0.5 * dot( kb, m1e + m1f);
kappa(1) = kappa1;
kappa(2) = kappa2;

end

%% rotate axis angle
function vNew = rotateAxisAngle( v, z, theta )

% This function outputs a vector "vNew" after rotating the input vector
% "v" by an angle "theta" about a vector "z".
% Example: If v = [1;0;0] (x-axis), z = [0;0;1] (z-axis), and theta=pi/2,
% then vNew = [0;1;0] (y-axis.

if (theta == 0) 
    vNew = v;
else
    c = cos(theta);
    s = sin(theta);
    vNew = c*v + s*cross(z,v) + dot(z,v) * (1.0-c) * z;
end
end

%% signed angle
function angle = signedAngle( u, v, n )
% "angle" is signed angle from vector "u" to vector "v" with axis "n"
w = cross(u,v);
angle = atan2( norm(w), dot(u,v) );
if (dot(n,w) < 0) 
    angle = -angle;
end
end

%% gradEb_hessEb
function [dF, dJ] = ...
    gradEb_hessEb(node0, node1, node2, ...
    m1e, m2e, m1f, m2f, ...
    kappaBar, l_k, EI)
%
% Inputs:
% node0: 1x3 vector - position of the node prior to the "turning" node
% node1: 1x3 vector - position of the "turning" node
% node2: 1x3 vector - position of the node after the "turning" node
%
% m1e: 1x3 vector - material director 1 of the edge prior to turning
% m2e: 1x3 vector - material director 2 of the edge prior to turning
% m1f: 1x3 vector - material director 1 of the edge after turning
% m2f: 1x3 vector - material director 2 of the edge after turning
%
% kappaBar: 1x2 vector - natural curvature at the turning node
% l_k: voronoi length (undeformed) of the turning node
% EI: scalar - bending stiffness
%
% Outputs:
% dF: 11x1  vector - gradient of the bending energy at node1.
% dJ: 11x11 vector - hessian of the bending energy at node1.

%% Computation of gradient of the two curvatures
gradKappa = zeros(11,2);

ee = node1 - node0;
ef = node2 - node1;

norm_e = norm(ee);
norm_f = norm(ef);

te = ee / norm_e;
tf = ef / norm_f;

% Curvature binormal
kb = 2.0 * cross(te, tf) / (1.0 + dot(te, tf));

chi = 1.0 + dot(te, tf);
tilde_t = (te + tf) / chi;
tilde_d1 = (m1e + m1f) / chi;
tilde_d2 = (m2e + m2f) / chi;

% Curvatures
kappa1 = 0.5 * dot( kb, m2e + m2f); % CHECKED
kappa2 = -0.5 * dot( kb, m1e + m1f); % CHECKED

% kappa1 = kappa(c, 1);
% kappa2 = kappa(c, 2);

Dkappa1De = 1.0 / norm_e * (-kappa1 * tilde_t + cross(tf,tilde_d2));
Dkappa1Df = 1.0 / norm_f * (-kappa1 * tilde_t - cross(te,tilde_d2));

Dkappa2De = 1.0 / norm_e * (-kappa2 * tilde_t - cross(tf,tilde_d1));
Dkappa2Df = 1.0 / norm_f * (-kappa2 * tilde_t + cross(te,tilde_d1));

gradKappa(1:3, 1) = -Dkappa1De;
gradKappa(5:7, 1) = Dkappa1De - Dkappa1Df;
gradKappa(9:11, 1) = Dkappa1Df;

gradKappa(1:3, 2) = -Dkappa2De;
gradKappa(5:7, 2) = Dkappa2De - Dkappa2Df;
gradKappa(9:11, 2) = Dkappa2Df;

gradKappa(4, 1) = -0.5 * dot(kb, m1e);
gradKappa(8, 1) = -0.5 * dot(kb, m1f);
gradKappa(4, 2) = -0.5 * dot(kb, m2e);
gradKappa(8, 2) = -0.5 * dot(kb, m2f);

%% Computation of hessian of the two curvatures
DDkappa1 = zeros(11, 11);
DDkappa2 = zeros(11, 11);

norm2_e = norm_e^2;
norm2_f = norm_f^2;

tt_o_tt = tilde_t' * tilde_t; % must be 3x3. tilde_t is 1x3
tmp = cross(tf, tilde_d2);
tf_c_d2t_o_tt = tmp' * tilde_t; % must be 3x3
tt_o_tf_c_d2t = tf_c_d2t_o_tt'; % must be 3x3
kb_o_d2e = kb' * m2e; % must be 3x3
d2e_o_kb = kb_o_d2e'; % must be 3x3

Id3 = eye(3);
D2kappa1De2 ...
    = 1.0 / norm2_e * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt - tt_o_tf_c_d2t) ...
    - kappa1 / (chi * norm2_e) * (Id3 - te'*te) ...
    + 1.0 / (4.0 * norm2_e) * (kb_o_d2e + d2e_o_kb);

tmp = cross(te, tilde_d2);
te_c_d2t_o_tt = tmp' * tilde_t;
tt_o_te_c_d2t = te_c_d2t_o_tt';
kb_o_d2f = kb' * m2f;
d2f_o_kb = kb_o_d2f';

D2kappa1Df2 ...
    = 1.0 / norm2_f * (2 * kappa1 * tt_o_tt + te_c_d2t_o_tt + tt_o_te_c_d2t) ...
    - kappa1 / (chi * norm2_f) * (Id3 - tf'*tf) ...
    + 1.0 / (4.0 * norm2_f) * (kb_o_d2f + d2f_o_kb);

D2kappa1DeDf ...
    = -kappa1/(chi * norm_e * norm_f) * (Id3 + te'*tf) ...
    + 1.0 / (norm_e*norm_f) * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt + ...
    tt_o_te_c_d2t - crossMat(tilde_d2));
D2kappa1DfDe = D2kappa1DeDf';

tmp = cross(tf, tilde_d1);
tf_c_d1t_o_tt = tmp'*tilde_t; % must be 3x3
tt_o_tf_c_d1t = tf_c_d1t_o_tt'; % must be 3x3
kb_o_d1e = kb'*m1e; % must be 3x3
d1e_o_kb = kb_o_d1e'; % must be 3x3

D2kappa2De2 ...
    = 1.0 / norm2_e * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt + tt_o_tf_c_d1t) ...
    - kappa2 / (chi * norm2_e) * (Id3 - te'*te) ...
    - 1.0 / (4.0 * norm2_e) * (kb_o_d1e + d1e_o_kb);

tmp = cross(te, tilde_d1);
te_c_d1t_o_tt = tmp'*tilde_t; % must be 3x3
tt_o_te_c_d1t = te_c_d1t_o_tt'; % must be 3x3
kb_o_d1f = kb'*m1f; % must be 3x3
d1f_o_kb =  kb_o_d1f'; % must be 3x3

D2kappa2Df2 ...
    = 1.0 / norm2_f * (2 * kappa2 * tt_o_tt - te_c_d1t_o_tt - tt_o_te_c_d1t) ...
    - kappa2 / (chi * norm2_f) * (Id3 - tf'*tf) ...
    - 1.0 / (4.0 * norm2_f) * (kb_o_d1f + d1f_o_kb); % must be 3x3

D2kappa2DeDf ...
    = -kappa2/(chi * norm_e * norm_f) * (Id3 + te'*tf) ...
    + 1.0 / (norm_e*norm_f) * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt - tt_o_te_c_d1t + crossMat(tilde_d1));
% must be 3x3
D2kappa2DfDe = D2kappa2DeDf'; % must be 3x3

D2kappa1Dthetae2 = -0.5 * dot(kb, m2e);
D2kappa1Dthetaf2 = -0.5 * dot(kb, m2f);
D2kappa2Dthetae2 =  0.5 * dot(kb, m1e);
D2kappa2Dthetaf2 =  0.5 * dot(kb, m1f);

D2kappa1DeDthetae ...
    = 1.0 / norm_e * (0.5 * dot(kb, m1e) * tilde_t - 1.0 / chi * cross(tf, m1e));
D2kappa1DeDthetaf ...
    = 1.0 / norm_e * (0.5 * dot(kb, m1f) * tilde_t - 1.0 / chi * cross(tf, m1f));
D2kappa1DfDthetae ...
    = 1.0 / norm_f * (0.5 * dot(kb, m1e) * tilde_t + 1.0 / chi * cross(te, m1e));
D2kappa1DfDthetaf ...
    = 1.0 / norm_f * (0.5 * dot(kb, m1f) * tilde_t + 1.0 / chi * cross(te, m1f));
D2kappa2DeDthetae ...
    = 1.0 / norm_e * (0.5 * dot(kb, m2e) * tilde_t - 1.0 / chi * cross(tf, m2e));
D2kappa2DeDthetaf ...
    = 1.0 / norm_e * (0.5 * dot(kb, m2f) * tilde_t - 1.0 / chi * cross(tf, m2f));
D2kappa2DfDthetae ...
    = 1.0 / norm_f * (0.5 * dot(kb, m2e) * tilde_t + 1.0 / chi * cross(te, m2e));
D2kappa2DfDthetaf ...
    = 1.0 / norm_f * (0.5 * dot(kb, m2f) * tilde_t + 1.0 / chi * cross(te, m2f));

% Curvature terms
DDkappa1(1:3, 1:3)  =   D2kappa1De2;
DDkappa1(1:3, 5:7)  = - D2kappa1De2 + D2kappa1DeDf;
DDkappa1(1:3, 9:11) =               - D2kappa1DeDf;
DDkappa1(5:7, 1:3)  = - D2kappa1De2                + D2kappa1DfDe;
DDkappa1(5:7, 5:7)  =   D2kappa1De2 - D2kappa1DeDf - D2kappa1DfDe + D2kappa1Df2;
DDkappa1(5:7, 9:11) =                 D2kappa1DeDf                - D2kappa1Df2;
DDkappa1(9:11, 1:3)  =                              - D2kappa1DfDe;
DDkappa1(9:11, 5:7)  =                                D2kappa1DfDe - D2kappa1Df2;
DDkappa1(9:11, 9:11) =                                               D2kappa1Df2;

% Twist terms
DDkappa1(4, 4)     =   D2kappa1Dthetae2;
DDkappa1(8, 8)     =   D2kappa1Dthetaf2;

% Curvature-twist coupled terms
DDkappa1(1:3, 4)   = - D2kappa1DeDthetae;
DDkappa1(5:7, 4)   =   D2kappa1DeDthetae - D2kappa1DfDthetae;
DDkappa1(9:11,4)   =                       D2kappa1DfDthetae;
DDkappa1(4, 1:3)   =   transpose(DDkappa1(1:3, 4));
DDkappa1(4, 5:7)   =   transpose(DDkappa1(5:7, 4));
DDkappa1(4, 9:11)  =   transpose(DDkappa1(9:11,4));

% Curvature-twist coupled terms
DDkappa1(1:3, 8)   = - D2kappa1DeDthetaf;
DDkappa1(5:7, 8)   =   D2kappa1DeDthetaf - D2kappa1DfDthetaf;
DDkappa1(9:11, 8)  =                       D2kappa1DfDthetaf;
DDkappa1(8, 1:3)   =   transpose(DDkappa1(1:3, 8));
DDkappa1(8, 5:7)   =   transpose(DDkappa1(5:7, 8));
DDkappa1(8, 9:11)  =   transpose(DDkappa1(9:11,8));

% Curvature terms
DDkappa2(1:3, 1:3) =   D2kappa2De2;
DDkappa2(1:3, 5:7) = - D2kappa2De2 + D2kappa2DeDf;
DDkappa2(1:3, 9:11) =               - D2kappa2DeDf;
DDkappa2(5:7, 1:3) = - D2kappa2De2                + D2kappa2DfDe;
DDkappa2(5:7, 5:7) =   D2kappa2De2 - D2kappa2DeDf - D2kappa2DfDe + D2kappa2Df2;
DDkappa2(5:7, 9:11)=                 D2kappa2DeDf                - D2kappa2Df2;
DDkappa2(9:11, 1:3)=                              - D2kappa2DfDe;
DDkappa2(9:11, 5:7)=                                D2kappa2DfDe - D2kappa2Df2;
DDkappa2(9:11, 9:11)=                                               D2kappa2Df2;

% Twist terms
DDkappa2(4, 4)     = D2kappa2Dthetae2;
DDkappa2(8, 8)     = D2kappa2Dthetaf2;

% Curvature-twist coupled terms
DDkappa2(1:3, 4)   = - D2kappa2DeDthetae;
DDkappa2(5:7, 4)   =   D2kappa2DeDthetae - D2kappa2DfDthetae;
DDkappa2(9:11,4)   =                       D2kappa2DfDthetae;
DDkappa2(4, 1:3)   =   transpose(DDkappa2(1:3, 4));
DDkappa2(4, 5:7)   =   transpose(DDkappa2(5:7, 4));
DDkappa2(4, 9:11)  =   transpose(DDkappa2(9:11,4));

% Curvature-twist coupled terms
DDkappa2(1:3, 8)   = - D2kappa2DeDthetaf;
DDkappa2(5:7, 8)   =   D2kappa2DeDthetaf - D2kappa2DfDthetaf;
DDkappa2(9:11,8)   =                       D2kappa2DfDthetaf;
DDkappa2(8, 1:3)   =   transpose(DDkappa2(1:3, 8));
DDkappa2(8, 5:7)   =   transpose(DDkappa2(5:7, 8));
DDkappa2(8,9:11)   =   transpose(DDkappa2(9:11,8));
    
%% Gradient of Eb
EIMat = [ EI 0; ...
    0 EI];
kappaVector = [kappa1 kappa2];
dkappaVector = kappaVector - kappaBar;
dF = gradKappa * EIMat * dkappaVector' / l_k;

%% Hessian of Eb
dJ = 1.0 / l_k * gradKappa * EIMat * transpose(gradKappa);
temp = 1.0 / l_k * dkappaVector * EIMat;
dJ = dJ + temp(1) * DDkappa1 + temp(2) * DDkappa2;

end

%% gradEt_hessEt
function [dF, dJ] = ...
    gradEt_hessEt(node0, node1, node2, ...
    theta_e, theta_f, refTwist, ...
    l_k, GJ)
%
% Inputs:
% node0: 1x3 vector - position of the node prior to the "twisting" node
% node1: 1x3 vector - position of the "twisting" node
% node2: 1x3 vector - position of the node after the "twisting" node
%
% theta_e: scalar - twist angle of the first edge
% theta_f: scalar - twist angle of the second (last) edge
%
% l_k: voronoi length (undeformed) of the turning node
% GJ: scalar - twisting stiffness
%
% Outputs:
% dF: 11x1  vector - gradient of the twisting energy at node1.
% dJ: 11x11 vector - hessian of the twisting energy at node1.

%% Computation of gradient of the twist
gradTwist = zeros(11,1);

ee = node1 - node0;
ef = node2 - node1;

norm_e = norm(ee);
norm_f = norm(ef);

norm2_e = norm_e^2;
norm2_f = norm_f^2;

te = ee / norm_e;
tf = ef / norm_f;

% Curvature binormal
kb = 2.0 * cross(te, tf) / (1.0 + dot(te, tf));

gradTwist(1:3) = -0.5 / norm_e * kb;
gradTwist(9:11) = 0.5 / norm_f * kb;
gradTwist(5:7) = -(gradTwist(1:3)+gradTwist(9:11));
gradTwist(4) = -1;
gradTwist(8) = 1;

%% Computation of hessian of twist
DDtwist = zeros(11, 11);

chi = 1.0 + dot(te, tf);
tilde_t = (te + tf) / chi;

D2mDe2 = -0.25 / norm2_e * ( kb' * (te + tilde_t) ...
    + (te + tilde_t)' * kb);
D2mDf2 = -0.25 / norm2_f  * ( kb' * (tf + tilde_t) ...
    + (tf + tilde_t)' * kb );
D2mDeDf = 0.5 / ( norm_e * norm_f ) * ( 2.0 / chi * crossMat( te ) - ...
    kb' * tilde_t );
D2mDfDe = D2mDeDf';

DDtwist(1:3,1:3) = D2mDe2;
DDtwist(1:3, 5:7) = -D2mDe2 + D2mDeDf;
DDtwist(5:7, 1:3) = -D2mDe2 + D2mDfDe;
DDtwist(5:7, 5:7) = D2mDe2 - ( D2mDeDf + D2mDfDe ) + D2mDf2;
DDtwist(1:3, 9:11) = -D2mDeDf;
DDtwist(9:11, 1:3) = -D2mDfDe;
DDtwist(9:11, 5:7) = D2mDfDe - D2mDf2;
DDtwist(5:7,9:11) = D2mDeDf - D2mDf2;
DDtwist(9:11,9:11) = D2mDf2;
    
%% Gradient of Et
integratedTwist = theta_f - theta_e + refTwist;
dF = GJ/l_k * integratedTwist * gradTwist;

%% Hessian of Eb
dJ = GJ/l_k * (integratedTwist * DDtwist + gradTwist*gradTwist');

end

%% gradEs_HessEs
function [dF, dJ] = ...
    gradEs_hessEs(node0, node1, ...
    l_k, EA)
%
% Inputs:
% node0: 1x3 vector - position of the first node
% node1: 1x3 vector - position of the last node
%
% l_k: reference length (undeformed) of the edge
% EA: scalar - stretching stiffness - Young's modulus times area
%
% Outputs:
% dF: 6x1  vector - gradient of the stretching energy between node0 and node 1.
% dJ: 6x6 vector - hessian of the stretching energy between node0 and node 1.

%% Gradient of Es
edge = (node1 - node0)'; % 3x1 edge vector
edgeLen = norm(edge);
tangent = edge / edgeLen;
epsX = edgeLen/l_k - 1;
dF_unit = EA * tangent * epsX;

dF = zeros(6,1);
dF(1:3) = - dF_unit;
dF(4:6) =   dF_unit;

%% Hessian of Es
Id3 = eye(3);
M = EA * ( ...
    (1/l_k - 1/edgeLen) * Id3 + ...
    1/edgeLen * (edge*edge')/ edgeLen^2 ...
    ); %Note edge * edge' must be 3x3

dJ = zeros(6,6);
dJ(1:3, 1:3) = M;
dJ(4:6, 4:6) = M;
dJ(1:3, 4:6) = - M;
dJ(4:6, 1:3) = - M;

end

