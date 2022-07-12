clc
clear all
close all
tic

Nx = 10;
Ny = 10;
Nt = 10;

Lx = 1;
Ly = 1;
T = 1;

% Mesh
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
t = linspace(0,  T, Nt);

% Mesh size
hx = x(2) - x(1);
hy = y(2) - y(1);
dt = t(2) - t(1);

% number of unknown at one step
N = Nx * Ny;

M = zeros(N, N);
B = zeros(N, 1);   % N rows, N columns
U_0 = zeros(N, 1); % N rows, 1 column
D = zeros(Nx, Ny);

syms sx sy real;
hamd = 1;
dx = diff(hamd, sx);
dy = diff(hamd, sy);

hamd = @(vx, vy) subs(hamd, {sx, sy}, {vx, vy});
dx = @(vx, vy) subs(dx, {sx, sy}, {vx, vy});
dy = @(vx, vy) subs(dy, {sx, sy}, {vx, vy});

X = repmat(x', 1, Ny);
Y = repmat(y, Nx, 1);

D = double(hamd(X, Y)); 
Dx = double(dx(X, Y));
Dy = double(dy(X, Y));

D_hx2 = D ./ (hx^2);
D_hy2 = D ./ (hy^2);
Dx_hx = Dx ./ hx;
Dy_hy = Dy ./ hy;

U_0 = reshape(u_0(X, Y), [], 1);

% loop over t-direction
for k = 1:Nt
    % interior points
    for i = 2:(Nx-1)
        for j = 2:(Ny-1)
            % convert the cell (i,j) to the nth grid point
            n = i + (j-1)*Nx;
           
            M(n, n)    = 2*D_hx2(i,j) + 2*D_hy2(i,j) + Dx_hx(i,j) + Dy_hy(i,j) + (1/dt);
            M(n, n-1)  = -D_hx2(i,j);
            M(n, n+1)  = -D_hx2(i,j) - Dx_hx(i,j);
            M(n, n-Nx) = -D_hy2(i,j);
            M(n, n+Nx) = -D_hy2(i,j) - Dy_hy(i,j);
           
            B(n) = f(U_0(n)) + U_0(n)/dt;
        end
    end
