clc
clear 
close all
tic

Nx = 50;
Ny = 50;
Nt = 50;

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
dx = diff(hamd, sx);  % case 1: d is continuous
dy = diff(hamd, sy);  % case 1: d is continuous



hamd = @(vx, vy) subs(hamd, {sx, sy}, {vx, vy});
dx = @(vx, vy) subs(dx, {sx, sy}, {vx, vy});  % case 1: d is continuous
dy = @(vx, vy) subs(dy, {sx, sy}, {vx, vy});  % case 1: d is continuous


X = repmat(x', 1, Ny);
Y = repmat(y, Nx, 1);

D = double(hamd(X, Y)); 
% Dx = double((hamd(X+hx, Y)-hamd(X, Y)/hx)); % case 2: d is uncontinuous
% Dy = double((hamd(X+hy, Y)-hamd(X, Y)/hy));   % case 2: d is uncontinuous

% 
Dx = double(dx(X, Y));  % case 1: d is continuous
Dy = double(dy(X, Y));  % case 1: d is continuous

D_hx2 = D ./ (hx^2);
D_hy2 = D ./ (hy^2);
Dx_hx = Dx ./ hx;
Dy_hy = Dy ./ hy;

U_0 = reshape(u_0(X, Y), [], 1);

for i = 2:(Nx-1)
    for j = 2:(Ny-1)
        % convert the cell (i,j) to the nth grid point
        n = i + (j-1)*Nx;

        M(n, n)    = 2*D_hx2(i,j) + 2*D_hy2(i,j) + Dx_hx(i,j) + Dy_hy(i,j) + (1/dt);
        M(n, n-1)  = -D_hx2(i,j);
        M(n, n+1)  = -D_hx2(i,j) - Dx_hx(i,j);
        M(n, n-Nx) = -D_hy2(i,j);
        M(n, n+Nx) = -D_hy2(i,j) - Dy_hy(i,j);
    end
end
% 
% % loop over t-direction
for k = 1:Nt
    % interior points
    for i = 2:(Nx-1)
        for j = 2:(Ny-1)
            B(n) = f(U_0(n)) + U_0(n)/dt;
        end
    end
    
   %Boundary conditions
    % BC:{u^n+1}_1j={u^n+1}_2j
    i=1;
    for j=1:Ny-1
        n=i+(j-1)*Nx;
        M(n,n)= 1;
        M(n, n+1)=-1;
        B(n)=0;
    end
    % BC:{u^n+1}_Nxj={u^n+1}_(Nx-1)j
    i=Nx;
    for j=2:Ny-1
        n=i+(j-1)*Nx;
         M(n,n)= 1;
        M(n, n-1)=-1;
        B(n)=0;
        
    end
    % BC:{u^n+1}_i1={u^n+1}_i2
    j=1;
    for i=2:Nx
        n=i+(j-1)*Nx;
         M(n,n)= 1;
        M(n, n+Nx)=-1;
        B(n)=0;
    end
    % BC:{u^n+1}_iNy={u^n+1}_iNy-1
    j=Ny;
    for i=1:Nx
        n=i+(j-1)*Nx;
        M(n,n)= 1;
        M(n, n-Nx)=-1;
        B(n)=0;
     
    end 
    
    U_1_vec= inv(M)*B;
    U_1 = reshape(U_1_vec, Nx, Ny);

    % reset the value of U_0   
    U_0 = U_1_vec;
    disp(sprintf('%d step', k));
 end
 toc

surf(x,y,U_1');
xlabel('x');
ylabel('y');
