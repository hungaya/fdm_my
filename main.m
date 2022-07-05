clc
clear all
close all

Nx=10;
Ny=10;
Nt=10;
Lx=1;
Ly=1;
T=1;


x = linspace(0,Lx,Nx);   %Mesh
y = linspace(0,Ly,Ny);
t = linspace(0,T, Nt);


hx=x(2)-x(1);     %Mesh size
hy=y(2)-y(1);
dt=t(2)-t(1);

N=Nx*Ny;       %number of unknown at one step
M=zeros(N,N);
B=zeros(N,1);  % N rows, N columns
U_0 = zeros(N,1);   % N rows, 1 columns


syms sx sy real;

hamd = 1;
dx = diff(hamd, sx);
dy = diff(hamd, sy);

hamd = @(vx, vy) subs(hamd, {sx,sy}, {vx, vy});
dx = @(vx,vy) subs(dx, {sx,sy}, {vx, vy});
dy = @(vx,vy) subs(dy, {sx,sy}, {vx, vy});
dx(1,2);

for i=1:Nx
    for j=1:Ny
         n=i+(j-1)*Nx;
         U_0(n,1)= u_0(x(i),y(j));
    end
end



for k=1:2 % loop over t-direction
    %Interior points
    for i=2:Nx-1 % loop over x-direction skipping 1  st and last grid points.
        for j=2:Ny-1
            n=i+(j-1)*Nx;
            M(n,n)= (2*hamd(x(i),y(j))/(hx^2)+ 2*hamd(x(i),y(j))/(hy^2)-1/dt);
            M(n,n-1)= -hamd(x(i),y(j))/(hx^2);
            M(n,n+1)= -hamd(x(i),y(j))/(hx^2)-dx(x(i),y(j))/hx; 
            M(n,n-Nx)= -hamd(x(i),y(j))/(hy^2);
            M(n,n+Nx)= -hamd(x(i),y(j))/(hy^2)-dy(x(i),y(j))/hx;
            B(n,1)   = f(U_0(n,1));
        end
    end

%Boundary conditions
    % BC:{u^n}_1j={u^n}_2j
    i=1;
    for j=1:Ny
        n=i+(j-1)*Nx;
        M(n,n)= 1;
        B(n,1)=U_0(2+(j-1)*Nx,1);
    end
    % BC:{u^n}_Nxj={u^n}_(Nx-1)j
    i=Nx;
    for j=1:Ny
        n=i+(j-1)*Nx;
        M(n,n)=1;
         B(n,1)=U_0(Nx-1+(j-1)*Nx,1);
        
    end
    % BC:{u^n}_i1={u^n}_i2
    j=1;
    for i=1:Nx
        n=i+(j-1)*Nx;
        M(n,n)=1;
        B(n,1)=U_0(i+(2-1)*Nx,1); 
    end
    % BC:{u^n}_i1={u^n}_i2
    j=Ny;
    for i=1:Nx
        n=i+(j-1)*Nx;
        M(n,n)=1;
        B(n,1)=U_0(i+(Ny-1-1)*Nx,1);
     
    end      
    
    
%calculate U_1
    disp(M);
    U_1_vec=dt *(B - M*(U_0));
    
  
   for i=1:Nx
    for j=1:Ny
        n=i+(j-1)*Nx;
        U_1(i,j)=U_1_vec(n);
    end
   end


surf(x,y,U_1');
xlabel('x')
ylabel('y')

%reset the value of U_0   
  for i=1:Nx
        for j=1:Ny
            n=i+(j-1)*Nx;
            U_0(n,1)=U_1_vec(n,1);
        end
  end
end
