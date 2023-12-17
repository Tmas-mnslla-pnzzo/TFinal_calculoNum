%constantes:
L=1;

#funciones

function u_1=u1(x,y)
  u_1=x*(1-x)*y*(1-y);
endfunction

function u_2=u2(x,y)
  u_2=((x-0.5)**2+(y-0.5)**2+(1/36))**2;
endfunction

function f_1=f1(x,y)
  f_1=2*x*y*y-2*x*y-y*y+y;
endfunction

function f_2=f2(x,y)
  f_2=2*y*x*x-2*y*x-x*x+x;
endfunction

function f_3=f3(x,y)
  f_3=4*(x**3)-6*(x**2)+4*x*(y**2)-4*x*y+(37/9)*x-2*(y**2)+2*y-(19/18);
endfunction

function f_4=f4(x,y)
  f_4=4*(y**3)-6*(y**2)+4*y*(x**2)-4*y*x+(37/9)*y-2*(x**2)+2*x-(19/18);
endfunction

function f_5=f5(x,y)
  f_5=2*(x*x+y*y-x-y);
endfunction

function f_6=f6(x,y)
  f_6=16*(x*x+y*y-x-y+(2*37)/(16*9));
endfunction

function f_0=f(x,y)
  f_0=u1(x,y)*f6(x,y)+u2(x,y)*f5(x,y)+2*(f1(x,y)*f3(x,y)+f2(x,y)*f4(x,y));
endfunction

function mT=matriz_T(ni)
  N=ni+2;
  a1=1;
  a2=-4;
  a3=1;
  M0=[-4 1 zeros(1,N-2)];
  MN=[zeros(1,N-2)  1 -4];
  M1=diag(a1*ones(1,N-1),1)+diag(a2*ones(1,N))  + diag(a3*ones(1,N-1),-1);
  M1=[M0;M1];
  M1(N+1,:)=MN;
  M1(2,:)=[];
  mT=M1;
endfunction 

function mi=matriz_i(ni)
  N=ni+2;
  a1=0;
  a2=1;
  a3=0;
  M1=diag(a1*ones(1,N-1),1)+diag(a2*ones(1,N))  + diag(a3*ones(1,N-1),-1);
  mi=M1;
endfunction 

function z=concat(m1,m2,m3,Ni)
  M7=diag(1*ones(1,Ni-1),1)+diag(2*ones(1,Ni))  + diag(1*ones(1,Ni-1),-1);
  L_T=[];
  for u=M7
    L_i=[];
    for v=1:length(u)
      if u(v)==2
        L_i=[L_i, m1];
      else
        if u(v)==1
          L_i=[L_i, m2];
        else
          L_i=[L_i, m3];
        endif
      endif
    end
    L_T=[L_T; L_i];
  end
  z=L_T;
endfunction 

#Parametros
N_n=40;

D_c=ones(N_n,N_n);
i0=round(N_n/2);
j0=round(N_n/2);
R=round(N_n/6);

[x_c,y_c]=meshgrid(1:N_n);
D_c((x_c-i0).^2+(y_c-j0).^2<R^2)=0;

Num=N_n-1;
N_v=Num-1;
N=N_v+2;
x = 0:L/(N-1):L;
Dx=L/(N-1);
x(end) = [];
x(1) = [];

N=N_v;

#Desarrollo

F=[];
U=[];
for k=x
  F_i=[];
  U_i=[];
  for j=x
    F_i=[F_i f(j,k)];
    U_i=[U_i u1(j,k)*u2(j,k)];
  end
  F=[F;F_i];
  U=[U;U_i];
end


ceros=zeros(1,N)';
ceros2=zeros(1,N+2);

L2=F'(:)*(Dx**2);

M_T=matriz_T(N-2);
M_I=matriz_i(N-2);
M_0=diag(zeros(1,N-1),1);
Ni=N;

M_T=concat(M_T,M_I,M_0,Ni);

M_T=inv(M_T);


Z=(M_T*L2);

Z=reshape (Z, N, N);
Z=horzcat(ceros,Z);
Z=horzcat(Z,ceros);
Z=vertcat(Z,ceros2);
Z=vertcat(ceros2,Z);

x_DF= 0:L/(N+1):L;
Z_DF=Z'(:);
D_c=D_c'(:);
for jg=1:1:(N_n**2)
  if D_c(jg)==0
    Z_DF(jg)=NaN;
  endif
end

Z_DF=reshape (Z_DF, N_n, N_n);

figure 3
surf(x_DF,x_DF,Z_DF); xlabel ("Distancia [m]"); ylabel ("Distancia [m]"); zlabel("Desplazamiento [m]");
view([45 45 130])
