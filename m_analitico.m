L=1;

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
N_n=50;
Num=N_n-1;
N_v=Num-1;
N=N_v+2;
x = 0:L/(N-1):L;
Dx=L/(N-1);
x(end) = [];
x(1) = [];

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

ceros=zeros(1,N_v)';
ceros2=zeros(1,N_n);


Z_A=horzcat(U',ceros); Z_A=horzcat(ceros,Z_A);
Z_A=vertcat(Z_A,ceros2); Z_A=vertcat(ceros2,Z_A);
x_A=[0 x 1];

figure 2
waterfall(x_A,x_A,Z_A)
xlabel ("Distancia [m]"); ylabel ("Distancia [m]"); zlabel("Desplazamiento [m]");
view([45 45 130])
