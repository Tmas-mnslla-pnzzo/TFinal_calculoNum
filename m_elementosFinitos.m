pkg load io

%parametros
L=1;
Cn=8;
r=8;
R=1/6;
nodos_x=3*r-2;

%funciones

u1 = @(x,y) x.*(1-x).*y.*(1-y);
u2 = @(x,y) ((x-0.5).^2+(y-0.5).^2+(1/36)).^2;
f1  = @(x,y) x.*y.*y.*2-x.*y*2-y.*y+y;
f2  = @(x,y) y.*x.*x.*2-y.*x.*2-x.*x+x;
f3  = @(x,y) (x.^3).*4-(x.^2).*6+x.*(y.^2).*4-x.*y.*4+(37/9).*x-(y.^2).*2+y.*2-(19/18);
f4  = @(x,y) (y.^3).*4-(y.^2).*6+y.*(x.^2).*4-y.*x.*4+(37/9).*y-(x.^2).*2+x.*2-(19/18);
f5  = @(x,y) (x.*x+y.*y-x-y).*2;
f6  = @(x,y) (x.*x+y.*y-x-y+(2*37)/(16*9)).*16;
f    = @(x,y) u1(x,y).*f6(x,y)+u2(x,y).*f5(x,y)+2.*(f1(x,y).*f3(x,y)+f2(x,y).*f4(x,y));

function A=area(p1,p2,p3)
  a1=p1(1);
  a2=p1(2);
  b1=p2(1);
  b2=p2(2);
  c1=p3(1);
  c2=p3(2);
  a=sqrt((a1-c1)**2+(c2-a2)**2);
  b=sqrt((b1-a1)**2+(a2-b2)**2);
  c=sqrt((b1-c1)**2+(c2-b2)**2);
  s=(a+b+c)/2;
  A=s*(s-a)*(s-b)*(s-c);
  A=sqrt(A);
endfunction

function A=area2(a1,a2,b1,b2,c1,c2)
  a=sqrt((a1-c1)**2+(c2-a2)**2);
  b=sqrt((b1-a1)**2+(a2-b2)**2);
  c=sqrt((b1-c1)**2+(c2-b2)**2);
  s=(a+b+c)/2;
  A=s*(s-a)*(s-b)*(s-c);
  A=sqrt(A);
endfunction

function inM=M_i(p1,p2,p3,k)
  a1=p1(1);
  a2=p1(2);
  b1=p2(1);
  b2=p2(2);
  c1=p3(1);
  c2=p3(2);
  inM=[a1,a2,1;b1,b2,1;c1,c2,1];
  inM=inv(inM);
  if k==1
   inM=inM*[1 0 0]';
  elseif k==2
    inM=inM*[0 1 0]';
  elseif k==3
    inM=inM*[0 0 1]';
  endif
endfunction

function grad=grad_nodo(p1,p2,p3)
  area=area(p1,p2,p3);
  
  C1=M_i(p1,p2,p3,1);
  c_i=[C1(1),C1(2)];
  
  C2=M_i(p1,p2,p3,2);
  c_j=[C2(1),C2(2)];
 
  C3=M_i(p1,p2,p3,3);
  c_k=[C3(1),C3(2)];
 
  grad=[c_i*c_i',c_i*c_j',c_i*c_k'];
  grad=[grad; c_i*c_j',c_j*c_j',c_j*c_k'];
  grad=[grad; c_i*c_k',c_j*c_k',c_k*c_k'];
  grad=area*grad;
endfunction

function Ub=ubicar(n1,n2,n3,S,M)
  M(n1,n1)=M(n1,n1)+S(1,1);
  M(n1,n2)=M(n1,n2)+S(1,2);
  M(n1,n3)=M(n1,n3)+S(1,3);
  M(n2,n2)=M(n2,n2)+S(2,2);
  M(n2,n1)=M(n2,n1)+S(2,1);
  M(n2,n3)=M(n2,n3)+S(2,3);
  M(n3,n3)=M(n3,n3)+S(3,3);
  M(n3,n1)=M(n3,n1)+S(3,1);
  M(n3,n2)=M(n3,n2)+S(3,2);
  Ub=M;
endfunction

function w=omega(x,y,a,b,c)
  w=a*x+b*y+c;
endfunction

%generar circulo
C=[];
for ci=0:1:(Cn-2)
  C= [C; R*cos((ci+1)*(90/Cn)*(pi/180)) R*sin((ci+1)*(90/Cn)*(pi/180))];
end
C=[[1/6 0];C];
C=[C;[0 1/6]];
Xr4=C(:,1)+0.5;
Yr4=C(:,2)+0.5;

%cuadrados
S=0:(R)/(r-1):R;
S(1)=[];
Sx=(R)*ones(1,size(S)(2));
Dx=R:(R)/(r-1):0.5;

[Xr1,Yr1]=meshgrid(Dx,[0 S]);
[Xr2,Yr2]=meshgrid([0 S],Dx);
[Xr3,Yr3]=meshgrid(Dx,Dx);
[Xr5,Yr5]=meshgrid([0 S],[0 S]);

Yr1=Yr1+0.5; Xr1=Xr1+0.5;
Yr2=Yr2+0.5; Xr2=Xr2+0.5;
Yr3=Yr3+0.5; Xr3=Xr3+0.5;
Yr5=Yr5+0.5; Xr5=Xr5+0.5;

Xr=horzcat(Xr2,Xr3); Xrr=horzcat(Xr5,Xr1); Xr=vertcat(Xr,Xrr);

Yr=vertcat(Yr2,Yr5); Yrr=vertcat(Yr3,Yr1); Yr=horzcat(Yr,Yrr);


%generacion de nodos

Px=[Xr1(:) ;Xr2(:) ;Xr3(:); Xr4(:)];
Py=[Yr1(:) ;Yr2(:) ;Yr3(:) ;Yr4(:)];
P=[Px,Py];P=unique(P(:,1:2),'rows');
TF_1=delaunay(P(:,1),P(:,2));

%eliminacion trianglos
idx = find(ismember(P, C+0.5,'rows'));
nk=nchoosek(idx,3);
G=ismember(TF_1,nk);
indexes=find(ismember(G,[1,1,1],'rows'));
TF_1(indexes(:),:)=[];

%triplot(TF_1,Px,Py)

%triangulacion provicional
TF_2=[];
tf_ii=1;
for tfcc=1:1:(nodos_x-1)
  TF2_i=[];
  for tfc=1:1:(nodos_x-1)
    TF2_i=[TF2_i; tf_ii tf_ii+1 nodos_x+tf_ii+1];
    TF2_i=[TF2_i; tf_ii nodos_x+tf_ii nodos_x+tf_ii+1];
    tf_ii=tf_ii+1;
  end
  TF_2=[TF_2; TF2_i];
  tf_ii=tf_ii+1;
end

h=(1-0.5)/(nodos_x-1);
Xp=0.5:h:1;
[Xp,Yp]=meshgrid(Xp,Xp);
Pp=horzcat(Yp(:),Xp(:));

%desarrollo del metodo
TF=TF_1;

M=zeros(size(P)(1),size(P)(1));

for indez=1:1:size(TF)(1)
  n_1=TF(indez,:)(1);
  n_2=TF(indez,:)(2);
  n_3=TF(indez,:)(3);
  p_1=[P(n_1,:)];
  p_2=[P(n_2,:)];
  p_3=[P(n_3,:)];
  S=grad_nodo(p_1,p_2,p_3);
  M=ubicar(n_1,n_2,n_3,S,M);
end

B=[];
for s_i=1:1:size(P)(1)
  [file,row]=find(TF==s_i);
  suma=0;
  for t_i=1:1:size(TF(file,:))(1)
    
    P_I=TF(file,:)(t_i,:);
    
    if P_I(1)==s_i
      p11=P(P_I(1),:);
      p22=P(P_I(2),:);
      p33=P(P_I(3),:);
    else
      if P_I(2)==s_i
        p11=P(P_I(2),:);
        p22=P(P_I(1),:);
        p33=P(P_I(3),:);
      else
        p11=P(P_I(3),:);
        p22=P(P_I(1),:);
        p33=P(P_I(2),:);
      endif
    endif
    
    a=M_i(p11,p22,p33,1)(1);
    b=M_i(p11,p22,p33,1)(2);
    c=M_i(p11,p22,p33,1)(3);
    
    a1=p11(1);
    a2=p11(2);
    b1=p22(1);
    b2=p22(2);
    c1=p33(1);
    c2=p33(2);
    
    J=[(b1-a1) (c1-a1); (b2-a2) (c2-a2)];
    J2=inv(J);
   
    fun = @(u,v) -f(((b1-a1).*u+(c1-a1).*v+a1),((b2-a2).*u+(c2-a2).*v+a2)).*(a.*((b1-a1).*u+(c1-a1).*v+a1)+b.*((b2-a2).*u+(c2-a2).*v+a2)+c);
    ymx = @(u) 1 - u;
    inte1 = integral2(fun,0,1,0,ymx);
    suma=suma+abs(det(J))*inte1;
  endfor
  B=[B;suma];
end

[FX,RX]=find(P==1);
B(unique(FX))=0;

for i=1:1:size(unique(FX))(1)
  M(unique(FX)(i),:)=0;
  M(:,unique(FX)(i))=0;
  M(unique(FX)(i),unique(FX)(i))=1;
end

K=inv(M);
U=K*B;
Z_EF=horzcat(P,U);

mx_y=Z_EF((find(Z_EF(:,3)==max(Z_EF(:,3)))),:)(1);
ind_mx_y=find(Z_EF(:,2)==mx_y);
L1_mx_y=Z_EF(ind_mx_y,:)(:,1);
L1_mx_y2=L1_mx_y-0.5; L1_mx_y2(end)=[];
L1_mx_y=vertcat(L1_mx_y2,L1_mx_y);

L2_mx_y=Z_EF(ind_mx_y,:)(:,3);
L2_mx_y2=fliplr(L2_mx_y'); L2_mx_y2(end)=[];
L2_mx_y=vertcat(L2_mx_y2',L2_mx_y);

%cargar palnillas

cel_nodos={'Nodo','X','Y'};
cel_elementos={'Elemento','N1','N2','N3'};

for nod=1:1:size(P)(1)
  cel_nodos(nod+1,:)={nod,P(nod,1),P(nod,2)};
end

for elemen=1:1:size(TF)(1)
  cel_elementos(elemen+1,:)={elemen,TF(elemen,1),TF(elemen,2),TF(elemen,3)};
end

xlswrite('nodos.xlsx',cel_nodos);
xlswrite('elementos.xlsx',cel_elementos);

%graficar

%plot(x_DF,max(Z_DF)); hold on
%plot(L1_mx_y,L2_mx_y); hold on
%plot(x_A,max(Z_A)); xlabel ("Distancia [m]"); ylabel ("Desplazamiento[m]"); legend('Diferencias Finitas','Elementos Finitos','Función Analitica')


dir=[0 0 1];

%{
figure 1
tr1=trimesh(TF,Z_EF(:,1),Z_EF(:,2)-0.5,Z_EF(:,3)); 
rotate(tr1,dir,270)
hold on
tr2=trimesh(TF,Z_EF(:,1),Z_EF(:,2),Z_EF(:,3)); 
rotate(tr2,dir,360)
hold on
tr3=trimesh(TF,Z_EF(:,1)+1,Z_EF(:,2)+1,Z_EF(:,3)); 
rotate(tr3,dir,180)
hold on
tr4=trimesh(TF,Z_EF(:,1),Z_EF(:,2),Z_EF(:,3)); 
rotate(tr4,dir,90)
xlabel ("Distancia [m]"); ylabel ("Distancia [m]"); zlabel("Desplazamiento [m]");
view([45 45 130])

