#!/usr/bin/octave
clear
#Entry for 1D 2D or 3D
do
 D=input("Enter 1 for 1D problem, 2 for 2D, 3 for 3D= ");
until (D>=1 && D<=3)	#Check entry is between 1 and 3

E=input("Enter number of elements= ");	#Entry for number of elements
N=input("Enter number of nodes= ");	#Entry for number of nodes
K=zeros(N*D);	#Intialize global stiffness matrix to zero
e=zeros(E,2);	#Intialize nodes entries to zero

#Entry for node numbers
for i=1:E
 printf("Enter node number for element %u \n", i);
 e(i,1)=input("Enter for 1st node= ");
 e(i,2)=input("Enter for 2nd node= ");
endfor

#Entry for Area and Modulus
c=input("Enter if all elements have same Area & Modulus (y/n)? ","s"); 
if(c=="y")	#Entry only once if 'y' is entered
 printf("Enter Area & Modulus \n"); 
 A(1:E)=input("Area = ");
 E1(1:E)=input("Modulus = ");
else
 for i=1:E	#Entry for individual elements
  printf("Enter Area & Modulus for element %u \n", i);
  A(i)=input("Area = ");
  E1(i)=input("Modulus = ");
 endfor
endif

#Entry for Length and Cx, Cy, Cz
if(D==1)	#Cx & Length set for 1D system
 Cx=zeros(E,1); #Intialize Cx entries to zero
 Cx(1:E)=1; #Intialize all Cx entries to zero
 for i=1:E #Length entries
  printf ("For element %u ", i)
  do	#safety feature
   s=input("enter 1 for Cartesian coordinates or 2 for Length");	#Choice between Coordinates or Length
  until (s==1 || s==2)
  if (s==1) #Coordinate Entry
   x1=input("x in 1st node= ");	#1st node coordinate x
   x2=input("x in 2nd node= ");	#2nd node coordinate x
   L(i)=((x2-x1)^2)^.5;	#Length calculation from coordinates
  elseif (s==2)	#Length and Angle Entry
   L(i)=input("Length = ");
  endif
  k(i)=A(i)*E1(i)/L(i);
 endfor
 
 elseif (D==2)	#Cx, Cy & Length x, y set for 2D system
 for i=1:E
  printf ("For element %u ", i)
  do	#safety feature
   s=input("enter 1 for Cartesian coordinates or 2 for Length and angles= ");	#Choice between Coordinates or Length & angles
  until (s==1 || s==2)
  if (s==1) #Coordinate Entry
   x1=input("x in 1st node= ");	#1st node coordinate x
   y1=input("y in 1st node= ");	#1st node coordinate y
   x2=input("x in 2nd node= "); #2nd node coordinate x
   y2=input("y in 2nd node= "); #2nd node coordinate y
   L(i)=((x2-x1)^2+(y2-y1)^2)^.5;	#Length calculation from coordinates
   Cx(i)=(x2-x1)/L(i);	#cos x
   Cy(i)=(y2-y1)/L(i);	#cos y
  elseif (s==2)	#Length and Angle Entry
   L(i)=input("Length = ");
   thx=input("Theta-x with respect to global axis = ");
   thy=90-thx;
   Cx(i)=cosd(thx);
   Cy(i)=cosd(thy);
  endif
  k(i)=A(i)*E1(i)/L(i);
 endfor
else  #Cx, Cy, Cz set for 3D system
 for i=1:E
  printf ("For element %u ", i)
  do	#safety feature
   s=input("enter 1 for Cartesian coordinates or 2 for Length and angles= ");	#Choice between Coordinates or Length & angles
  until (s==1 || s==2)
  if (s==1) #Coordinate Entry
   x1=input("x in 1st node= ");	#1st node coordinate x
   y1=input("y in 1st node= ");	#1st node coordinate y
   z1=input("z in 1st node= ");	#1st node coordinate z
   x2=input("x in 2nd node= ");	#2nd node coordinate x
   y2=input("y in 2nd node= ");	#2nd node coordinate y
   z2=input("z in 2nd node= ");	#2nd node coordinate z
   L(i)=((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^.5;	#Length calculation from coordinates
   Cx(i)=(x2-x1)/L(i);	#cos x
   Cy(i)=(y2-y1)/L(i);	#cos y
   Cz(i)=(z2-z1)/L(i);	#cos z
  elseif (s==2)	#Length and Angle Entry
   L(i)=input("Length = ");
   do	#Entry of values of angles with checking
    thx=input("Theta-x with respect to global axis = ");
    thy=input("Theta-y with respect to global axis = ");
    thz=input("Theta-z with respect to global axis = ");
    Cx(i)=cosd(thx);
    Cy(i)=cosd(thy);
    Cz(i)=cosd(thz);
    l=(Cx(i)^2+Cy(i)^2+Cz(i)^2)^.5;
    if(l!=1)	#Error output if 'l' is not equal to 1
     printf("Error! Please enter angles again\n");
	endif
   until(l==1) #Checking Done
  endif
  k(i)=A(i)*E1(i)/L(i);
 endfor
endif

#Calclating lamda for each dimentional system and then calculate local stiffness and added top global stiffness
for i=1:E	
 if (D==1)	#Calclating lamda for 1D
  lamda=Cx(i)^2;
 elseif (D==2)	#Calclating lamda for 2D
  lamda=[Cx(i)^2, Cx(i)*Cy(i); Cx(i)*Cy(i), Cy(i)^2];
 else 	#Calclating lamda for 3D
  lamda=[Cx(i)^2, Cx(i)*Cy(i), Cx(i)*Cz(i); Cx(i)*Cy(i), Cy(i)^2, Cy(i)*Cz(i);  Cx(i)*Cz(i),  Cy(i)*Cz(i), Cz(i)^2];
 endif
 
 k1=zeros(N*D);	#Intializing local stiffness matrix to zero
 a=e(i,1)*D-D+1;	#Initial Row/Column value 1
 b=e(i,1)*D; #End Row/Column Value 1
 a1=e(i,2)*D-D+1;	#Initial Row/Column value 2
 b1=e(i,2)*D;	#End Row/Column Value 2
 k1(a:b,a:b)=lamda;	#Building local stiffness matrix
 k1(a1:b1,a1:b1)=lamda;
 k1(a:b,a1:b1)=-lamda;
 k1(a1:b1,a:b)=-lamda;
 k1=k(i)*k1;	#Multiplying to stiffness of element
 K=K+k1;	#Adding local stiffness matrix to global stiffness matrix
 k2(:,:,i)=k1;	#keeps track of local stiffness matrix for each element (unused code)
endfor

K1=K;	#Storing Values of K in K1 for future use
N1=zeros(N,1);	#Stores values of Node numbers with restraints later on
F2=zeros(N*D,1);	#Calculates forces resulting from displacement later on
disp=zeros(N*D,1);	#Stores dispalcements later on

i1=input("Enter the number of nodes with restraints= ")
n1=zeros(i1*D,1);	#keeps track of Column numbers with restraints

if(D==1)
 for i=1:i1
  printf("Enter the restraint node number %u =", i);
  n1(i)=input("");
  disp(n1(i))=input("Displacement = ");
  F2=F2-K(:,n1(i))*disp(n1(i));
  endfor
  K1([n1],:)=[];
  K1(:,[n1])=[];
endif
if(D==2)
 for i=1:i1
  printf("Enter the restraint node number %u =", i);
  N1(i)=input("");
  cx=input("Is it restrainted x direction?(y/n) ","s");
  if(cx=="y")
   n1(i*D-1)=N1(i)*D-1;
   disp(N1(i)*D-1)=input("Displacement = ");
   F2=F2-K(:,N1(i)*D-1)*disp(N1(i)*D-1);
  endif
  cy=input("Is it restrainted y direction?(y/n) ","s");
  if(cy=="y")
   n1(i*D)=N1(i)*D;
   disp(N1(i)*D)=input("Displacement = ");
   F2=F2-K(:,N1(i)*D)*disp(N1(i)*D);
  endif
 endfor
 n1(~any(n1,2),:)=[];
 K1([n1],:)=[];
 K1(:,[n1])=[];
endif
if(D==3)
 for i=1:i1
  printf("Enter the restraint node number %u =", i);
  N1(i)=input("");
  cx=input("Is it restrainted x direction?(y/n) ","s");
  if(cx=="y")
   n1(i*D-2)=N1(i)*D-2;
   disp(N1(i)*D-2)=input("Displacement = ");
   F2=F2-K(:,N1(i)*D-2)*disp(N1(i)*D-2);
  endif
  cy=input("Is it restrainted y direction?(y/n) ","s");
  if(cy=="y")
    n1(i*D-1)=N1(i)*D-1;
    disp(N1(i)*D-1)=input("Displacement = ");
    F2=F2-K(:,N1(i)*D-1)*disp(N1(i)*D-1);
  endif
  cz=input("Is it restrainted z direction?(y/n) ","s");
  if(cz=="y")
   n1(i*D)=N1(i)*D;
   disp(N1(i)*D)=input("Displacement = ");
   F2=F2-K(:,N1(i)*D)*disp(N1(i)*D);
  endif
 endfor
 n1(~any(n1,2),:)=[];
 K1([n1],:)=[];
 K1(:,[n1])=[];
endif 

N2=zeros(N,1);
F=zeros(N*D,1);
n2=transpose(linspace(1,N*D,N*D));
f1=input("Enter the number of nodes with forces= ")


if(D==1)
 for i=1:f1
  printf("Enter the force node number %u =", i);
  N2(i)=input("");
  F(N2(i))=input("Enter the force =");
 endfor
 
endif
if(D==2)
 for i=1:f1
 printf("Enter the force node number %u =", i);
  N2(i)=input("");
  cx=input("Is force in x direction?(y/n) ","s");
  if(cx=="y")
   F(N2(i)*D-1)=input("Enter the force =");
  endif
  cy=input("Is force in y direction?(y/n) ","s");
  if(cy=="y")
   F(N2(i)*D)=input("Enter the force =");
  endif
 endfor
 
endif
if(D==3)
 for i=1:f1
 printf("Enter the force node number %u =", i);
  N2(i)=input("");
  cx=input("Is force in x direction?(y/n) ","s");
  if(cx=="y")
   F(N2(i)*D-2)=input("Enter the force =");
  endif
  cy=input("Is force in y direction?(y/n) ","s");
  if(cy=="y")
   F(N2(i)*D-1)=input("Enter the force =");
  endif
   cz=input("Is force in y direction?(y/n) ","s");
  if(cz=="y")
   F(N2(i)*D)=input("Enter the force =");
  
  endif
 endfor
endif
F1=F;
F1=F1+F2;
F1([n1],:)=[];
n2([n1],:)=[];
d=zeros(N*D,1);
d1=inv(K1)*F1;
d([n2],:)=d1;
d([n1],:)=disp([n1],:);

F=K*d

i=transpose(linspace(1, N*D-D+1, N));
j=transpose(linspace(2, N*D-D+2, N));
k3=transpose(linspace(2, N*D-D+3, N));
if(D>=1)
m=d(i);
save Dx.txt m
m=F(i);
save Fx.txt m
endif
if(D>=2)
m=d(j);
save Dy.txt m
m=F(j);
save Fy.txt m
endif
if(D==3)
m=d(k3);
save Dz.txt m
m=F(k3);
save Fz.txt m
endif

for i=1:E
 if(D==1)
  stress(i) =E1(i)/L(i)*[-Cx(i),Cx(i)]*d([e(i,1)*D,e(i,2)*D]);
 elseif(D==2)
  stress(i) =E1(i)/L(i)*[-Cx(i),-Cy(i), Cx(i), Cy(i)]*d([e(i,1)*D-1:e(i,1)*D,e(i,2)*D-1:e(i,2)*D]);
 elseif(D==3)
  stress(i) =E1(i)/L(i)*[-Cx(i),-Cy(i), -Cz(i), Cx(i), Cy(i), Cz(i)]*d([e(i,1)*D-2:e(i,1)*D,e(i,2)*D-2:e(i,2)*D]);
 endif
endfor

save stress.txt stress

for i=1:E
 printf("Stress of element %u ", i);
stress(i)
 endfor
