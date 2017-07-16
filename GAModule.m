function GAmodule

Generations=500;
Mutation=0.05;

format short g
close all
clc
sq=input('Please specify which steel quality You intend to use [N/mm^2] (235/275/355/420/460)! ');
nb=input('Please specify the number of support nodes! ');
nl=input('Please specify the number of nodes affected by a load! ');
nn=input('Please specify the total number of nodes! ');
fixedcoord=zeros(nl+nb,2);
for i=1:nb
    fixedcoord(i,1)=input(['Please specify the X-coordinate of support node ',num2str(i),'! (dm) ']);
    fixedcoord(i,2)=input(['Please specify the Y-coordinate of support node ',num2str(i),'! (dm) ']);
end

for i=nb+1:nb+nl
    fixedcoord(i,1)=input(['Please specify the X-coordinate of load node ',num2str(i-nb),'! (dm) ']);
    fixedcoord(i,2)=input(['Please specify the Y-coordinate of load node ',num2str(i-nb),'! (dm) ']);
end

xspan=max(fixedcoord(:,1))-min(fixedcoord(:,1));
yspan=max(fixedcoord(:,2))-min(fixedcoord(:,2));

if nn<nb+nl
    nn=nb+nl;
end

f=zeros(2*nn,1);

for i=nb+1:nb+nl
    f(2*i-1)=input(['For load node ',num2str(i-nb),' please specify the X-resultant of the load! (N) ']);
    f(2*i)=input(['For load node ',num2str(i-nb),' please specify the Y-resultant of the load! (N) ']);
end

Dof=[1:2:2*nn-1;2:2:2*nn]';
bc=[1:2*nb;zeros(1,2*nb)]';

Top=input('Please select method for topology optimization! (Press g for ground structure method, r for reduced method) ','s');

if Top=='g'
    neft=factorial(nn)/(2*factorial(nn-2));
    nevt=0;
    A=[];B=[];TOPGSM=[];
    for i=1:nn-1
        A=[A,i*ones(1,nn-i)];
        B=[B,(i+1:nn)];
    end
    for i=1:neft
        TOPGSM=[TOPGSM;i Dof(A(i),:) Dof(B(i),:)];
    end
else
    TOPGSM=0;
    neft=2*nn-3;
    maxnevt=factorial(nn)/(2*factorial(nn-2))-neft;
    nevt=input(['Please specify the number of elements with variable topology (0-',num2str(maxnevt),', ',num2str(round(maxnevt/10)),' recommended!) ']);
end
Population=input('Please specify the size of the initial population (150-1000 recomended)! ');

disp(' ')
disp('Calculating...')
tic

% Fabricational dimensions (in SI units):
% Element thickness:
t=1e-3*[2.5 3.0 3.2 4.0 4.9 5.0 2.5 3.0 3.2 4.0 4.9 5.0 6.0 6.3 3.0 3.2 4.0 4.9 5.0 6.0 6.3 8.0 3.0 3.2 3.6 4.0 4.9 5.0 6.0 6.3 7.1 8.0 3.2 3.6 4.0 4.9 5.0 6.0 6.3 7.1 8.0 3.2 3.6 4.0 4.9 5.0 5.6 6.0 6.3 7.1 8.0 3.6 4.0 4.9 5.0 5.6 6.0 6.3 7.1 8.0 3.6 4.0 4.9 5.0 5.6 6.0 6.3 7.1 8.0 10.0 4.0 4.9 5.0 5.6 6.0 6.3 7.1 8.0 8.8 10.0 12.0 12.5 4.9 5.0 5.6 6.0 6.3 7.1 8.0 8.8 10.0 12.0 12.5 4.9 5.0 5.6 6.0 6.3 7.1 8.0 8.8 10.0 12.0 12.5 16.0 5.0 5.6 6.0 6.3 7.1 8.0 8.8 10.0 12.0 12.5 14.2 16.0 5.0 5.6 6.0 6.3 7.1 8.0 8.8 10.0 12.0 12.5 14.2 16.0 5.0 5.6 6.0 6.3 7.1 8.0 8.8 10.0 12.0 12.5 14.2 16.0 5.0 5.6 6.0 6.3 7.1 8.0 8.8 10.0 12.0 12.5 14.2 16.0 6.0 6.3 7.1 8.0 8.8 10.0 12.0 12.5 14.2 16.0 6.0 6.3 7.1 8.0 8.8 10.0 12.0 12.5 14.2 16.0 8.0 8.8 10.0 12.0 12.5 14.2 16.0 8.0 8.8 10.0 12.0 12.5 14.2 16.0 20.0 19.0 22.0 25.0 22.0 25.0 12.0 16.0 19.0 22.0 25.0 28.0 32.0 12.0 16.0 19.0 22.0 25.0 28.0 32.0 36.0 16.0 19.0 22.0 25.0 28.0 32.0 36.0 40.0 25.0 28.0 32.0 36.0 40.0 25.0 28.0 32.0 36.0 40.0];


% Cross-sectional width:
w=1e-3*[40 40 40 40 40 40 50 50 50 50 50 50 50 50 60 60 60 60 60 60 60 60 70 70 70 70 70 70 70 70 70 70 76.2 76.2 76.2 76.2 76.2 76.2 76.2 76.2 76.2 80 80 80 80 80 80 80 80 80 80 90 90 90 90 90 90 90 90 90 100 100 100 100 100 100 100 100 100 100 120 120 120 120 120 120 120 120 120 120 120 120 140 140 140 140 140 140 140 140 140 140 140 150 150 150 150 150 150 150 150 150 150 150 150 160 160 160 160 160 160 160 160 160 160 160 160 180 180 180 180 180 180 180 180 180 180 180 180 200 200 200 200 200 200 200 200 200 200 200 200 250 250 250 250 250 250 250 250 250 250 250 250 260 260 260 260 260 260 260 260 260 260 300 300 300 300 300 300 300 300 300 300 350 350 350 350 350 350 350 400 400 400 400 400 400 400 400 350 350 350 400 400 450 450 450 450 450 450 450 500 500 500 500 500 500 500 500 550 550 550 550 550 550 550 550 600 600 600 600 600 700 700 700 700 700];
% Moment of inertia:
I=1e-8*[8.54 9.78 10.20 11.80 13.20 13.40 17.50 20.20 21.20 25.00 28.50 28.90 32.00 32.80 36.20 38.20 45.40 52.50 53.30 59.90 61.60 69.70 59.00 62.30 68.60 74.70 87.20 88.50 101.00 104.00 112.00 120.00 81.50 89.90 98.00 115.00 117.00 133.00 138.00 149.00 160.00 95.00 105.00 114.00 135.00 137.00 149.00 156.00 162.00 176.00 189.00 152.00 166.00 196.00 200.00 218.00 230.00 238.00 260.00 281.00 212.00 232.00 275.00 279.00 306.00 323.00 336.00 367.00 400.00 462.00 410.00 489.00 498.00 547.00 579.00 603.00 663.00 726.00 779.00 852.00 958.00 982.00 793.00 807.00 891.00 944.00 984.00 1086.00 1195.00 1287.00 1416.00 1609.00 1653.00 984.00 1002.00 1106.00 1174.00 1223.00 1352.00 1491.00 1608.00 1773.00 2023.00 2080.00 2430.00 1225.00 1353.00 1437.00 1499.00 1659.00 1831.00 1978.00 2186.00 2502.00 2576.00 2809.00 3028.00 1765.00 1954.00 2077.00 2168.00 2404.00 2661.00 2880.00 3193.00 3677.00 3790.00 4154.00 4504.00 2445.00 2710.00 2883.00 3011.00 3345.00 3709.00 4021.00 4471.00 5171.00 5336.00 5872.00 6394.00 4861.00 5399.00 5752.00 6014.00 6701.00 7455.00 8107.00 9055.00 10556.00 10915.00 12094.00 13267.00 6491.00 6788.00 7567.00 8423.00 9164.00 10242.00 11954.00 12365.00 13714.00 15061.00 10080.00 10547.00 11775.00 13128.00 14305.00 16026.00 18777.00 19442.00 21637.00 23850.00 21129.00 23055.00 25884.00 30435.00 31541.00 35211.00 38942.00 31857.00 34798.00 39128.00 46130.00 47839.00 53526.00 59344.00 71535.00 43360.00 48360.00 52890.00 74710.00 82150.00 65430.00 84070.00 97060.00 109200.00 120600.00 131200.00 144100.00 90750.00 117100.00 135500.00 153000.00 169400.00 184900.00 204000.00 221500.00 157700.00 183000.00 207100.00 230000.00 251600.00 278600.00 303500.00 326500.00 303400.00 332700.00 369400.00 403700.00 435500.00 494100.00 543500.00 606200.00 665400.00 721200.00];
% c:
c=w-2*t;
% Cross-sectional area:
A=w.^2.-c.^2;

% available=[A(1:3:150)' c(1:3:150)' t(1:3:150)' I(1:3:150)'];
available=[A' c' t' I'];
zerobar=[eps*ones(30,1) zeros(30,3)];
available=[zerobar;available];
% Number of available elements:
na=size(available,1);
% Number of bits required in the different segments of the bit string:
numbitEp=ceil(log2(na));
if nn == nb+nl
    numbitCoordY=0;
    numbitCoordX=0;

else
    numbitCoordX=ceil(log2(xspan));
    numbitCoordY=ceil(log2(yspan));
end
if Top == 'g'
    numbitEdof=0;
else
    numbitEdof=ceil(log2(nn));
end

% Creating an initial population:
Initpop=round(rand(Population,(nn-nb-nl)*(numbitCoordX+numbitCoordY)+(neft+nevt)*numbitEp+2*nevt*numbitEdof));

options = gaoptimset('PopulationType', 'bitString',...
    'FitnessLimit',0,...
    'InitialPopulation', Initpop,...
    'PlotFcns', {@gaplotbestf2},...
    'Generations', Generations,...
    'PopulationSize', Population,...
    'StallGenLimit', Inf,...
    'StallTimeLimit', Inf,...
    'SelectionFcn', @selectiontournament,...
    'FitnessScalingFcn', @fitscalingrank,...
    'EliteCount', 2,...
    'CrossoverFraction', 0.9,...
    'CrossoverFcn', @crossoverscattered,...
    'MutationFcn', {@mutationuniform, Mutation}...
    );


E=210e9;
[x,fval]=ga(@(x) penalty2(x,available,bc,Dof,E,f,fixedcoord,na,nb,neft,nevt,nl,nn,numbitCoordX,numbitCoordY,numbitEdof,numbitEp,sq,TOPGSM,xspan,yspan),(nn-nb-nl)*(numbitCoordX+numbitCoordY)+(neft+nevt)*numbitEp+2*nevt*numbitEdof,options);
Time=toc;

[Coord,Edof,Ep]=bintranslate(available,x,Dof,fixedcoord,na,nb,neft,nevt,nl,nn,numbitCoordX,numbitCoordY,numbitEdof,numbitEp,xspan,yspan);
if Edof==0;
    Edof=TOPGSM;
end
[Ed,Ex,Ey,v1,v2,w,wrong]=FEM2(bc,Coord,Dof,E,Edof,Ep,f,neft,nevt,sq,TOPGSM);
if wrong==1
    disp('Unable to find a feasible solution!')
else


% Removing zero-bars from the FEM parameters before the plot:
for i=1:neft+nevt
    if v2(i,1)==0
        P(i)=0;
    else
        P(i)=i;
    end
end
J=find(P);

Ex2=zeros(length(J),2);
Ey2=zeros(length(J),2);
Ed2=zeros(length(J),4);
Edof2=zeros(length(J),5);
Ep2=zeros(length(J),4);
v22=zeros(length(J),2);

for i=1:length(J)
    Ex2(i,:)=Ex(J(i),:);
    Ey2(i,:)=Ey(J(i),:);
    Ed2(i,:)=Ed(J(i),:);
    Edof2(i,:)=Edof(J(i),:);
    Ep2(i,:)=Ep(J(i),:);
    v22(i,:)=v2(J(i),:);
end


% Displaying the calculated truss:
figure
eldraw2(Ex2,Ey2,[1 3 1],Edof2)
axis([0 xspan/10+1 -1 yspan/10+1]);
legend(['Elapsed time: ',num2str(Time),' seconds; Fitness: ',num2str(fval),'; Actual weight: ',num2str(w),' Lbs.'])

figure
eldraw2(Ex2,Ey2,[1 3 1],Edof2)
axis([0 xspan/10+1 -1 yspan/10+1]);
[sfac]=scalfact2(Ex2,Ey2,Ed2,0.1);
eldisp2(Ex2,Ey2,Ed2,[2 1 1],sfac)
pltscalb2(sfac,[5e-2 10 8]);


% Printing the results:
disp(['Element no.','      ','Dimensions in mm (w*w*t)','      ','Start coordinates(m)','      ','End coordinates (m)','      ','Stress in % of EC3'])
for i=1:length(J)
    disp(['      ',num2str(J(i)),'      ',num2str((Ep2(i,2)+2*Ep2(i,3))*1000),'x',num2str((Ep2(i,2)+2*Ep2(i,3))*1000),'x',num2str(Ep2(i,3)*1000),'      ',num2str(Ex2(i,1)),',',num2str(Ey2(i,1)),'      ',num2str(Ex2(i,2)),',',num2str(Ey2(i,2)),'      ',num2str(v22(i,1)/v22(i,2)*100)])
end

disp(' ')
disp(' ')
disp('Element stresses in MPa:')
for i=1:length(J)
disp([num2str(J(i)),' ',num2str((v22(i,1)/Ep2(i,1))/1e6)])
end
disp(' ')
disp(' ')
disp('Displacements:')
disp(['Horizontally','      ',num2str(v1(1,1)*1000),'      ','mm,','      ',num2str(100*v1(1,1)/v1(1,2)),'% of the limit'])
disp(['Vertically','      ',num2str(v1(2,1)*1000),'      ','mm,','      ',num2str(100*v1(2,1)/v1(2,2)),'% of the limit'])
disp(' ')

disp(' ')
disp(['Weight of the structure is ',num2str(w),'Lbs or ',num2str(w*0.45359237),'kg'])
end