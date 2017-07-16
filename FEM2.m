function [Ed,Ex,Ey,v1,v2,w,wrong]=FEM2(bc,Coord,Dof,E,Edof,Ep,f,neft,nevt,sq,TOPGSM)

density=7850;
w=0;
wrong=0;
Ed=0;
Ex=0;
Ey=0;
v1=0;
v2=0;
ne=neft+nevt;

% element coordinates:
[Ex,Ey]=coordxtr(Edof,Coord,Dof,2);

% Checking if double nodes exists:
Ucoord=unique(Coord,'rows');

if size(Ucoord,1) ~= size(Coord,1)
    wrong=1;
    return
end

% Looking for elements with identical topology, and elements without
% any spatial extension:
if TOPGSM == 0
    PD=Edof(:,2).*Edof(:,4);
    ND=unique(PD);
    for i=1:ne
        if PD(i)==Edof(i,2)^2
            wrong=1;
            return
        end
    end
    if length(PD) ~= length(ND)
        wrong=1;
        return
    end
end

% Stiffness matrix of the structure:
K=zeros(length(f));
for i=1:ne
    Ke=bar2e(Ex(i,:),Ey(i,:),[E Ep(i,1)]);
    K=assem(Edof(i,:),K,Ke);
end
if rcond(K) == NaN
    wrong=1;
    return
end

% Stability check:
if abs(det(K))<1e-9
    wrong=1;
    return
end


% Calculating the element displacements:
[a,r]=solveq(K,f,bc);
Ed=extract(Edof,a);


% Calculating max displacements and limits and storing the info in v1:
% (All odd-numbered DOFs are x-directions)
Xdispl=[];
Ydispl=[];
for i=1:2:length(a)-1
    Xdispl=[Xdispl a(i)];
end
for i=2:2:length(a)
    Ydispl=[Ydispl a(i)];
end


% Lx (only valid for simple span structures):
Lx=max(max(Ex)); Ly=max(max(Ey)); v1=zeros(2,2);
v1(1,1)=max(abs([max(Xdispl) min(Xdispl)]));
v1(2,1)=max(abs([max(Ydispl) min(Ydispl)]));

% Horizontal limit
v1(1,2)=Ly/250;

% Vertical limit
v1(2,2)=Lx/250;


% Calculating cross-sectional class and effective area if needed:
epsilon=sqrt(235/sq);
Ep2=zeros(ne,2);
lambdap=zeros(ne,1);
ro=zeros(ne,1);
for i=1:ne
    if Ep(i,1) == eps
        Ep2(i,:)=[0 123];
    else
        lambdap(i)=(Ep(i,2)/Ep(i,3))/(28.4*epsilon*2);
        ro(i)=(lambdap(i)-0.055*4)/lambdap(i)^2;
        if ro(i)>1
            ro(i)=1;
        end
        if Ep(i,2)/Ep(i,3) <= 42*epsilon
            Ep2(i,1)=Ep(i,1);
            Ep2(i,2)=123;
        else
            Ep2(i,1)=Ep(i,1)*ro(i);
            Ep2(i,2)=4;
        end
    end
end


% Calculating the weight of the structure:
for i=1:ne
    % Calculating element length:
    Le(i)=sqrt((Ex(i,1)-Ex(i,2))^2+(Ey(i,1)-Ey(i,2))^2);
    w=w+Ep(i,1)*Le(i)*density/0.45359237;
end


% Calculating buckling stability:
if sq == 235||sq == 275||sq == 355||sq == 420
    alpha=0.21;
else
    alpha=0.13;
end

Ncr=zeros(ne,1);
Nbrd=zeros(ne,1);
lambda=zeros(ne,1);
Phi=zeros(ne,1);
Chi=zeros(ne,1);
Ned=zeros(ne,1);
v2=zeros(ne,2);


for i=1:ne
    % Critical load:
    Ncr(i)=pi^2*E*Ep(i,4)/Le(i)^2;
    if Ncr(i) == 0
        v2(i,:)=[0 1];
    else
        % Buckling reduction factor Chi:
        lambda(i)=sqrt((Ep2(i,1)*sq*1e6)/Ncr(i));
        Phi(i)=0.5*(1+alpha*(lambda(i)-0.2)+lambda(i)^2);
        Chi(i)=1/(Phi(i)+sqrt(Phi(i)^2-lambda(i)^2));
        if Chi(i)>1
            Chi(i)=1;
        end
        % Calculating element forces (tensile=positive) and constraints:
        Ned(i)=bar2s(Ex(i,:),Ey(i,:),[E Ep(i,1)],Ed(i,:));
        if lambda(i)<=0.2||abs(Ned(i))/Ncr(i)<=0.04
            Nbrd(i)=Ep(i,1)*sq*1e6;
        else
            Nbrd(i)=Chi(i)*Ep2(i,1)*sq*1e6;
        end
        % Storing current forces and constraints in v2:
        if Ned(i)<0
            v2(i,:)=[Ned(i) -Nbrd(i)];
        else
            v2(i,:)=[Ned(i) Ep(i,1)*sq*1e6];
        end
    end
end