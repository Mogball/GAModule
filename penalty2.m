function [fitness]=penalty2(chromosome,available,bc,Dof,E,f,fixedcoord,na,nb,neft,nevt,nl,nn,numbitCoordX,numbitCoordY,numbitEdof,numbitEp,sq,TOPGSM,xspan,yspan)

[Coord,Edof,Ep]=bintranslate(available,chromosome,Dof,fixedcoord,na,nb,neft,nevt,nl,nn,numbitCoordX,numbitCoordY,numbitEdof,numbitEp,xspan,yspan);
if Edof==0;
    Edof=TOPGSM;
end
[Ed,Ex,Ey,v1,v2,w,wrong]=FEM2(bc,Coord,Dof,E,Edof,Ep,f,neft,nevt,sq,TOPGSM);
P=zeros(3+neft+nevt,1);
if wrong==1
    P(1)=1e20;
    fitness=w+sum(P);
    return
else
    P(1)=0;
    dispfactors=v1(:,1)./v1(:,2);
    stressfactors=v2(:,1)./v2(:,2);
    for i=1:2
        if dispfactors(i)>1
            P(i+1)=dispfactors(i)*1e8;
        else
            P(i+1)=0;
        end


    end
    for i=1:neft+nevt
        if stressfactors(i)>1
            P(i+3)=stressfactors(i)*1e8;
        else
            P(i+3)=0;
        end
    end
    fitness=w+sum(P);
end