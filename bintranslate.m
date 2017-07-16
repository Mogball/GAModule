function [Coord,Edof,Ep]=bintranslate(available,chromosome,Dof,fixedcoord,na,nb,neft,nevt,nl,nn,numbitCoordX,numbitCoordY,numbitEdof,numbitEp,xspan,yspan)

format short
% Extracting nodal coordinates
if numbitCoordX==0
    Coord=1e-1*fixedcoord;
else
    chromosome1=chromosome(1:numbitCoordX*(nn-nb-nl));
    chromosome2=chromosome(1+numbitCoordX*(nn-nb-nl):numbitCoordX*(nn-nb-nl)+numbitCoordY*(nn-nb-nl));


    variablecoord=zeros(nn-(nb+nl),2);
    for i=numbitCoordX:numbitCoordX:numbitCoordX*(nn-nb-nl)
        if round(bin2dec(num2str(chromosome1(i-(numbitCoordX-1):i))))>xspan
            variablecoord(i/numbitCoordX,1)=xspan-1;
        else
            variablecoord(i/numbitCoordX,1)=round(bin2dec(num2str(chromosome1(i-(numbitCoordX-1):i))));
        end
    end

    for i=numbitCoordY:numbitCoordY:numbitCoordY*(nn-nb-nl)
        if round(bin2dec(num2str(chromosome2(i-(numbitCoordY-1):i))))>yspan
            variablecoord(i/numbitCoordY,2)=yspan;
        else
            variablecoord(i/numbitCoordY,2)=round(bin2dec(num2str(chromosome2(i-(numbitCoordY-1):i))));
        end
    end
    Coord=1e-1*[fixedcoord;variablecoord];
end
chromosome3=chromosome(1+(numbitCoordX+numbitCoordY)*(nn-nb-nl):numbitEp*(neft+nevt)+(numbitCoordX+numbitCoordY)*(nn-nb-nl));
chromosome4=chromosome(1+length(chromosome)-2*nevt*numbitEdof:length(chromosome)-nevt*numbitEdof);
chromosome5=chromosome(1+length(chromosome)-nevt*numbitEdof:length(chromosome));

if numbitEdof == 0
    Edof=0;
else
    add=0.01*[1:50]';
    Coord1=Coord(:,1)+add(1:size(Coord,1));

    sortcoord=sort(Coord1);

    % Creating a basic structure
    A=zeros(neft,2);B=zeros(neft,2);
    A(neft,:)=Dof(find(Coord1==sortcoord(end-1),1),:);
    B(1,:)=Dof(find(Coord1==sortcoord(2),1),:);
    for i=2:2:neft-1
        A(i-1,:)=Dof(find(Coord1==sortcoord(i/2),1),:);
        A(i,:)=Dof(find(Coord1==sortcoord(i/2),1),:);
        B(i,:)=Dof(find(Coord1==sortcoord((i+4)/2),1),:);
        B(i+1,:)=Dof(find(Coord1==sortcoord((i+4)/2),1),:);
    end


    % Extracting the elements with variable topology
    if nevt==0
        Edof=[(1:neft)' A B];
    else
        C=zeros(nevt,2);D=zeros(nevt,2);
        for i=numbitEdof:numbitEdof:numbitEdof*nevt
            if ceil(bin2dec(num2str(chromosome4(i-(numbitEdof-1):i))))>nn
                C(i/numbitEdof,:)=Dof(nn,:);
            elseif ceil(bin2dec(num2str(chromosome4(i-(numbitEdof-1):i))))==0
                C(i/numbitEdof,:)=Dof(1,:);
            else
                C(i/numbitEdof,:)=Dof(ceil(bin2dec(num2str(chromosome4(i-(numbitEdof-1):i)))),:);
            end
            if ceil(bin2dec(num2str(chromosome5(i-(numbitEdof-1):i))))>nn
                D(i/numbitEdof,:)=Dof(nn,:);

            elseif ceil(bin2dec(num2str(chromosome5(i-(numbitEdof-1):i))))==0
                D(i/numbitEdof,:)=Dof(1,:);
            else
                D(i/numbitEdof,:)=Dof(ceil(bin2dec(num2str(chromosome5(i-(numbitEdof-1):i)))),:);
            end
        end

        Edof=[(1:neft)' A B;(1+neft:neft+nevt)' C D];
    end
end
% Extracting the element properties
Ep=zeros(neft+nevt,4);
for i=numbitEp:numbitEp:(neft+nevt)*numbitEp
    if ceil(bin2dec(num2str(chromosome3(i-(numbitEp-1):i))))>na
        Ep(i/numbitEp,:)=available(na,:);
    elseif ceil(bin2dec(num2str(chromosome3(i-(numbitEp-1):i))))==0
        Ep(i/numbitEp,:)=available(1,:);
    else
        Ep(i/numbitEp,:)=available(ceil(bin2dec(num2str(chromosome3(i-(numbitEp-1):i)))),:);
    end
end