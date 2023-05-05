function [CMatrix,weight]=microstates2cmatrix(onunits,I,J)
%only spikes occurred within 20ms of another spike is counted as a link
diffJ=diff(J);
elimind=find(diffJ>=10);
for n=1:length(elimind)
    if diffJ(elimind(n)+1)==0
    elimind=[elimind elimind(n)+1];
    end
end

I(elimind)=[];J(elimind)=[];

CMatrix=zeros(size(onunits,1),size(onunits,1));
weight=zeros(size(onunits,1),size(onunits,1));
emptyind=[];
for k=1:max(J)
    CMind{k}=I(J==k);
end
for k=1:length(CMind)
if isempty(CMind{k})==1
    emptyind=[emptyind,k];end
end
CMind(emptyind)=[];
for k=1:(length(CMind)-1)
CMatrix(CMind{k},CMind{k+1})=1;
weight(CMind{k},CMind{k+1})=J(k+1)-J(k)+1;
end

% removing self-connections
for n=1:size(onunits,1)
    CMatrix(n,n)=0;
end

end



