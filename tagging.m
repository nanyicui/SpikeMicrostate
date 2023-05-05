function [S1CMatrix,WCMatrix,S2CMatrix,S1nonemptystatesrate,S2nonemptystatesrate]=tagging(S1CMatrix,WCMatrix,S2CMatrix)

% Tagging CMatrices
for a=1:length(WCMatrix)
    WCMatrix{2,a}=a;
end
% Unique Tagging CMatrices
for n=1:(length(WCMatrix)-1)
    for m=1:length(WCMatrix)
        
        if isequal(WCMatrix(1,n),WCMatrix(1,m));
            
            WCMatrix{2,m}=WCMatrix{2,n};
        end
    end
end

% Plotting histogram of tags
%hist(cell2mat(WCMatrix(2,:)),0.5:1:(length(WCMatrix)+0.5))
title('Histogram of Unique microstates in W')
xlabel('Unique tags ID')
% Tagging CMatrix in another episode accordingly given the existing tags

for n=1:length(WCMatrix)
    for m=1:length(S1CMatrix)
        if isequal(WCMatrix{1,n},S1CMatrix{1,m});
            S1CMatrix{2,m}=WCMatrix{2,n};
        end
    end
end
%hist(cell2mat(S1CMatrix(2,:)),0.5:1:(max(cell2mat(S1CMatrix(2,:)))+0.5))
title('Histogram of Unique microstates in S1')
xlabel('Unique tags ID')
for n=1:length(WCMatrix)
    for m=1:length(S2CMatrix)
        if isequal(WCMatrix{1,n},S2CMatrix{1,m});
            S2CMatrix{2,m}=WCMatrix{2,n};
        end
    end
end
%hist(cell2mat(S2CMatrix(2,:)),0.5:1:(max(cell2mat(S2CMatrix(2,:)))+0.5))
title('Histogram of Unique microstates in S2')
xlabel('Unique tags ID')


S1nonemptystates=0;
for n=1:length(S1CMatrix)
    S1nonemptystates=S1nonemptystates+isempty(S1CMatrix{2,n});end
S1nonemptystatesrate=S1nonemptystates/length(S1CMatrix);
S2nonemptystates=0;
for n=1:length(S2CMatrix)
    S2nonemptystates=S2nonemptystates+isempty(S2CMatrix{2,n});end
S2nonemptystatesrate=S2nonemptystates/length(S2CMatrix);

end