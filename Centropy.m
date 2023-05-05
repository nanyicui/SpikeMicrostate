% CMatrix entropy calculation
%convert cell to mat
for n=1:length(S1CMatrix)  
   S1CMatrixm(:,:,n)=S1CMatrix{n};
end
for n=1:length(S2CMatrix)  
   S2CMatrixm(:,:,n)=S2CMatrix{n};
end
for n=1:length(WCMatrix)  
   WCMatrixm(:,:,n)=WCMatrix{n};
end

%using MIToolbox
S1CMatrixm=reshape(S1CMatrixm,14*14,length(S1CMatrixm));
S1CMatrixmsum=sum(S1CMatrixm,2);S1CMatrixwco=[];
for n=1:14*14
    S1CMatrixwco=[S1CMatrixwco; n*ones(S1CMatrixmsum(n),1)];
end
HS1CMatrix=h(S1CMatrixwco);
S2CMatrixm=reshape(S2CMatrixm,14*14,length(S2CMatrixm));
S2CMatrixmsum=sum(S2CMatrixm,2);S2CMatrixwco=[];
for n=1:14*14
    S2CMatrixwco=[S2CMatrixwco; n*ones(S2CMatrixmsum(n),1)];
end
HS2CMatrix=h(S2CMatrixwco);
WCMatrixm=reshape(WCMatrixm,14*14,length(WCMatrixm));
WCMatrixmsum=sum(WCMatrixm,2);WCMatrixwco=[];
for n=1:14*14
    WCMatrixwco=[WCMatrixwco; n*ones(WCMatrixmsum(n),1)];
end
HWCMatrix=h(WCMatrixwco);
%motif entropy calculation
HS1motifmat=0;
S1pmotifmat=sum(S1motifmat,3)/sum(sum(sum(S1motifmat)));
for m=1:13
    for n=1:14
        if S1pmotifmat(m,n)==0
            S1pmotifmat(m,n)=1;
        end
        HS1motifmat=HS1motifmat-S1pmotifmat(m,n)*log(S1pmotifmat(m,n));
    end
end
