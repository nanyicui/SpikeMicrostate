function [off_periods, sum_windowsweep]=classifier_onoff_mssim(Units,swpwinsize,centreunq)
% %% window sweeping all units using a uniform distribution
% 
% for n=((swpwinsize+1)/2):(size(Units,2)-(swpwinsize-1)/2)
%     
%     windowcontour(:,:,n)=convn(Units(:,(n-(swpwinsize-1)/2):(n+(swpwinsize-1)/2)),ones(size(Units,1),swpwinsize));
%     
% end

%% window sweeping all units using a Gaussian distribution
%swpwinsize must be positive odd integer
for n=((swpwinsize+1)/2):(size(Units,2)-(swpwinsize-1)/2)
    
    windowcontour(:,:,n)=convn(Units(:,(n-(swpwinsize-1)/2):(n+(swpwinsize-1)/2)),repmat(pdf('Normal',-((swpwinsize-1)/2):((swpwinsize-1)/2),0,1.2),size(Units,1),1));
    
end

mu_contour=mean(windowcontour,3); % mean matrix of windowsweeping as kernel
%% multidimensional deconvolution

for m=1:size(mu_contour,2)
de_mu_contour(:,m)=deconv(mu_contour(:,m),ones(size(Units,1),1));
end

for o=1:size(de_mu_contour,1)
deconv_mu_contour(o,:)=deconv(de_mu_contour(o,:),ones(1,(size(de_mu_contour,2)+1)/2));
end

%% summing values in window by multiplication of deconv_mu_contour
for n=((swpwinsize+1)/2):(size(Units,2)-(swpwinsize-1)/2)
    
    windowsweep(:,:,n)=deconv_mu_contour'*Units(:,(n-(swpwinsize-1)/2):(n+(swpwinsize-1)/2));
    
end

sum_windowsweep=squeeze(sum(sum(windowsweep,1),2))';
unq_windowsweep=unique(sum_windowsweep);

off_periods=sum_windowsweep<centreunq(1);%criteria bottom ~%, 1 is OFF;

off_periods=smooth(off_periods,'rloess');off_periods=off_periods>0.5;


end