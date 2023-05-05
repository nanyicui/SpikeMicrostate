function [onperiodstartend]=findonstartend(off_periods)
off_periods(end)=0;
seq=find(off_periods==1);
sequ=diff(seq);
seque= sequ~=1;
offend=seq(seque);offend=[offend; seq(end)];offend(1)=[];
seq0=find(off_periods==0);
sequ0=diff(seq0);
seque0= sequ0~=1;
offstart=seq0(seque0)+1;
if offstart(2)<offend(1)
    offstart(1)=[];
elseif offstart(2)==offend(1)
    offstart(2)=[];
end
offstart=offstart-4;
offend=offend+4;
onend=offstart-1;
onstart=offend+1;onstart=[1;onstart];onstart(end)=[];
onperiodstartend=[onstart,onend];
end