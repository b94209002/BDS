function [minLL maxLL minRL maxRL minLH maxLH minRH maxRH] = BDS_bilinear_compute_rnage(sh);

srt = circshift(sh,[1,1]);
srb = circshift(sh,[1,-1]);
sr  = circshift(sh,[1,0]);
sl  = circshift(sh,[-1,0]);
slt = circshift(sh,[-1,1]);
slb = circshift(sh,[-1,-1]);
st  = circshift(sh,[0,1]);
sb  = circshift(sh,[0,-1]);

minLL = min(slb,sb,sl,sh);maxLL = max(slb,sb,sl,sh);
minRL = min(sb,srb,sr,sh);maxRL = max(sb,srb,sr,sh);
minLH = min(st,slt,sl,sh);maxLH = max(st,slt,sl,sh);
minRH = min(st,srt,sr,sh);maxRH = max(st,srt,sr,sh);