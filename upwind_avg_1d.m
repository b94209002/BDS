function c = upwind_avg_1d(c)
c =(circshift(c,-1) + c);