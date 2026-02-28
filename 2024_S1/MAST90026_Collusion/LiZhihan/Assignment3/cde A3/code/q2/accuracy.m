function p=accuracy(t,N,K)
[~,s31]=CN(t,N,K);
[~,s32]=CN(t,N,2*K);
[~,s33]=CN(t,N,4*K);
p3=(log(abs(s33-s32))-log(abs(s32-s31)))/log(0.5);
[~,s41]=TRBFD(t,N,K);
[~,s42]=TRBFD(t,N,2*K);
[~,s43]=TRBFD(t,N,4*K);
p4=(log(abs(s43-s42))-log(abs(s42-s41)))/log(0.5);
p=[p3,p4];






