function ERROR=er(t,N,K)
ERROR=zeros(1,2);
[~,ERROR(1)]=CN(t,N,K);
[~,ERROR(2)]=TRBFD(t,N,K);