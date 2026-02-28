function [JR, bR] = rectangle_map(P1, P2, P3, P4)
 JR=zeros(2,2);
 bR=zeros(2,1);
 JR(1,1)=P2(1)-P1(1);
 JR(1,2)=P4(1)-P1(1);
 JR(2,1)=P2(2)-P1(2);
 JR(2,2)=P4(2)-P1(2);
 bR(1,1)=P1(1);
 bR(2,1)=P1(2);


