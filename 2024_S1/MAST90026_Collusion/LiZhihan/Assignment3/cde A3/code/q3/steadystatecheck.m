L=10;N=100;
[t1,solution1]=sol(L,N);
[t2,solution2]=sol2(L,N);
[t3,solution3]=sol3(L,N);
% Filter the results to include only t > 60
t_filtered1 = t1(t1 > 50);
solution_filtered1 = solution1(t1 > 60, :);
[m1,n1]=size(solution_filtered1);
error1=norm(solution_filtered1(:)-1,inf);
disp(error1)

t_filtered2 = t2(t2 > 50);
solution_filtered2 = solution2(t2 > 60, :);
[m2,n2]=size(solution_filtered2);
error2=norm(solution_filtered2(:)-1,inf);
disp(error2)

t_filtered3 = t3(t3 > 50);
solution_filtered3 = solution3(t3 > 60, :);
[m3,n3]=size(solution_filtered3);
error3=norm(solution_filtered3(:)-1,inf);
disp(error3)

