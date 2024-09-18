function [x,N] = Jacobi2(eps,A)
%20211081037-俞昊然
%A为第二题用数据做出的矩阵
s1=size(A);
s1=s1(1);
%设定初始点
b=A*zeros(s1,1);
x0=ones(s1,1);
x=x0;
iter=0;normg=norm(A*x-b)
while normg>eps && iter<40000
iter=iter+1; 
N(iter)=normg;
    for i=1:n^2
        s=A(i,1:end)*x1(1:end);
        s=s-A(i,i)*x1(i);
        %计算当前迭代中，第i个点的新迭代值
        x(i)=(b(i)-s)/A(i,i);
    end
   %存放当前迭代点
   x1=x;
    %计算当前迭代点残量
    normg=norm(A*x-b);
    
    fprintf('当前迭代次数为%d 当前迭代点的参量为 %d\n',iter,normg)
end
end