function [x,N] = Jacobi(eps,n)
%该函数为Jacobi迭代方法，eps为终止准则的检验精度
%20211081037-俞昊然
%生成方程组的系数矩阵，n为系数矩阵中T矩阵子块的阶数
[A,b] = Ematrix(n);
%设定初始点
x0=zeros(n^2,1);
x=x0;x1=x;
iter=0;normg=norm(A*x-b);N=[];
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
