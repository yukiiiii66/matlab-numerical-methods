function [x,N] = Gauss(eps,n)
%该函数为Gauss迭代方法，eps为终止准则的检验精度
%20211081037-俞昊然
%生成方程组的系数矩阵，n为系数矩阵中T矩阵子块的阶数
[A,b] = Ematrix(n);
%设定初始点
x0=zeros(n^2,1);
x=x0;
iter=0;normg=norm(A*x-b);
N=[];
while normg>eps && iter<40000
iter=iter+1;
 N(iter)=normg;
    for i=1:n^2
        s=0;d=0;
        %考虑在当前迭代次数中，已生成的新迭代点
        if i>1
                s=A(i,1:i-1)*x(1:i-1);
        else
        end
        %考虑在当前迭代次数中，还未生成的迭代点
        if i<n^2
                d=A(i,i+1:end)*x(i+1:end);
        else
        end
        %计算当前迭代中，第i个点的新迭代值
        x(i)=(b(i)-(s+d))/A(i,i);
    end
    %计算当前迭代点残量
    normg=norm(A*x-b);
    fprintf('当前迭代次数为%d 当前迭代点的参量为 %d\n',iter,normg)
end
end