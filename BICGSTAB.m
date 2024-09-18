function [x,N,I] = BICGSTAB(A,m,n,e)
%此为BICGSTAB算法，A是系数矩阵，m为行，n为列，e是判断收敛的精度
%先设置初始变量
x_star=ones(n,1);
b=A*x_star;
x=zeros(n,1);
r=b-A*x;
%设置rx,使(r,rx)即内积不为0即可
%(若A对称，rx=r，则为CG算法)
rx=rand(n,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=0;
for j=0:1:100000
    %更新p
    if j==0
        p=r;
    else
        p=r+beta*(p-w*A*p);
    end
    %更新alpha
    alpha=dot(r,rx)/dot(A*p,rx);
    %更新向量t
    t=r-alpha*A*p;
    %计算w
    w=t'*A*t/dot((A*t),(A*t));
    %更新x
    x=x+alpha*p+w*t;
    r0=r;%存储旧的r
    %更新残量
    r=t-w*A*t;
    normg=norm(r);
    fprintf('第 %d 次迭代，使用BICGSTAB算法计算的残量为 %f \n',j,normg);
    %数据记录
    N(j+1)=normg;
    I=j+1;
    %判断条件
    if normg<e
        break;
    end
    beta=dot(r,rx)/dot(rx,r0)*alpha/w;
    
end
end

