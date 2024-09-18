function [x,N,I] = GPBICG(A,m,n,e)
%此为GPBICG算法，A是系数矩阵，m为行，n为列，e是判断收敛的精度

% 先进行初始化
x_star=ones(n,1);
b=A*x_star;
x=zeros(n,1);
r=b-A*x;
%设置rx,使(r,rx)即内积不为0即可
%(若A对称，rx=r，则为CG算法)
rx=rand(n,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=zeros(n,1);
w=zeros(n,1);
beta=0;
for j=0:1:n
    %更新p
    if j==0
        p=r;
    else
        p=r+beta*(p-u);
    end
    %更新alpha
    alpha=dot(r,rx)/dot(A*p,rx);
    %更新向量y
    y=t-r-alpha*w+alpha*A*p;
    %更新向量t
    t0=t; %存储上一个t向量
    t=r-alpha*A*p;
    %计算系数ζ，η
    if j==0
        zeta=dot(t,A*t)/dot(A*t,A*t);
        eta=0;
    else
        zeta=((y'*y)*dot(t,A*t)-(y'*t)*dot(y,A*t))/(dot(A*t,A*t)*(y'*y)-(dot(y',A*t))^2);
        eta=(dot(A*t,A*t)*(y'*t)-dot(y,A*t)*dot(t,A*t))/(dot(A*t,A*t)*(y'*y)-(dot(y',A*t))^2);
    end
    %计算向量u
    if j==0
        u=zeta*A*p;
    else
        u=zeta*A*p+eta*(t0-r+beta*u);
    end
    %计算向量v
    if j==0
        v=zeta*r;
    else
        v=zeta*r+eta*v-alpha*u;
    end
    %更新x
    x=x+alpha*p+v;
    r0=r;%存储旧的r
    %更新残量
    r=t-eta*y-zeta*A*t;
    normg=norm(r);
    %更新beta
    beta=(dot(rx,r)/dot(rx,r0))*(alpha/zeta);
    %更新小Ω
    w=A*t+beta*A*p;
    fprintf('第 %d 次迭代，使用GPBICG算法计算的残量为 %f \n',j,normg);
    %判断条件
     %数据记录
    N(j+1)=normg;
    I=j+1;
    if normg<e
        break;
    end
      
end
end

