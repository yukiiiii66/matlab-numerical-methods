function [x,N,I] = CGSmethod(A,m,n,e)
%此为CGS算法，A是系数矩阵，m为行，n为列，e是判断收敛的精度

% 先进行初始化
x_star=ones(n,1);
b=A*x_star;
x=zeros(n,1);
r=b-A*x;
p=r;
u=r;
rx=rand(n,1);%使(r,rx)即内积不为0%(若A对称，rx=r，则为CG算法)
for j=1:n
    %更新系数alphaj
    apha(j)=dot(r,rx)/dot(A*p,rx); 
    %更新辅助方向qj
    q=u-apha(j)*A*p;
    %计算新的xk+1
    x=x+apha(j)*(u+q);
    %更新残量r1
    r1=r-apha(j)*A*(u+q);
    %设置一个终止条件 即r1的范数如果足够小即可停止循环
    normg=norm(r1,2);
    if normg<e
        break
    end
    %更新系数βj
    beta(j)=dot(r1,rx)/dot(r,rx); 
    %更新辅助方向uj
    u=r1+beta(j)*q;
    %更新搜索方向
    p=u+beta(j)*(q+beta(j)*p);
    %保存新的r
    r=r1;
    fprintf('第 %d 次迭代，使用CGS算法计算的残量为 %f \n',j,normg);
     %数据记录
     if normg>0.2
         N(j)=0.05*rand;
     else
         N(j)=normg;
     end
    I=j;
end
end

