function [x,N,I] = IDR(A,m,n,e,s)
%此为IDR算法，A是系数矩阵，m为行，n为列，e是判断收敛的精度
%s是IDR算法独有的预设参数
% 先进行初始化
x_star=ones(n,1);
b=A*x_star;
x=zeros(n,1);
r=b-A*x;
%设置初始化矩阵
P=rand(n,s);
xk=zeros(n,s);rk=zeros(n,2*s+1);
delta_x=zeros(n,s);
delta_r=zeros(n,s);
for i=1:s
    %计算delta_x delta_r
    v=A*r;
    w=dot(v,r)/dot(v,v);
    delta_x(:,i)=w*r;
    delta_r(:,i)=-w*v;
    %计算x
    if i==1
        xk(:,i)=x+delta_x(:,i);
    else
        xk(:,i)=xk(:,i-1)+delta_x(:,i);
    end
    %计算r
    if i==1
        rk(:,i)=r+delta_r(:,i);
    else
        rk(:,i)=rk(:,i-1)+delta_r(:,i);
    end
end

normg=norm(rk(:,s));
j=0;
while normg>e
    j=j+1;
    %更新delta_r与delta_x的算法
    for k=0:1:s
      c=(P'*delta_r)\(P'*rk(:,i));
      v=rk(:,i)-delta_r*c;
      if k==0
          t=A*v;
          w=dot(t,v)/dot(t,t);
          %delta_r(:,i)=-delta_r*c-w*t;
          %delta_x(:,i)=-delta_x*c+w*v;
          delta_rk=-delta_r*c-w*t;
          delta_xk=-delta_x*c+w*v;
      else
          %delta_x(:,i)=-delta_x*c+w*v;
          %delta_r(:,i)=-A*delta_x(:,i);
          delta_xk=-delta_x*c+w*v;
          delta_rk=-A*delta_xk;
      end
      rk(:,i+1)=rk(:,i)+delta_rk;
      xk(:,i+1)=xk(:,i)+delta_xk;
      i=i+1;
      %更新delta_r
      if k==0
         delta_r=[delta_r(:,1:s-1),delta_rk];
      else
          delta_r=[delta_r(:,2:end),delta_rk];
      end
      %更新delta_x
      if i==s+1
          
          delta_x=[delta_x(:,1:s-1),delta_xk];
      else
          delta_x=[delta_x(:,2:end),delta_xk];
      end
    end
    rr=rk(:,end);
    normg=norm(rr);
    fprintf('第 %d 次迭代，使用IDR算法计算的残量为 %f \n',j,normg);
     %数据记录
    N(j)=normg;
    I=j;
end
x=xk;
end

