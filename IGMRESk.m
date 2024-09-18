function [x,N,iter] = IGMRESk(A,m,n,e,k)
%此为重启的广义最小残量法，A是系数矩阵，m为行，n为列，e是判断收敛的精度，k是重启动参数
%先设置初始变量
x_star=ones(n,1);
b=A*x_star;
x=zeros(n,1);
r=b-A*x;
beta=norm(r,2);
V=zeros(n,k);
V(:,1)=r/beta;
H=zeros(k+1,k);
normg=1;iter=0;

while normg>e
    iter=iter+1;
    r=b-A*x;
    beta=norm(r,2);
    V(:,1)=r/beta;
%进行Arnoldi过程
    for j=1:k
      for i=1:j
        H(i,j)=dot(V(:,i),A*V(:,j));
      end
      w=A*V(:,j);
      for i=1:j
        w=w-H(i,j)*V(:,i);
      end
      H(j+1,j)=norm(w,2);
      if H(j+1,j)==0
        k=j;
        break;
      end
      if j<k
         V(:,j+1)=w/H(j+1,j);
      end
    end
%结束本轮Arnoldi过程
  e1=zeros(k+1,1);
  e1(1)=beta;
  %最小二乘问题求解计算线性组合系数向量
  y=(H'*H)\(H'*e1);
  %更新xk
  xk=x+V*y;
  normg=norm(A*xk-b);
  fprintf('第 %d 次迭代，使用GMRES(m)算法计算的残量为 %f \n',iter,normg);
  N(iter)=normg;
  x=xk;
end
end

