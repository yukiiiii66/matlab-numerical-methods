function [xk] = GMRRS(eps,A)
%此方法为稀疏非对称的非奇异矩阵的广义最小残量法
%思路是先对Krylov子空间的基作单位正交化
%然后对教材p155页的一个方程组求最小二乘解yk进而由xk=x0+Vk*yk将当前子空间的最优解xk解出
s=size(A);
n=s(1);
xk=randn(n,1);
x0=xk;
%
b=A*zeros(n,1);
iter=1;
%先用残量表达式定义循环
r=A*xk-b;
normg=norm(r);
normg0=normg;
b0=norm(b);
%对krylov子空间的基底进行正交化
[V,hk] = Arnoldi_chong(A,r,300,1e-5);
rank(V)
pause(2)
s=zeros(n,1);
tao=zeros(n,1);

while normg>eps&&iter<300
    iter=iter+1;
    %计算Vk
    
    Vk=V(1:end,1:iter);
    Hk=Vk'*A*Vk;
   
    e2=zeros(1,iter);
    e2(iter)=hk(iter-1);
    Hk=[Hk;e2];
    
    %e1=zeros(iter+1,1);e1(1)=normg0;
    %最小二乘法使用QR分解正交化方法求解
    %Givens

   for i=1:iter
       
    
      if Hk(i+1,i)==0&&i~=1
        c=1;s(i)=0;
        tao(i)=(-1)^(i-1)*c;
        for j=1:i-1
            tao(i)=tao(i)*s(j);
        end
        return;
      end
      %givens
      
      if abs(Hk(i+1,i))>abs(Hk(i,i))
          t=Hk(i,i)/Hk(i+1,i);s(i)=1/sqrt(1+t^2);c=s(i)*t;
      else
        t=Hk(i+1,i)/Hk(i,i);c=1/sqrt(1+t^2);s(i)=c*t;
      end
      B=Hk(i,:);C=Hk(i+1,:);
      %QR分解
      Hk(i,:)=c*B+s(i)*C;
      Hk(i+1,:)=-s(i)*B+c*C;
      %计算tao
      if i==1
           tao(1)=b0*c;
      else
           tao(i)=(-1)^(i-1)*c;
        for j=1:i-1
            tao(i)=tao(i)*s(j);
        end
      end
   end
    Rk=Hk(1:iter,1:iter);
    y=inv(Rk)*tao(1:iter);
    
    %使用回代法求解 R是上三角
    %y=zeros(iter,1);
    %y(iter)=b0(iter)/R(iter,iter);
    %for i=iter-1:-1:1
     %   for k=i+1:iter
      %      y(i)=b0(i)-y(k)*R(i,k);
       % end
        %y(i)=y(i)/R(i,i);
    %end
    xk=x0+Vk*y;
    normg=norm(A*xk-b);
    fprintf('当前迭代次数为 %d 残量为 %f\n',iter,normg)
end 
   
    
end

