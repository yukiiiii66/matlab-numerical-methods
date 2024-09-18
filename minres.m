function [x,N,I] = minres(A,m,n,e,k)
%A是系数矩阵，m为行，n为列，e是判断收敛的精度，k是设置的最大迭代次数
% 先进行初始化
x_star=ones(n,1);
b=A*x_star;
x=0.0005*rand(n,1)+ones(n,1);
r=b-A*x;
beta=norm(r,2);
g=zeros(n,1);
g(1)=beta;
V=zeros(n,k);
V(:,1)=r/beta;
 c=zeros(k,1);
 s=zeros(k,1);
 p=zeros(n,k);
%H=zeros(k+1,k);
for j=1:k
   
    %Lanczos过程
    H(j,j)=dot(V(:,j),A*V(:,j));
    if j==1
        w=A*V(:,j)-H(j,j)*V(:,j);
    else
        w=A*V(:,j)-H(j,j)*V(:,j)-H(j-1,j)*V(:,j-1);
    end
    H(j+1,j)=norm(w,2);
    H(j,j+1)=H(j+1,j);
    if j<k
      V(:,j+1)=w/H(j+1,j);
    end
    %Givens变换
    %由于Hk一般情况是三对角阵，需要注意的是前两轮循环 每一次做givens变换的要求不一样
    if j-2==-1
        %不做循环
        c(j)=H(j,j)/sqrt(H(j,j)^2+H(j+1,j)^2);
        s(j)=(H(j+1,j)/H(j,j))*c(j);
        H(j,j)=c(j)*H(j,j)+s(j)*H(j+1,j);
        H(j+1,j)=0;
        %计算方向pj
        p(:,1)=V(:,1)/H(1,1);
        
           
        
    elseif j-2==0   
        %做一次循环
            h=H(1,j);
            H(1,j)=c(1)*H(1,j)+s(1)*H(2,j);
            H(2,j)=-s(1)*h+c(1)*H(2,j);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        c(j)=H(j,j)/sqrt(H(j,j)^2+H(j+1,j)^2);
        s(j)=(H(j+1,j)/H(j,j))*c(j);
        H(j,j)=c(j)*H(j,j)+s(j)*H(j+1,j);
        H(j+1,j)=0;
        %计算方向pj
        p(:,2)=(V(:,2)-H(1,2)*p(:,1))/H(2,2);
    else
        %做两次循环
        for i=j-2:j-1
            h=H(i,j);
            H(i,j)=c(i)*H(i,j)+s(i)*H(i+1,j);
            H(i+1,j)=-s(i)*h+c(i)*H(i+1,j);
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%继续对下副对角线Givens变换
        c(j)=H(j,j)/sqrt(H(j,j)^2+H(j+1,j)^2);
        s(j)=(H(j+1,j)/H(j,j))*c(j);
        H(j,j)=c(j)*H(j,j)+s(j)*H(j+1,j);
        H(j+1,j)=0;
        %计算方向pj
        p(:,j)=(V(:,j)-H(j-2,j)*p(:,j-2)-H(j-1,j)*p(:,j-1))/H(j,j);
        
    end
    %对g进行更新
    g(j:j+1)=[c(j),s(j);-s(j),c(j)]*[g(j);0];    
    %更新x
    x=x+g(j)*p(:,j);
    %计算新的残量
    normg=norm(A*x-b,2);
    N(j)=normg;
    I=j;
    if normg<e
        break;
    end
    fprintf('第 %d 次迭代，使用MINRES算法计算的残量为 %f \n',j,normg);    
end
end

