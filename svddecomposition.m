function [U,Z,V] = svddecomposition(A)
%20211081037 2103 俞昊然
%   此函数用于对一个秩为r的mxn阶矩阵进行svd分解，其具有分解式A=UZV'，其中Z为2x2分块阵，左上角分块为一个对角阵
n=size(A);
m=n(1);n=n(2);
S=A'*A;
%分类讨论
if m>=n
%考虑m>=n的情况
    [V,Z]=eig(S);
    Z=sqrt(Z);
%对特征值大小和对应的特征向量进行冒泡排序
    for i=1:n
        for k=1:n-i
            if Z(i,i)<Z(k+i,k+i)
              s=Z(i,i);
              q=V(:,i);
              Z(i,i)=Z(k+i,k+i);
              V(:,i)=V(:,k+i);
              Z(k+i,k+i)=s;
              V(:,k+i)=q;
            end
        end
    end
    Z1=Z;Z1=inv(Z1);
    U=A*V*Z1;
%考虑m不等于n情况
    if m>n
    %扩充矩阵U的正交基，使用施密特正交化方法
       U=[U,zeros(m,m-n)];
       for  i=n+1:m
          U(i-n,i)=1;
          for j=1:i-1
              U(:,i)=U(:,i)-((U(:,i)'*U(:,j))/(U(:,j)'*U(:,j)))*U(:,j);
          end
          U(:,i)=U(:,i)/norm(U(:,i));
       end
       Z=[Z;zeros(m-n,n)];
    else
    end
else 
%考虑m<n  的情况
      S=A*A';
      [U,Z]=eig(S);
      Z=sqrt(Z);
      for i=1:m
        for k=1:m-i
            if Z(i,i)<Z(k+i,k+i)
              s=Z(i,i);
              q=U(:,i);
              Z(i,i)=Z(k+i,k+i);
              U(:,i)=U(:,k+i);
              Z(k+i,k+i)=s;
              U(:,k+i)=q;
            end
        end
      end
     r=rank(Z);
     V=zeros(n,r);
     for i=1:r
         V(:,i)=A'*U(:,i)/Z(i,i);
     end 
      %由于m<n还需要扩充V的基,使用施密特正交化方法
       V=[V,zeros(n,n-r)];
       for  i=r+1:n
          V(i-r,i)=1;
          for j=1:i-1
              V(:,i)=V(:,i)-((V(:,i)'*V(:,j))/(V(:,j)'*V(:,j)))*V(:,j);
          end
          V(:,i)=V(:,i)/norm(V(:,i));
       end
       Z=[Z,zeros(m,n-m)];
end
  
        
end

