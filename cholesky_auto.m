function [x] = cholesky_auto(A,b)
%20211081037 2103 俞昊然
%此为Cholesky分解求解方程组自动程序，方程组系数矩阵要求对称正定
n=size(A);
n=n(1);
L=zeros(n);
L(1,1)=sqrt(A(1,1));
%计算第一列下方元素
for i=2:n
    L(i,1)=A(i,1)/L(1,1);
end
for i=2:n-1
    %计算对角元素
    L(i,i)=sqrt(A(i,i)-sum((L(i,1:i-1)).*L(i,1:i-1)));
    %计算各列对角元下方元素
    for j=i+1:n
        L(j,i)=A(j,i);
        for k=1:i-1
        L(j,i)=L(j,i)-L(i,k)*L(j,k);
        end
        L(j,i)=L(j,i)/L(i,i);
    end
end
L(n,n)=sqrt(A(n,n)-sum((L(n,1:n-1)).*L(n,1:n-1)));
y=zeros(n,1);
    %前代法解y
    y(1)=b(1)/L(1,1);
    for i=2:n
        for k=1:i-1
            b(i)=b(i)-y(k)*L(i,k);
        end
        y(i)=b(i)/L(i,i);    
    end
    x=zeros(n,1);
    %后代法解x
    L=L';
    x(n)=y(n)/L(n,n);
    for i=n-1:-1:1
        for k=i+1:n
            y(i)=y(i)-x(k)*L(i,k);
        end
        x(i)=y(i)/L(i,i);
    end

end

