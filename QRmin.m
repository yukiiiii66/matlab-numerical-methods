function [x,A,r] = QRmin(m,n,b)
%20211081037 2103 俞昊然
%本方法为QR正交化方法求最小二乘问题
%m，n为随机矩阵的行数与列数，b是常向量，且b维数=m
if m<n
    error('你所输入的问题一定有解，不适合寻找最小二乘解');
end
%随机生成矩阵A
A=randn(m,n);
A1=A;
r=rank(A);
if r<n
    fprintf('你所输入的矩阵为秩亏损型');
else
    %对A进行QR分解 
    for i=1:n
        %调用householder变换
        [v,beta]=Householder1(A(i:m,i));
        A(i:m,i:n)=A(i:m,i:n)-(beta*v)*(v'*A(i:m,i:n));
        d(i)=beta;
        
        A(i+1:m,i)=v(2:m-i+1);
        
    end
    R1=A(1:n,1:n);
    %提取Q,R矩阵
    Q=eye(m,m);
    for i=1:n
        v=[1; A(i+1:m,i)];
        H=eye(m+1-i)-d(i)*v*v';
        H=[eye(i-1),zeros(i-1,m+1-i);zeros(m+1-i,i-1),H];
        Q=Q*H;
    end
    Q1=Q(:,1:n);
    b1=b;
    b=Q1'*b;
    %后代法求解
    x=zeros(n,1);
    x(n)=b(n)/R1(n,n);
    for i=n-1:-1:1
        for k=i+1:n
            b(i)=b(i)-x(k)*R1(i,k);
        end
        x(i)=b(i)/R1(i,i);
    end
    %计算残量
    A=A1;
    r=norm(A'*(A*x-b1));
    
    
end