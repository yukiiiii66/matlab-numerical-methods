function  [Q,R] = HESSEN_QR(A)
%20211081037 2103 俞昊然
%此为专门针对海森伯格矩阵的QR分解，使用的是Givens变换
%首先验证A是否已经海森伯格化，并确定是不是非退化情形
[n,~]=size(A);
for i=1:n-2
    for j=i+2:n
    if abs(A(j,i))>1e-3
        error('你所输入的不是上海森伯格矩阵')
    end
    end
    if A(i+1,i)==0
        error('你所输入的是退化情形')
    end
end
%进行givens变换以QR分解
Q=eye(n);
for k=1:n-1
    sigma=sqrt(A(k,k)^2+A(k+1,k)^2);
    if abs(A(k+1,k))>=abs(A(k,k))
        t=A(k,k)/A(k+1,k);
        s=1/sqrt(1+t^2);
        c=s*t;
    else
       t=A(k+1,k)/A(k,k);
        c=1/sqrt(1+t^2);
        s=c*t;
    end
    D= A(k,:);
    A(k,:)=c*A(k,:)+s*A(k+1,:);
     A(k+1,:)=-s*D+c*A(k+1,:);
    H=eye(n);H(k,k)=c;H(k+1,k+1)=c;H(k,k+1)=s;H(k+1,k)=s;
    Q=Q*H;
end
R=A;

end

