function [x,T] = LUcol_solve(L,U,P1,b)
%20211081037 2103 俞昊然
%通过对列选主元的LU分解来求解x
%b为列向量
tic;
%启动计时器
n=size(b);
n=n(1);
m=size(L);
m=m(1);
if m~=n
    fprintf('输入参数出现错误确保b是一个行向量')
else
    b=P1*b;
    y=zeros(n,1);
     %使用前代法求解y
    y(1)=b(1)/L(1,1);   
    for i=2:n
        for k=1:i-1
            b(i)=b(i)-y(k)*L(i,k);
        end
        y(i)=b(i)/L(i,i);    
    end
    x=zeros(n,1);
    %使用后代法求解x
    x(n)=y(n)/U(n,n);
    for i=n-1:-1:1
        for k=i+1:n
            y(i)=y(i)-x(k)*U(i,k);
        end
        x(i)=y(i)/U(i,i);
    end
    T = toc;
%终止计时器
end

