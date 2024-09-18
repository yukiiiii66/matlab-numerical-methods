function [A] = hessenbergup(A)
%A为n阶方阵对其上海森伯格化
%20211081037 2103 俞昊然
n=size(A,1);
for k=1:n-2
    [v,beta] = Householder2(A(k+1:n,k));
    H=eye(n-k)-beta*v*v';
    A(k+1:n,k:n)=H*A(k+1:n,k:n);
    A(1:n,k+1:n)=A(1:n,k+1:n)*H;
end

