function [x,A,r] = SVDmin(m,n,b)
%20211081037 2103 俞昊然
%本方法为奇异值方法求最小二乘问题
%m，n为随机矩阵的行数与列数，b是常向量，且b维数=m
if m<n
    error('你所输入的问题一定有解，不适合寻找最小二乘解');
end
%随机生成矩阵A
A=randn(m,n);
r=rank(A);
if r<n
    fprintf('你所输入的矩阵为秩亏损型')'
else
    [U,Z,V] = svddecomposition(A);
    B=U'*b;
    B1=B(1:n);
    Z=Z(1:n,1:n);
    x=V*inv(Z)*B1;
     %计算残量
    r=norm((A'*(A*x-b)));
end