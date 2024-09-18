function [namb,x,iter] = mifa(A,u,eps,k)
%20211081037 俞昊然 数计2103
%此为幂法求解最大特征值的程序，其中namb为最大特征值，x为其对应的特征向量，我们要求A的最大特征值是半单的
%A为求特征值的矩阵，u为初始向量，eps为收敛精度，k为最大迭代步数
%设定初始值
e=1;
iter=0;
[~,index]=max(abs(u));
namb1=u(index);
u=u./namb1;
%开始迭代求解
while e>eps&&iter<k
    v=A*u;
    [~,index]=max(abs(v));
    namb=v(index);
    u=v./namb;
    %计算残量，并更新namb1，更新步数
    e=abs(namb1-namb);namb1=namb;iter=iter+1;
end
x=u;
end