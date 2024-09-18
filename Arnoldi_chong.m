function [V,hk] = Arnoldi_chong(A,r,k,eps)
%此程序为对Krylov子空间一组基底进行Arnoldi正交化的过程,考虑到我们处理的矩阵阶数很大，使用重正交化
%A一般是一个稀疏非对称非奇异的方阵
%A为指定的的线性方程组稀疏矩阵，r=Ax0-b为初始向量的残量，k是需要的正交化向量个数，最大为A的阶数
s=size(A);
n=s(1);
V=zeros(n,k);
%构造第一个正交基
v1=r/norm(r);
%对后续向量作Arnoldi过程
V(1:end,1)=v1;
for j=1:k-1
    w=A*V(1:end,j);
    h=zeros(j+1,1);
    for i=1:j
        h(i)=V(1:end,i)'*w;
        w=w-h(i)*V(1:end,i);
    end
    %使用重正交化方法
    for i=1:j
        s=V(1:end,i)'*w;
         h(i)=(i)+s;
         w=w-s*V(1:end,i);
    end
    h(j+1)=norm(w);
    hk(j)=h(j+1);
    %如果出现为0结果说明找到了一个不变子空间
    if h(j+1)<eps
        break;
    end
    V(1:end,j+1)=w./h(j+1);
end

end

