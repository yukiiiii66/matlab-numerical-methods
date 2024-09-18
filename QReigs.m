function [B] = QReigs(A,eps,k)
%A为n阶方阵，eps收敛精度，k是最大迭代次数
%此为特征值计算的带原点qr方法
%先将A上Hessenberg化
[B] = hessenbergup(A);
H=B;
[n,n]=size(H);
m=n;
e=1;iter=0;
while m>1
    
  while(abs(H(m,m-1))>=eps)&&(iter<k)
      %求二阶矩阵的特征值
      T=eig_2(H(m-1:m,m-1:m));
      if isreal(T)
          %计算位移量
          [~,index]=min([abs(H(m,m)*[1;1]-T)]);
          %进行位移后再QR迭代
          [Q,R]=qr(H-T(index)*eye(m));
          %更新H
          H=R*Q+T(index)*eye(m);
      else
          %直接取H（m，m）为位移量
          [Q,R]=qr(H-H(m,m)*eye(m));
          H=R*Q+H(m,m)*eye(m);
      end
      iter=iter+1;
      if iter==10;
          disp('以右下角元为位移量所得的第十次迭代结果为')
          H
      end
      if iter==20;
          disp('以右下角元为位移量所得的第二十次迭代结果为')
          H
      end
  end
  if (abs(H(m,m-1))<eps)
      B(1:m,1:m)=H;
      m=m-1;
      %选取新的B的1：m-1子块
      H=B(1:m,1:m);
  end
  %上述是一个标准的选子块进行qr迭代的流程但是需要考虑出现复根的情况做如下的改动
  %考虑H右下角块复根情形
  if (m>2)&&(abs(H(m,m-1))>=eps)&&(abs(H(m-1,m-2))<eps)
      T=eig_2(H(m-1:m,m-1:m));%两块有复根的情况
      B(1:m,1:m)=H(1:m,1:m);
      B(m-1,m-1)=T(1);
      B(m,m)=T(2);
      m=m-2;
      H=B(1:m,1:m);%提取B的m-2阶子式
  end
  if (m==2)&&(abs(H(m,m-1))>=eps)%特殊情形只剩最后2阶
      T=eig_2(H(m-1:m,m-1:m));
      B(1:m,1:m)=H(1:m,1:m);
      B(m-1,m-1)=T(1);
      B(m,m)=T(2);
      if isreal(T(1))
          B(m,m-1)=0;
      end
      m=m-2;
  end          
end
B;
end

