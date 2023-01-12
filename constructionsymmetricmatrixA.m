n=input('Number of subintervals=');
d=input('dimension of the problem');
h=1/n; x=0:h:1; y=0:h:1; z=0:h:1;
% the objetive is to obtain the symmetric matrix A and vector F
% we are going to use the Kroenecker product for carrying out this task
if d==2
   tic
   That=(1/h^2)*(diag([h^2,4*ones(1,n-1),h^2])+diag([0,-1*ones(1,n-2),0],1)+diag([0,-1*ones(1,n-2),0],-1));
   Ihat=-(1/h^2)*diag([0,ones(1,n-1),0]);
   I=eye(n+1);
   % for defining the blocks A_{1,1} and A_{n+1,n+1}
   A=kron(diag([1,zeros(1,n-1),1]),I);
   % for defining the blocks A_{i,i} of the main diagonal i=2,...,n
   A=A+kron(diag([0,ones(1,n-1),0]),That);
   % for defining the blocks A_{i,i-1} of the 1st low diagonal i=2,...,n
   A=A+kron(diag([0,ones(1,n-2),0],-1),Ihat);
   % for defining the blocks A_{i,i+1} of the 1st upper diagonal i=2,...,n
   A=A+kron(diag([0,ones(1,n-2),0],1),Ihat);
   F=[];
   %interior points
   for i=2:n
       for j=2:n
           F(i+(j-1)*(n+1))=(x(i)^2+y(j)^2)*sin(x(i)*y(j));
       end
   end
   % boundaries conditions
   for i=1:n+1
       F(i)=sin(x(i)*y(1));
       F(i+n*(n+1))=sin(x(i)*y(n+1));
   end
   for j=1:n+1
       F(1+(j-1)*(n+1))=sin(x(1)*y(j));
       F(j*(n+1))=sin(x(n+1)*y(j));
   end
   % contribution of the boundaries conditions to the internal points
   for i=2:n
       F(i+(n+1))=F(i+(n+1))+(1/h^2)*sin(x(i)*y(1));
       F(i+(n-1)*(n+1))=F(i+(n-1)*(n+1))+(1/h^2)*sin(x(i)*y(n+1));
   end
   for j=2:n
       F(2+(j-1)*(n+1))=F(2+(j-1)*(n+1))+(1/h^2)*sin(x(1)*y(j));
       F(n+(j-1)*(n+1))=F(n+(j-1)*(n+1))+(1/h^2)*sin(x(n+1)*y(j));
   end
   F=F'; % vertical vector.
 %the exact solution will be
   for i=1:n+1
       for j=1:n+1
           uex(i+(j-1)*(n+1),1)=sin(x(i)*y(j));
       end
   end
   toc
elseif d==3
    tic
    T3=diag([h^2,6*ones(1,n-1),h^2])+diag([0,-1*ones(1,n-2),0],1)+diag([0,-1*ones(1,n-2),0],-1);
    Ihat=-diag([0,ones(1,n-1),0]);
    % to define the blocks corresponding to k=1 and k=n+1
    A=kron(kron(diag([1,zeros(1,n-1),1]),eye(n+1)),eye(n+1));
    % to define the blocks of the diagonal
    B=kron(diag([1,zeros(1,n-1),1]),eye(n+1))+(1/h^2)*kron(diag([0,ones(1,n-1),0]),T3)+(1/h^2)*kron(diag([0,ones(1,n-2),0],-1),Ihat)+(1/h^2)*kron(diag([0,ones(1,n-2),0],1),Ihat);
    % to define the diagonal itself (expect k=1 and k=n+1)
    A=A+kron(diag([0,ones(1,n-1),0]),B);
    % to define the contributions of u_{i,j,k-1} and u_{i,j,k+1}
    C=(1/h^2)*kron(diag([0,ones(1,n-1),0]),Ihat);
    % to define the blocks above and below the diagonal
    A=A+kron(diag([0,ones(1,n-2),0],-1),C)+kron(diag([0,ones(1,n-2),0],1),C);
    %interior points
    F=[];
    for i=2:n
        for j=2:n
            for k=2:n
                F(i+(j-1)*(n+1)+(k-1)*(n+1)^2,1)=(x(i)^2*y(j)^2+x(i)^2*z(k)^2+y(j)^2*z(k)^2)*sin(x(i)*y(j)*z(k));
            end
        end
    end
    % boundary points
    for i=1:n+1
        for j=1:n+1
            %k=1,k=n+1
            F(i+(j-1)*(n+1))=sin(x(i)*y(j)*z(1));
            F(i+(j-1)*(n+1)+n*(n+1)^2)=sin(x(i)*y(j)*z(n+1));
        end
    end
    for i=1:n+1
        for k=1:n+1
            %j=1,j=n+1
            F(i+(k-1)*(n+1)^2)=sin(x(i)*y(1)*z(k));
            F(i+n*(n+1)+(k-1)*(n+1)^2)=sin(x(i)*y(n+1)*z(k));
        end
    end
    for j=1:n+1
        for k=1:n+1
            %i=1,i=n+1
            F(1+(j-1)*(n+1)+(k-1)*(n+1)^2)=sin(x(1)*y(j)*z(k));
            F(n+1+(j-1)*(n+1)+(k-1)*(n+1)^2)=sin(x(n+1)*y(j)*z(k));
        end
    end
    % contributions of the boundary points for the internal points
    for i=2:n
        for j=2:n
            % k=1,k=n+1 contribution to k=2, k=n
            F(i+(j-1)*(n+1)+(n+1)^2)=F(i+(j-1)*(n+1)+(n+1)^2)+(1/h^2)*sin(x(i)*y(j)*z(1));
            F(i+(j-1)*(n+1)+(n-1)*(n+1)^2)=F(i+(j-1)*(n+1)+(n-1)*(n+1)^2)+(1/h^2)*sin(x(i)*y(j)*z(n+1));
        end
    end
    for i=2:n
        for k=2:n
            % j=1,j=n+1 contribution to j=2 j=n
            F(i+(n+1)+(k-1)*(n+1)^2)=F(i+(n+1)+(k-1)*(n+1)^2)+(1/h^2)*sin(x(i)*y(1)*z(k));
            F(i+(n-1)*(n+1)+(k-1)*(n+1)^2)=F(i+(n-1)*(n+1)+(k-1)*(n+1)^2)+(1/h^2)*sin(x(i)*y(n+1)*z(k));
        end
    end
    for j=2:n
        for k=2:n
            %i=1,i=n+1 contribution to i=2 i=n
            F(2+(j-1)*(n+1)+(k-1)*(n+1)^2)=F(2+(j-1)*(n+1)+(k-1)*(n+1)^2)+(1/h^2)*sin(x(1)*y(j)*z(k));
            F(n+(j-1)*(n+1)+(k-1)*(n+1)^2)=F(n+(j-1)*(n+1)+(k-1)*(n+1)^2)+(1/h^2)*sin(x(n+1)*y(j)*z(k));
        end
    end
    %the exact solution will be
   for i=1:n+1
       for j=1:n+1
           for k=1:n+1
               uex(i+(j-1)*(n+1)+(k-1)*(n+1)^2,1)=sin(x(i)*y(j)*z(k));
           end
       end
   end
   toc
else
    disp('Dimension must be d=2 or d=3')
end

