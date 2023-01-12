function [u,r,mu]=SSOR(A,F,w,tol)
tic
% This function uses the symmetric successive overrelaxation method
% to solve the system A*u=f
% it also gives us the vector r to measure the errors
% the asymptotic rate of convergencem mu
% and the CPU time required to solve the linear system using SSOR as a BIM.
n=length(F);
U=[]; %matrix colum k represent iteration vector u^{k+1}
U(:,1)=zeros(length(F),1); %zero initial guess for u^{h}.
R(:,1)=F; %matrix -> each column 
k=1;
r(1)=1;
phi=[];
while norm(R(:,k))/norm(F)>tol
    k=k+1;
    % the forward part is:
    U(1,k)=F(1,1)-A(1,2:n)*U(2:n,k-1);
    U(1,k)=(1-w)*U(1,k-1)+w*U(1,k); 
    for i=2:length(F)
        % the damping is performed component wise as well    
        U(i,k)=(1-w)*U(i,k-1)+w*(F(i,1)-A(i,1:i-1)*U(1:i-1,k)-A(i,i+1:n)*U(i+1:n,k-1))/A(i,i); 
    end
    % the backward part is:
    phi(:,k)=U(:,k);
    U(n,k)=F(n,1)-A(n,1:n-1)*phi(1:n-1,k);
    U(n,k)=(1-w)*phi(n,k)+w*U(n,k);
    for i=n-1:-1:1
        % the damping is performed component wise as well
        U(i,k)=(1-w)*phi(i,k)+w*(F(i,1)-A(i,1:i-1)*phi(1:i-1,k)-A(i,i+1:n)*U(i+1:n,k))/A(i,i);
    end
    R(:,k)=F-A*U(:,k);
    r(k)=norm(R(:,k))/norm(F);
end
r=r';
% to compute the  asymptotic rate of convergence we proceed as follows:
counter=0;
for i=length(r)-5:length(r)-1
    counter=counter+(r(i+1)/r(i));
end
mu=(counter/5);
scatter(1:k,r,'b')
set(gca,'yscale','log')
hold on
semilogy(1:k,r,'b')
ylabel('$ \frac{\|r_{k}\|_{2}}{\|F\|_{2}}$','Interpreter','Latex')
xlabel('Number of iterations')
grid
u=U(:,end);
toc
end

