function [u,r] = SSORpreconditionerCG(A,F,Mprime,tol)
% This function uses the SSOR BIM as a preconditioner for 
% the Conjugate Gradient method.
n=length(F);
% Column 1 correspond to iteration 0
U(:,1)=zeros(n,1); R(:,1)=F;
k=1; Z=[]; P=[]; beta=[]; alpha=[];
while (norm(R(:,k))/norm(F))>tol
    Z(:,k)=Mprime*R(:,k);
    if k==1
        P(:,2)=Z(:,1); %P=[0,p1,p2,..]
    else
        beta(k+1)=(R(:,k)'*Z(:,k))/(R(:,k-1)'*Z(:,k-1));
        P(:,k+1)=Z(:,k)+beta(k+1)*P(:,k);
    end
    alpha(k+1)=(R(:,k)'*Z(:,k))/(P(:,k+1)'*A*P(:,k+1));
    U(:,k+1)=U(:,k)+alpha(k+1)*P(:,k+1);
    R(:,k+1)=R(:,k)-alpha(k+1)*A*P(:,k+1);
    r(k+1)=norm(R(:,k+1))/norm(F);
    k=k+1;
end
r=r';
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
    





