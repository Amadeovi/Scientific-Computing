function [u,r]=BSSOR(A,F,w,tol)
n=length(F);
U=[]; 
R(:,1)=F; U(:,1)=zeros(length(F),1);
k=1; r(1)=1; phi=[];
% Block SSOR
% p=1 so only one big block
while norm(R(:,k))/norm(F)>tol
    k=k+1;
    % forward step
    U(:,k)=F(:,1)-A*U(:,k-1);
    U(:,k)=(1-w)*U(:,k-1)+w*U(:,k);
    phi(:,k)=U(:,k);
    % backward step
    U(:,k)=F(:,1)-A*phi(:,k);
    U(:,k)=(1-w)*phi(:,k)+w*U(:,k);
    R(:,k)=F-A*U(:,k);
    r(k)=norm(R(:,k))/norm(F);
end

% p=2

    U(1:,k)


scatter(1:k,r,'b')
set(gca,'yscale','log')
hold on
semilogy(1:k,r,'b')
ylabel('$ \frac{\|r_{k}\|_{2}}{\|F\|_{2}}$','Interpreter','Latex')
xlabel('Number of iterations')
u=U(:,end);





    % forward step
    U(:,k)=F(:,1)-A*U(:,k-1);
    U(:,k)=(1-w)*U(:,k-1)+w*U(:,k);
    phi(:,k)=U(:,k);
    % backward step
    U(:,k)=F(:,1)-A*phi(:,k);
    U(:,k)=(1-w)*phi(:,k)+w*U(:,k);
    R(:,k)=F-A*U(:,k);
    r(k)=norm(R(:,k))/norm(F);

end


