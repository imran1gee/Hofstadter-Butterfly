clear all
tic
qmax = 40;
for i = 1:1:qmax                                             
    for j = 1:1:qmax                                         
        A(i,j) = i/j;                                        
    end                                                      
end 
B(:) = unique(reshape(A',1,[]));                             
for i = 1:1:length(B)                                        
    if B(i)<1
        C1(i) = B(i);
    end
end                                                          
C2 = nonzeros(C1');                                          
[p,q1] = numden(sym(C2));                                                                           
Q(:) = double(q1(:));                                        
P(:) = double(p(:));  
l = 19;


for i = 1:1:length(Q)
    D = Var(P(i),Q(i),l);
    D = D(:)/norm(D);
    alpha(i) = (P(i)/Q(i));
    DE(i,:) = D(:);
    clear D
end
N1 = -1:2/(length(DE(1,:))-1):1;
N = N1/max(N1);
contourf(alpha,N,DE')
A1{1} = DE;
A1{2} = alpha;
A1{3} = N;
save('testG1.mat','A1')
toc



function D = Var(p,q,l)
        x = 0:(1/l)*2*pi/q:2*pi/q;
        y = 0:(1/l)*2*pi:2*pi;
        for i = 1:1:length(x)
            for j = 1:1:length(y)
                E(i,j,:) = eig(Ham(x(i),y(j),p,q));
            end
        end
        EF1 = sort(reshape(E,[],1));
        LL = (3/2)*min(min(min(E(:,:,:))));
        UL = (3/2)*max(max(max(E(:,:,:))));
        st1 = LL:(UL-LL)/(length(EF1)-1):UL;
        add = zeros(1,length(EF1)+length(st1));
        for i = 1:1:length(st1)
            add(i) = EF1(i);
            add(i+length(st1)) = st1(i);
        end
        add = sort(add);
        g = 0.01;
        count = 0;
        for i = 1:2*q:length(add)
            count =  count + 1;
            D(count) = (1/1)*(g/l^2)*sum(1./((add(i)-EF1(:)).^2+g^2));
        end
end


function H = Ham(x,y,p,q)
         H = zeros(2*q,2*q);
         for j = 1:2:2*q-1
             H(j,j+1) = (exp(1i*y)+exp(-1i*(2*pi*(p/q)*(j+1)/2)));
         end
         for j = 2:2:2*q-1
             H(j,j+1) = exp(-1i*(2*pi*(p/q)*(j+2)/2))*exp(-1i*x);
         end
         H(1,2*q) = H(1,2*q)+exp(1i*x);
         H = H+ctranspose(H);
end