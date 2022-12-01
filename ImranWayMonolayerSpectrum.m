clear all
p=1;
q=2;
alpha = p/q;
pts=50;

dk = 2*pi/(q*pts);
kx = 0:dk:2*pi/q;
ky = 0:dk:2*pi;
Eig = zeros(length(kx),length(ky),2*q);
tic
for i = 1:1:length(kx)
    for j = 1:1:length(ky)
        Eig(i,j,:) = eig(Ham(kx(i),ky(j),alpha,q));
    end
end
for i = 1:1:2*q
    mesh(kx,ky,Eig(:,:,i)')
    EG(i,j) =  -max(max(Eig(:,:,i)))+min(min(Eig(:,:,i+1)))
    hold on
end
box on
xlabel('kxa')
ylabel('kya')
zlabel('E[meV]')
hold off
toc


function Ha = Ham(x,y,alpha,q)
    Hf = zeros(2*q,2*q);
    for j = 1:2:2*q-1
        Hf(j,j+1) = (exp(1i*y)+exp(-1i*(2*pi*(alpha)*(j+1)/2)));
    end
    for j1 = 2:2:2*q-2
        Hf(j1,j1+1) = exp(-1i*(2*pi*(alpha)*(j1+2)/2))*exp(-1i*x);
    end
    Hf(1,2*q) = Hf(1,2*q)+exp(1i*(x));
    Ha = Hf+ctranspose(Hf);
end