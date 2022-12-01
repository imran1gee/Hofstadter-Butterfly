clear all
p=1;
q=3;
m = 2;
alpha = p/q;
pts=20;

dk = 2*pi/(q*pts);
kx = 0:dk:2*pi/q;
ky = 0:dk:2*pi;
kxx = dk:dk:2*pi/q;
kyy = dk:dk:2*pi;
Eig = zeros(length(kx),length(ky),2*q);
Vec = zeros(length(kx),length(ky),2*q,2*q);
U1 = zeros(length(kx)-1,length(ky)-1);
U2 = U1;
U3 = U1;
U4 = U1;
UU1 = U1;
UU2 = U1;
UU3 = U1;
UU4 = U1;
% BF = U1;

for i = 1:1:length(kx)
    for j = 1:1:length(ky)
        [V, D] = eig(Ham(kx(i),ky(j),alpha,q));
        for ii = 1:1:2*q
            Eig(i,j,ii) = D(ii,ii);
            Vec(i,j,ii,:) = V(:,ii);
            nom(i,j,ii) = det(ctranspose(squeeze(Vec(i,j,ii:q,:)))*squeeze(Vec(i,j,ii:q,:)));
        end
    end
end

% for i = 1:1:2*q
%     mesh(kx,ky,Eig(:,:,i)')
%     hold on
% end
% hold off

% m = 3;

for i = 1:1:length(kx)-1
    for j = 1:1:length(ky)-1
        UU1(i,j) =(ctranspose(squeeze(Vec(i,j,m,:)))*squeeze(Vec(i+1,j,m,:)));
        if UU1(i,j)~=0
           U1(i,j) = ((UU1(i,j))/norm(UU1(i,j)));
        end
        UU2(i,j) =(ctranspose(squeeze(Vec(i+1,j,m,:)))*squeeze(Vec(i+1,j+1,m,:)));
        if UU2(i,j)~=0
            U2(i,j) = ((UU2(i,j))/norm(UU2(i,j)));
        end
        UU3(i,j) =(ctranspose(squeeze(Vec(i+1,j+1,m,:)))*squeeze(Vec(i,j+1,m,:)));
        if UU3(i,j)~=0
           U3(i,j) = ((UU3(i,j))/norm(UU3(i,j)));
        end
        UU4(i,j) =(ctranspose(squeeze(Vec(i,j+1,m,:)))*squeeze(Vec(i,j,m,:)));
        if UU4(i,j)~=0
            U4(i,j) = ((UU4(i,j))/norm(UU4(i,j)));
        end
            BF(i,j) =(imag(log(U1(i,j)*U2(i,j)*U3(i,j)*U4(i,j))));
    end
end
figure 
contourf(kx(1:end-1),ky(1:end-1),BF')

Chern =(sum(sum(BF(:,:))))/(2*pi);


Chern


function Ha = Ham(x,y,alpha,q)
    Hf = zeros(2*q,2*q);
    for j = 1:2:2*q-1
%         Hf(j,j+1) = (exp(1i*(4*pi/3*(alpha)*(j+1)/2+y))+exp(-1i*(2*pi/3*(alpha)*(j+1)/2)));
        Hf(j,j+1) = (exp(1i*(y))+exp(-1i*(2*pi*(alpha)*(j+1)/2)));
    end
    for j1 = 2:2:2*q-2
%         Hf(j1,j1+1) = exp(-1i*(2*pi/3*(alpha)*(j1+2)/2));
        Hf(j1,j1+1) = exp(-1i*(2*pi*(alpha)*(j1+2)/2))*exp(-1i*(x));
    end
%     Hf(1,2*q) = exp(-1i*(2*pi/3*(alpha)*(q)-q*x));
    Hf(1,2*q) = Hf(1,2*q)+exp(1i*(x));
    Ha = Hf+ctranspose(Hf);
end