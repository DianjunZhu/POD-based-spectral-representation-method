function [ Um, a, b, h ] = POD_fft( Sn, ncut, dw, num)
% Input: Sn - EPSD matrix; ncut - truncation number of POD; dw - frequency
% resolution; num -  number of stochastic ground motions.
% Output: Um - stochastic ground motions; a - eigenvectors; b - eigenvalues;
% h - principal coordinates.
nw=size(Sn,1);% The number of frequency-domain samples
nt=size(Sn,2);% The number of time-domain samples
Um=zeros(nt,num);
G=sqrt(Sn);
R=G*G'/nt;% Covariance matrix
[a,b]=eigs(R,ncut);% Eigen decomposition
h=a'*G;
U=lhsnorm(zeros(nw,1),eye(nw), 2*num)';% Random number
e=U(:,1:num)+1i*U(:,num+1:2*num);
for i=1:num
    Y1=a*sqrt(2*dw).*e(:,i);
    Y2 = [Y1(1,:);Y1(2:end,:)/2;fliplr(conj(Y1(1:end,:)))/2];
    af=ifft(Y2, 'symmetric')*length(Y2);
    Um(:,i)=real(sum(h'.*af,2));
end
end