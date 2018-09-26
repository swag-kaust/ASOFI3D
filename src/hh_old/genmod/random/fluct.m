function v=fluct(nx,nz,a,dh,vm,sigma)

x=(-nx/2:nx/2)*dh; z=(-nz/2:nz/2)*dh;

for k=1:nx,
   for l=1:nz,
		r=sqrt(x(k)*x(k)+z(l)*z(l));
      A(k,l)=sigma*sigma*exp(-r/a);  % exponential
   end
end


AMP=sqrt(abs(fft2(A)));

rand('seed',0);
AMP=sqrt(nx*nz*dh*dh)*AMP.*exp(-i*pi*(2*rand(size(AMP))-1));

d=real(ifft2(AMP));

v=vm+d;

break

