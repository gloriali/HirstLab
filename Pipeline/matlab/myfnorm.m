function d = myfnorm(X,Y,a)

ni=abs(a(3)); mi=a(1); si=abs(a(2));

fgauss=@(n,m,s,z)n/(s*sqrt(2*pi)) * exp(-(z-m).^2/(2*s^2) );
d = sum((Y - fgauss(ni,mi,si,X)).^2 );

end
