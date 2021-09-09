load ZrsatDif_2D
%load M3


timesteps=12
IMAX=250;
JMAX=300;
KMAX=250;




X=50000/IMAX:50000/IMAX:50000;
Y=60000/JMAX:60000/JMAX:60000;
Z=50000/IMAX:50000/IMAX:50000;
DX=50000/IMAX;
DY=60000/JMAX;
DZ=50000/IMAX;

%for I=1:IMAX
% TOPO_PROFILE(I)=0.01*((IMAX/2)-I)^2;
%end
%discretize

%t2=TOPO_PROFILE(:,2);
%t3=interp(t2,2);

%for I=1:IMAX
% for J=1:JMAX
%  C(J,I)=round(M3(I)/DY);
% end 
%end

for t=1:timesteps
%PAD=200;
figure;
A=reshape(ZrsatDif_2D((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
%B=squeeze(A(:,:,KMAX/2));
%C=B';




contourf(A,25)
colorbar
time=num2str(t);
title(['time= ',time'])
M(t)=getframe;
max(max(A))
%B(1:JMAX+PAD,1:IMAX)=0.0;

%for I=1:IMAX
% for J=1:JMAX
%  B(J+C(J,I)+PAD/2,I)=A(J,I);
% end
%end

%XX=X;
%YY=DY:DY:DY*(JMAX+PAD);
%YYY=DY:DY:DY*(JMAX+PAD-150);
%CMODE=mode(mode(C));
%YYYY=YYY(:)-DY.*(JMAX+PAD-150+CMODE);
%Y5=YYYY./1000;
%Y6=Y5+6;
%X5=XX./1000;
%contourf(X5,Y6,B(151:JMAX+PAD,:),20)
%contourf(X,Y,A)
end

%for I=1:IMAX
% for J=1:JMAX
%   if (C(J,I)<600)
%      D(J,I)=0.0;
%      T(J,I)=600.;
%    else
%     D(J,I)=1.0;
%     T(J,I)=C(J,I);
%    end
% end
%end

%pcolor(D)

%D1=D(160:200,60:140);
%T1=T(160:200,60:140);

%pcolor(D1)


%D2=D1;
%T2=T1;

%for I=58:80
% for J=30:40
%   D2(J,I)=0.0;
%   T2(J,I)=600;
% end 
%end

%D3=interp2(D2,2);
%T3=interp2(T2,2);

%size(D3)
%pcolor(D3)

%dlmwrite('bc_filter',D3,'delimiter',' ')
%dlmwrite('temp_filter',T3,'delimiter',' ')

