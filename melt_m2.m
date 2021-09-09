load MAFIC_MELT
load M3


t=5
IMAX=1000;
JMAX=1000;




X=50000/IMAX:50000/IMAX:50000;
Y=60000/JMAX:60000/JMAX:60000;
DX=50000/IMAX;
DY=60000/JMAX;

%for I=1:IMAX
% TOPO_PROFILE(I)=0.01*((IMAX/2)-I)^2;
%end
%discretize

%t2=TOPO_PROFILE(:,2);
%t3=interp(t2,2);

for I=1:IMAX
 for J=1:JMAX
  C(J,I)=round(M3(I)/DY);
 end 
end


PAD=200;

A=reshape(MAFIC_MELT((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[IMAX JMAX])';

B(1:JMAX+PAD,1:IMAX)=0.0;

for I=1:IMAX
 for J=1:JMAX
  B(J+C(J,I)+PAD/2,I)=A(J,I);
 end
end

XX=X;
YY=DY:DY:DY*(JMAX+PAD);
YYY=DY:DY:DY*(JMAX+PAD-150);
CMODE=mode(mode(C));
YYYY=YYY(:)-DY.*(JMAX+PAD-150+CMODE);
Y5=YYYY./1000;
Y6=Y5+6;
X5=XX./1000;
contourf(X5,Y6,B(151:JMAX+PAD,:),20)
%contourf(X,Y,A)
axis equal 
axis tight
colorbar

set( gca                       , ...
    'FontWeight' , 'bold'      , ...
    'FontName'   , 'Helvetica' );

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , -60:10:0, ...
  'XTick'       , 0:10:50, ...
  'LineWidth'   , 1         );


hTitle=title('Temperature');
hXLabel=xlabel( ' km ')
hYLabel=ylabel( 'km ')

filename='temp_ldm';
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 filename
print(filename,'-dtiff')
print(filename,'-depsc2')
