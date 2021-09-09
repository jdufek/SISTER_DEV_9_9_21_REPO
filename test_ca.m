CaO = bin2mat('MELT_CAO','double','b');
SiO2= bin2mat('MELT_SIO2','double','b');
MF=bin2mat('MELT_FRACTION','double','b'); 


L=length(MF)

for i=1:L
 if (MF(i)>0.0)
  if (CaO<.0000001)
    i
  end
 end
end




