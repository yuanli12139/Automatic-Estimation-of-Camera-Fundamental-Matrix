function [ php,Pp_corrected,Xh,Xpi_corrected ] = getFrPV( PV )
% Get updated Pp and Xpi from the Parameter Vector 
php=PV(1:11);
Pp_corrected=reshape(deParamtrz(php),4,3)';

n=length(PV(12:end))/3;
Xh=reshape(PV(12:end),3,n);

Xpi_corrected=zeros(4,n);
for pt=1:n
    
    Xpi_corrected(:,pt)=deParamtrz(Xh(:,pt));

end

