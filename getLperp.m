function [ l_perp ] = getLperp( xh,xhp,F )

lp=F*xh;

ap=lp(1);bp=lp(2);cp=lp(3);

xp=xhp(1);yp=xhp(2);wp=xhp(3);

l_perp=[-bp*wp,ap*wp,bp*xp-ap*yp]';

end

