%$Id: area_integral.m,v 1.2 2006/10/07 00:59:56 taolu Exp $
%$Revision: 1.2 $
%$Author: taolu $
%$Date: 2006/10/07 00:59:56 $
function y=area_integral(fem,f)
[trig_area,a2,a3,a4]=pdetrg(fem.mesh.p,fem.mesh.t);
f_mid=pdeintrp(fem.mesh.p,fem.mesh.t,transpose(f));
%f_mid=(f(:,fem.mesh.t(1,:))+f(:,fem.mesh.t(2,:))+f(:,fem.mesh.t(3,:)))/3;
y=sum(f_mid.*(ones(size(f_mid,1),1)*trig_area),2);

