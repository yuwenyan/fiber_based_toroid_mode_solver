%$Id: fiber_fem_analytic.m,v 1.11 2006/10/27 23:55:21 taolu Exp $
%$Revision: 1.11 $
%$Author: taolu $
%$Date: 2006/10/27 23:55:21 $
%function [analytic_results,fem0]=fiber_fem_analytic(parm_str)
function [analytic_results,fem0]=fiber_fem_solver_graded(parm_str)

% [x_grids,y_grids]=meshgrid(parm_str.mesh_X,parm_str.mesh_Y);
% p_grids=[transpose(x_grids(:)); transpose(y_grids(:))];
analytic_results.Field=[];
V_num=parm_str.k0_a*sqrt(parm_str.n_guide^2-parm_str.n_sub^2);
num_modes=parm_str.num_modes;%round(1.2*V_num^2/2);

flclear fem

% COMSOL version
clear vrsn_fiber
vrsn_fiber.name = 'COMSOL 3.3';
vrsn_fiber.ext = '';
vrsn_fiber.major = 0;
vrsn_fiber.build = 405;
vrsn_fiber.rcs = '$Name:  $';
vrsn_fiber.date = '$Date: 2006/10/27 23:55:21 $';
fem.version = vrsn_fiber;

% Geometry
calc_radius_fiber=5*parm_str.core_radius;
g1_fiber=circ2(parm_str.core_radius,'base','center','pos',[0,0]);
g2_fiber=circ2(calc_radius_fiber,'base','center','pos',[0,0]);

clear s_fiber
s_fiber.objs={g1_fiber,g2_fiber};
s_fiber.name={'C1_fiber','C2_fiber'};
s_fiber.tags={'g1_fiber','g2_fiber'};

fem.draw=struct('s',s_fiber);
fem.geom=geomcsg(fem);
% geomplot(fem,'pointmode','off');

% Constants
fem.const = {'bend_R',parm_str.toroid_R,...
    'n_guide',parm_str.n_guide,...
    'n_sub',parm_str.n_sub,...
    'e0',8.8542e-12,...
    'u0',4*pi*1e-7,...
    'lambda',parm_str.lambda};

% Initialize mesh
fem.mesh=meshinit(fem, ...
                  'hauto',5);

% Refine mesh
% fem.mesh=meshrefine(fem, ...
%                     'mcase',0, ...
%                     'rmethod','regular');

% fem.mesh=meshrefine(fem, ...
%                     'mcase',0, ...
%                     'rmethod','regular');

[trig_area,a2,a3,a4]=pdetrg(fem.mesh.p,fem.mesh.t);

% Application mode 1
clear appl_fiber
appl_fiber.mode.class = 'PerpendicularWaves';
appl_fiber.module = 'RF';
appl_fiber.gporder = 4;
appl_fiber.cporder = 2;
appl_fiber.assignsuffix = '_rfwv';
clear prop_fiber
prop_fiber.elemdefault='Lag2';
prop_fiber.inputvar='lambda';
prop_fiber.comps='2';
appl_fiber.prop = prop_fiber;
clear bnd_fiber
bnd_fiber.type = {'H0','cont'};
bnd_fiber.ind = [1,1,2,2,1,2,2,1];
appl_fiber.bnd = bnd_fiber;
clear equ_fiber
equ_fiber.n = {'sqrt(abs(n_sub^2+n_sub^2-n_sub^2*bend_R^2/(bend_R+x)^2))','sqrt(abs(n_guide^2+n_guide^2-n_guide^2*bend_R^2/(bend_R+x)^2))'};
equ_fiber.ind = [1,2];
appl_fiber.equ = equ_fiber;
appl_fiber.var = {'lambda0',num2str(parm_str.lambda,8)};
fem.appl{1} = appl_fiber;
fem.frame = {'ref'};
fem.border = 1;
clear units_fiber;
units_fiber.basesystem = 'SI';
fem.units = units_fiber;

% Subdomain settings
clear equ
equ.ind = [1,2];
equ.dim = {'Hx','Hy','hz'};
equ.expr = {'nindex',{'n_sub', 'n_guide'},...
    'incore',{0,1},...
    'outcore',{1,0}};
fem.equ = equ;

% Scalar expressions
    fem.expr = {'Poynting_x','0.5*(real(Ey*conj(Hz)-Ez*conj(Hy)))',... % Srad S=E x H
      'Poynting_y','0.5*(real(Ez*conj(Hx)-Ex*conj(Hz)))',...
      'Poynting_z','0.5*(real(Ex*conj(Hy)-Ey*conj(Hx)))',...
      'intensity','0.5*abs(real(Ex*conj(Hy)-Ey*conj(Hx)))',...
      'ElecMag','sqrt(Ex*conj(Ex)+Ey*conj(Ey)+Ez*conj(Ez))',...
      'ElecMag2','sqrt(Ex*conj(Ex)+Ey*conj(Ey))',...
      'ElecEnDens','0.5*(Ex*conj(Ex)+Ey*conj(Ey)+Ez*conj(Ez))*nindex*nindex*e0',...
      'ElecTransSqrd','abs(Ex)^2+abs(Ey)^2',...
      'MagEnDens','0.5*(Hx*conj(Hx)+Hy*conj(Hy)+Hz*conj(Hz))*u0',...
      'EnDens','ElecEnDens+MagEnDens',...
      'Poynting_r','Poynting_x*(x/sqrt(x^2+y^2))+Poynting_y*(y/sqrt(x^2+y^2))',...
      'Poynting_rho','Poynting_x*(x/sqrt(x^2+y^2))+Poynting_y*(y/sqrt(x^2+y^2))'};

% Solution form
fem.solform = 'weak';

% Multiphysics
fem=multiphysics(fem);

% Extend mesh
fem.xmesh=meshextend(fem);

shift_val=-1j*(2*pi*parm_str.n_guide/parm_str.lambda); %%%%%?????

% Solve problem
fem.sol=femeig(fem, ...
               'solcomp',{'Hy','Hx'}, ...
               'outcomp',{'Hy','Hx'}, ...
               'neigs',num_modes, ...
               'shift',shift_val);

% Save current fem structure for restart purposes

fem0=fem;
analytic_results.Neff2=postinterp(fem,'neff_rfwv',[0;0],'solnum','all');
[neff_Val,neff_indx]=sort(real(analytic_results.Neff2),'descend');
Neff_indx=[];
analytic_results.num_modes=length(neff_indx);
for count=1:length(neff_indx)
    Pin=postint(fem,'Pozav_rfwv*incore','solnum',neff_indx(count));
    Pout=postint(fem,'Pozav_rfwv*outcore','solnum',neff_indx(count));
    ratio=Pin/Pout;
    if (ratio>1)
        Neff_indx(end+1)=neff_indx(count);
    end
end
analytic_results.Neff=analytic_results.Neff2(Neff_indx);

field_str={'Ex','Ey','Ez','Hx','Hy','Hz'};
Field2=zeros(6,size(fem.mesh.p,2));
power=zeros(length(Neff_indx),1);
Prad=zeros(length(Neff_indx),1);
int_factor=zeros(length(Neff_indx),1);
for count=1:length(Neff_indx)
    % Plot solution
    for count_field=1:length(field_str)
        Field2(count_field,:)=postinterp(fem,field_str{count_field},fem.mesh.p,'solnum',Neff_indx(count));
    end
    %     analytic_results.Field{count}=Field2/sqrt(sum(area_integral(fem0,conj(Field2([4 5],:)).*Field2([4 5],:))));
    %     int_factor(count)=postint(fem,'Pozav_rfwv','solnum',Neff_indx(count));
    int_factor(count)=sum(area_integral(fem0,conj(Field2([4 5],:)).*Field2([4 5],:)));
    analytic_results.Field{count}=Field2/sqrt(int_factor(count));
%     [power(count) Prad(count)]=calc_loss(parm_str,fem0,Neff_indx(count),analytic_results.Neff(count),int_factor(count));
%  if count==1
%     figure(100);clf;
%     postplot(fem, ...
%         'tridata',{'intensity','cont','internal','unit','W/m^2'}, ...
%         'trimap','Rainbow', ...
%         'solnum',Neff_indx(count), ...
%         'title',['mode_num:' int2str(count) '  neff_rfwv   Surface: [W/m^2]'], ...
%         'axis',[-1.5E-5,1.5E-5,-1.5E-5,1.5E-5]);
%  end
%     pause
   
end
analytic_results.power=power;
analytic_results.Prad=Prad;
analytic_results.int_factor=int_factor;
% analytic_results.int_factor2=int_factor2;
analytic_results.Neff_indx=Neff_indx;
analytic_results.V_num=V_num;
% disp(['fiber neff: ' num2str(transpose(analytic_results.Neff))]);


% % Plot solution
% for count=1:length(Neff_indx)
%     figure(100+count);clf;
%     postplot(fem, ...
%         'tridata',{'Pozav_rfwv','cont','internal','unit','W/m^2'}, ...
%         'trimap','Rainbow', ...
%         'solnum',Neff_indx(count), ...
%         'title',['mode_num:' int2str(count)  '  neff=' int2str(neff_Val(count))  '   Surface: Ex]'], ...
%         'axis',[-calc_radius_fiber,calc_radius_fiber,-calc_radius_fiber,calc_radius_fiber]);
%     
% %     figure(200+count);clf;
% %     postplot(fem, ...
% %         'tridata',{'Ey','cont','internal','unit','W/m^2'}, ...
% %         'trimap','Rainbow', ...
% %         'solnum',Neff_indx(count), ...
% %         'title',['mode_num:' int2str(count)  '  neff=' int2str(real(analytic_results.Neff2(count)))  '   Surface: Ey'], ...
% %         'axis',[-calc_radius_fiber,calc_radius_fiber,-calc_radius_fiber,calc_radius_fiber]);
% end

% % Plot solution
% figure(100);clf;
% postplot(fem, ...
%          'tridata',{'Pozav_rfwv','cont','internal','unit','W/m^2'}, ...
%          'trimap','Rainbow', ...
%          'solnum',Neff_indx(1), ...
%          'title',['mode_num:' int2str(1) 'neff_rfwv(1)=1.442448    Surface: Power flow, time average, z component [W/m^2]'], ...
%          'axis',[-calc_radius_fiber,calc_radius_fiber,-calc_radius_fiber,calc_radius_fiber]);
%      
%      
%      for count=1:length(field_str)
%          figure(10+count);clf;
%          postplot(fem, ...
%              'tridata',{field_str{count},'cont','internal','unit','W/m^2'}, ...
%              'trimap','Rainbow', ...
%              'solnum',Neff_indx(1), ...
%              'title',['mode_num:' int2str(1) 'neff_rfwv(1)=1.442448    Surface: ' field_str{count} ], ...
%              'axis',[-calc_radius_fiber,calc_radius_fiber,-calc_radius_fiber,calc_radius_fiber]);
%      end
%      
%      figure(17);clf;
%      postplot(fem, ...
%          'tridata',{'imag(Ez)','cont','internal','unit','W/m^2'}, ...
%          'trimap','Rainbow', ...
%          'solnum',Neff_indx(1), ...
%          'title',['mode_num:' int2str(1) 'neff_rfwv(1)=1.442448    Surface: imag(Ez)' ], ...
%          'axis',[-calc_radius_fiber,calc_radius_fiber,-calc_radius_fiber,calc_radius_fiber]);
%      
%     figure(18);clf;
%      postplot(fem, ...
%          'tridata',{'imag(Hz)','cont','internal','unit','W/m^2'}, ...
%          'trimap','Rainbow', ...
%          'solnum',Neff_indx(1), ...
%          'title',['mode_num:' int2str(1) 'neff_rfwv(1)=1.442448    Surface: imag(Hz)' ], ...
%          'axis',[-calc_radius_fiber,calc_radius_fiber,-calc_radius_fiber,calc_radius_fiber]);


