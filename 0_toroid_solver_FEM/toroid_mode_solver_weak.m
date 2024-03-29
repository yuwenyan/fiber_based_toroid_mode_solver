function y=toroid_mode_solver_weak(toroid_struct)
%toroid_struct=[];
% COMSOL Multiphysics Model M-file
% Generated by COMSOL 3.4 (COMSOL 3.4.0.250, $Date: 2010/01/28 21:23:08 $)
% Some geometry objects are stored in a separate file.
% The name of this file is given by the variable 'flbinaryfile'.

flclear fem

% COMSOL version
clear vrsn
vrsn.name = 'COMSOL 3.4';
vrsn.ext = '';
vrsn.major = 0;
vrsn.build = 250;
vrsn.rcs = '$Name:  $';
vrsn.date = '$Date: 2010/01/28 21:23:08 $';
fem.version = vrsn;


major_R=toroid_struct.toroid_R;
minor_R=toroid_struct.core_radius;
wnd_sz_X=10.*minor_R;
wnd_sz_Y=10.*minor_R;

n_core=toroid_struct.n_guide;
n_cladding=toroid_struct.n_sub;
lambda_m=toroid_struct.lambda;
wnd_corner_X=-wnd_sz_X/2+major_R;
wnd_corner_Y=-wnd_sz_Y/2;

if isfield(toroid_struct,'wafer_thickness')
    wafer_thickness=toroid_struct.wafer_thickness;
else
    wafer_thickness=0.5e-6;
end
if isfield(toroid_struct,'num_modes')
    num_modes=toroid_struct.num_modes;
else
    num_modes=10;
end
if isfield(toroid_struct,'M_guess')
    M_Guess=toroid_struct.M_guess;
else
    M_Guess=round(2*pi*major_R*n_core/lambda_m);
end
%end

c_m_per_sec=299792458;  % light speed
f_Hz=c_m_per_sec/lambda_m;

% Constants
fem.const = {'c',c_m_per_sec, ...
    'k','2*pi/c', ... % cbar = k
    'fc','k^2', ...
    'alpha','1.0', ...
    'M',M_Guess, ...
    'omega',2*pi*f_Hz, ...
    'delta_e','0.0', ...
    'e1',[num2str(n_core) '^2*(1+delta_e)'], ...
    'e2',n_cladding^2, ...
    'delta_eperp1','0*1e-3', ...
    'eperp1','9.2725*(1+delta_eperp1)', ...
    'delta_epara1','0*1e-3', ...
    'epara1','11.3486*(1+delta_epara1)', ...
    'eperp2','1.0', ...   % water=76.75 @ 303K, air=1.0
    'epara2','1.0', ...   % water
    'e0','8.8542e-12',...
    'u0','pi*4e-7',...
    'e_293K_alumina','9.8', ...
    'eperp_4K_sapph_UWA','9.2725', ...
    'epara_4K_sapph_UWA','11.3486', ...
    'eperp_293K_sapph','9.407', ...
    'epara_293K_sapph','11.62', ...
    'eperp_4K_sapph_NPL','9.2848', ...
    'epara_4K_sapph_NPL','11.3660', ...
    'n_silica','1.4457', ...
    'n_AlGaAs','3.36', ...
    'n_core',n_core,...
    'n_cladding',n_cladding,...
    'mf','2.374616e14', ...
    'ttgH','1', ...
    'ttgE','0', ...
    'rectangle_mf','2.376629e14', ...
    'circle_mf','2.374616e14', ...
    'mixing_angle','45', ...
    'cMW','sin(mixing_angle * pi /180)', ...
    'cEW','cos(mixing_angle * pi /180)', ...
    'tngM','1', ...
    'tngE','0'};

clear draw

% Geometry
g_membrane=rect2(major_R-wnd_corner_X,wafer_thickness,'base','corner','pos',[wnd_corner_X,-wafer_thickness/2]);
g_ring=ellip2(minor_R,minor_R,'base','center','pos',{major_R,0},'rot','0');
g_toroid=geomcomp({g_membrane,g_ring},'ns',{'R2','E1'},'sf','R2+E1','edge','all');
g_background=rect2(wnd_sz_X,wnd_sz_Y,'base','corner','pos',[wnd_corner_X,wnd_corner_Y]);

% Analyzed geometry
clear s
s.objs={g_toroid,g_background}; %
s.name={'CO1','R1'};
s.tags={'g_toroid','g_background'}; %

fem.draw=struct('s',s);
fem.geom=geomcsg(fem);

% Initialize mesh
fem.mesh=meshinit(fem, ...
    'hauto',5);

%     % Refine mesh
%     fem.mesh=meshrefine(fem, ...
%                         'mcase',0, ...
%                         'rmethod','regular');
% Refine mesh
fem.mesh=meshrefine(fem, ...
    'mcase',0, ...
    'rmethod','regular');
%     % Refine mesh
%     fem.mesh=meshrefine(fem, ...
%                     'mcase',0, ...
%                     'boxcoord',[major_R-1.5*minor_R major_R+1.5*minor_R -1.5*minor_R 1.5*minor_R], ...
%                     'rmethod','regular');


% (Default values are not included)

% Application mode 1
clear appl
appl.mode.class = 'FlPDEW';
appl.dim = {'Hrad','Hazi','Haxi','Hrad_t','Hazi_t','Haxi_t'};
appl.sdim = {'r','z','z2'};
appl.name = 'Axisymmetric';
appl.gporder = 4;
appl.cporder = 2;
appl.assignsuffix = '_Axisymmetric';
clear bnd
bnd.constrf = {{'test(Hrad*nr+Haxi*nz)';'test(-Haxir+Hradz)';...
    'test(-(Hazi*nr-Hrad*M*nr-Haxi*M*nz+Hazir*nr*r+Haziz*nz*r)/r)'},...
    {0;'test(-Haxir+Hradz)';'test(-(Hazi*nr-Hrad*M*nr-Haxi*M*nz+Hazir*nr*r+Haziz*nz*r)/r)'},...
    'test(Hrad*nr+Haxi*nz)',{'test(Haxi*nr-Hrad*nz)';'test(Hazi)';'test(-(Haxi*M*nr+Hazi*nz-Hrad*M*nz-Haziz*nr*r+Hazir*nz*r)/r)'},...
    {'test(Haxi*nr-Hrad*nz)';'test(Hazi)'},{0;0;'test(-(Haxi*M*nr+Hazi*nz-Hrad*M*nz-Haziz*nr*r+Hazir*nz*r)/r)'},...
    {'test(-i*cMW*Hazi*k*mf+cEW*(Hazi*nr-Hrad*M*nr-Haxi*M*nz+Hazir*nr*r+Haziz*nz*r)/r)';...
    'test(-i*cEW*(-Haxir+Hradz)+cMW*k*mf*(Haxi*nr-Hrad*nz))'},0,{'test(-Hrad)'; ...
    'test(-Hazi)';'test(-Haxi)'}};
bnd.name ={'electric_wall','normal_D','tangential_H','magnetic_wall','normal_H', ...
    'tangential_D','radiation_match','null',''};
bnd.constr ={{'Hrad*nr+Haxi*nz';'-Haxir+Hradz';'-(Hazi*nr-Hrad*M*nr-Haxi*M*nz+Hazir*nr*r+Haziz*nz*r)/r'},...
    {0;'-Haxir+Hradz';'-(Hazi*nr-Hrad*M*nr-Haxi*M*nz+Hazir*nr*r+Haziz*nz*r)/r'}, ...
    'Hrad*nr+Haxi*nz',{'Haxi*nr-Hrad*nz';'Hazi';'-(Haxi*M*nr+Hazi*nz-Hrad*M*nz-Haziz*nr*r+Hazir*nz*r)/r'},...
    {'Haxi*nr-Hrad*nz';'Hazi'},{0;0;'-(Haxi*M*nr+Hazi*nz-Hrad*M*nz-Haziz*nr*r+Hazir*nz*r)/r'},...
    {'-i*cMW*Hazi*k*mf+(cEW*(Hazi*nr-Hrad*M*nr-Haxi*M*nz+Hazir*nr*r+Haziz*nz*r))/r';...
    '-i*cEW*(-Haxir+Hradz)+cMW*k*mf*(Haxi*nr-Hrad*nz)'},0,{'-Hrad';'-Hazi'; ...
    '-Haxi'}};
%     bnd.ind = [7,7,7,9,7,9,7,7,9,9,9,9];
bnd.ind = [3,3,3,8,3,8,3,3,8,8,8,8];
appl.bnd = bnd;
clear equ
equ.dweak = 'fc*r*(Haxitt*test(Haxi)+Hazitt*test(Hazi)+Hradtt*test(Hrad))';
equ.name ={'isotrop_diel_2','isotrop_diel_1','dielectric_0:vacuum','uniaxial_diel_1', ...
    'uniaxial_diel_2'};
equ.weak ={{'(-(Haziz*M*test(Haxi))+Hazir*test(Hazi)+Hazi*test(Hazir)-Hrad*M*test(Hazir)-Haxi*M*test(Haziz)-Hazir*M*test(Hrad)+(Haxi*M^2*test(Haxi)+(Hazi-Hrad*M)*(test(Hazi)-M*test(Hrad)))/r+r*((Haxir-Hradz)*test(Haxir)+Hazir*test(Hazir)+Haziz*test(Haziz)-Haxir*test(Hradz)+Hradz*test(Hradz)))/e2';...
    'alpha*(Hrad*test(Haxiz)-Hazi*M*test(Haxiz)-Haxiz*M*test(Hazi)-Hradr*M*test(Hazi)+Haxiz*test(Hrad)+Hradr*test(Hrad)+(-(Hrad*M*test(Hazi))+Hazi*M^2*test(Hazi)+Hrad*test(Hrad)-Hazi*M*test(Hrad))/r+Hrad*test(Hradr)-Hazi*M*test(Hradr)+r*(Haxiz*test(Haxiz)+Hradr*test(Haxiz)+Haxiz*test(Hradr)+Hradr*test(Hradr)))'},...
    {'(-(Haziz*M*test(Haxi))+Hazir*test(Hazi)+Hazi*test(Hazir)-Hrad*M*test(Hazir)-Haxi*M*test(Haziz)-Hazir*M*test(Hrad)+(Haxi*M^2*test(Haxi)+(Hazi-Hrad*M)*(test(Hazi)-M*test(Hrad)))/r+r*((Haxir-Hradz)*test(Haxir)+Hazir*test(Hazir)+Haziz*test(Haziz)-Haxir*test(Hradz)+Hradz*test(Hradz)))/e1';...
    'alpha*(Hrad*test(Haxiz)-Hazi*M*test(Haxiz)-Haxiz*M*test(Hazi)-Hradr*M*test(Hazi)+Haxiz*test(Hrad)+Hradr*test(Hrad)+(-(Hrad*M*test(Hazi))+Hazi*M^2*test(Hazi)+Hrad*test(Hrad)-Hazi*M*test(Hrad))/r+Hrad*test(Hradr)-Hazi*M*test(Hradr)+r*(Haxiz*test(Haxiz)+Hradr*test(Haxiz)+Haxiz*test(Hradr)+Hradr*test(Hradr)))'},...
    {'-(Haziz*M*test(Haxi))+Hazir*test(Hazi)+Hazi*test(Hazir)-Hrad*M*test(Hazir)-Haxi*M*test(Haziz)-Hazir*M*test(Hrad)+(Haxi*M^2*test(Haxi)+(Hazi-Hrad*M)*(test(Hazi)-M*test(Hrad)))/r+r*((Haxir-Hradz)*test(Haxir)+Hazir*test(Hazir)+Haziz*test(Haziz)-Haxir*test(Hradz)+Hradz*test(Hradz))';...
    'alpha*(Hrad*test(Haxiz)-Hazi*M*test(Haxiz)-Haxiz*M*test(Hazi)-Hradr*M*test(Hazi)+Haxiz*test(Hrad)+Hradr*test(Hrad)+(-(Hrad*M*test(Hazi))+Hazi*M^2*test(Hazi)+Hrad*test(Hrad)-Hazi*M*test(Hrad))/r+Hrad*test(Hradr)-Hazi*M*test(Hradr)+r*(Haxiz*test(Haxiz)+Hradr*test(Haxiz)+Haxiz*test(Hradr)+Hradr*test(Hradr)))'},...
    {'(-(epara1*Haziz*M*test(Haxi))+eperp1*Hazir*test(Hazi)+eperp1*Hazi*test(Hazir)-eperp1*Hrad*M*test(Hazir)-epara1*Haxi*M*test(Haziz)-eperp1*Hazir*M*test(Hrad))/(epara1*eperp1)+(epara1*Haxi*M^2*test(Haxi)+eperp1*Hazi*test(Hazi)-eperp1*Hrad*M*test(Hazi)-eperp1*Hazi*M*test(Hrad)+eperp1*Hrad*M^2*test(Hrad))/(epara1*eperp1*r)+(r*(epara1*(Haxir-Hradz)*test(Haxir)+eperp1*Hazir*test(Hazir)+epara1*Haziz*test(Haziz)-epara1*Haxir*test(Hradz)+epara1*Hradz*test(Hradz)))/(epara1*eperp1)';...
    'alpha*(Hrad*test(Haxiz)-Hazi*M*test(Haxiz)-Haxiz*M*test(Hazi)-Hradr*M*test(Hazi)+Haxiz*test(Hrad)+Hradr*test(Hrad)+(-(Hrad*M*test(Hazi))+Hazi*M^2*test(Hazi)+Hrad*test(Hrad)-Hazi*M*test(Hrad))/r+Hrad*test(Hradr)-Hazi*M*test(Hradr)+r*(Haxiz*test(Haxiz)+Hradr*test(Haxiz)+Haxiz*test(Hradr)+Hradr*test(Hradr)))'},...
    {'(-(epara2*Haziz*M*test(Haxi))+eperp2*Hazir*test(Hazi)+eperp2*Hazi*test(Hazir)-eperp2*Hrad*M*test(Hazir)-epara2*Haxi*M*test(Haziz)-eperp2*Hazir*M*test(Hrad))/(epara2*eperp2)+(epara2*Haxi*M^2*test(Haxi)+eperp2*Hazi*test(Hazi)-eperp2*Hrad*M*test(Hazi)-eperp2*Hazi*M*test(Hrad)+eperp2*Hrad*M^2*test(Hrad))/(epara2*eperp2*r)+(r*(epara2*(Haxir-Hradz)*test(Haxir)+eperp2*Hazir*test(Hazir)+epara2*Haziz*test(Haziz)-epara2*Haxir*test(Hradz)+epara2*Hradz*test(Hradz)))/(epara2*eperp2)';...
    'alpha*(Hrad*test(Haxiz)-Hazi*M*test(Haxiz)-Haxiz*M*test(Hazi)-Hradr*M*test(Hazi)+Haxiz*test(Hrad)+Hradr*test(Hrad)+(-(Hrad*M*test(Hazi))+Hazi*M^2*test(Hazi)+Hrad*test(Hrad)-Hazi*M*test(Hrad))/r+Hrad*test(Hradr)-Hazi*M*test(Hradr)+r*(Haxiz*test(Haxiz)+Hradr*test(Haxiz)+Haxiz*test(Hradr)+Hradr*test(Hradr)))'}};

equ.ind = [1,2];
appl.equ = equ;
fem.appl{1} = appl;
fem.sdim = {'r','z'};
fem.border = 1;
clear units;
units.basesystem = 'SI';
fem.units = units;

% Subdomain settings
clear equ
equ.ind = [1,2];
equ.dim = {'Hrad','Hazi','Haxi'};

% Subdomain expressions
equ.expr = {'erel',{'e2','e1'},...
    'core',{'delta_e','alpha'},...
    'incore',{0,1},...
    'outtoroid',{1,0},...
    'nindex',{'n_cladding','n_core'},...
    'beta_x',{'n_cladding*x','n_core*x'}};
%        'EnDensBackground',{'EnDens',0}};
fem.equ = equ;

% Scalar expressions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fem.expr = {'DivH','(Hrad-Hazi*M+(Haxiz+Hradr)*r)/r', ...
    'Drad','(Haxi*M-Haziz*r)/r/omega', ...
    'Dazi','(-Haxir+Hradz)/omega', ...
    'Daxi','(Hazi-Hrad*M+Hazir*r)/r/omega', ...
    'Erad','Drad/erel/e0', ...
    'Eazi','Dazi/erel/e0', ...
    'Eaxi','Daxi/erel/e0', ...
    'Brad','u0*Hrad', ...
    'Bazi','u0*Hazi', ...
    'Baxi','u0*Haxi', ...
    'ElecMag','sqrt(Erad*conj(Erad)+Eazi*conj(Eazi)+Eaxi*conj(Eaxi))', ...
    'ElecMag2','sqrt(Erad*conj(Erad)+Eaxi*conj(Eaxi))', ...
    'erp1','erel', ...
    'comment','1', ...
    'MagAziSqrd','imag(Hazi)^2', ...
    'MagTransSqrd','real(Haxi)^2+real(Hrad)^2', ...
    'ElecAziSqrd','real(Eazi)^2', ...
    'I_tot','(erel*r)*(abs(Eaxi)^2+abs(Erad)^2)',...
    'ElecTransSqrd','abs(Eaxi)^2+abs(Erad)^2',...
    'Poynting_rad','0.5*(real(Eazi*conj(Haxi)-Eaxi*conj(Hazi)))',... % Srad S=E x H
    'Poynting_axi','0.5*(real(Erad*conj(Hazi)-Eazi*conj(Hrad)))',... % Saxi
    'Poynting_azi','0.5*(real(Eaxi*conj(Hrad)-Erad*conj(Haxi)))',... %
    'intensity','0.5*abs(real(Eaxi*conj(Hrad)-Erad*conj(Haxi)))',... %
    'Poynting_r','Poynting_rad*(r/sqrt(r^2+z^2))+Poynting_axi*(z/sqrt(r^2+z^2))',...
    'Poynting','(Poynting_rad^2+Poynting_axi^2)^(1/2)',...
    'MagEnDens','Hrad*conj(Brad)+Hazi*conj(Bazi)+Haxi*conj(Baxi)',... % magnetic energy density
    'ElecField','(Erad*conj(Erad)+Eazi*conj(Eazi)+Eaxi*conj(Eaxi))^(1/2)',...% electric field
    'ElecEnDens','Erad*conj(Drad)+Eazi*conj(Dazi)+Eaxi*conj(Daxi)',...  % electric energy density
    'EnDens','MagEnDens+ElecEnDens'};   % energy density

% Descriptions
clear descr
descr.expr= {'MagTransSqrd','transverse component of magnetic fieldsquared','ElecTransSqrd',...
    'transverse component of electric fieldsquared','ElecAziSqrd','azimuthal component of electric fieldsquared',...
    'MagAziSqrd','azimuthal component of magnetic fieldsquared','Daxi','axial component of electric displacement',...
    'DivH','divergence ofmagnetic field (should be zero!)','Erad','radial component of electric fieldstrength','Drad',...
    'radial component of electricdisplacement','comment','elemental volume = 2 pi r d_r d_phi','Eaxi',...
    'axialcomponent of electric field strength','Dazi','azimuthal component of electricdisplacement',...
    'Eazi','azimuthal component of electric field strength'};
fem.descr = descr;

% Descriptions
descr = fem.descr;
descr.const= {'n_silica','refractive index of thermally grown silica (Fig B.2,p. 172 of Kippenberg''s thesis)',...
    'eperp1','relative permittivity ofuniaxial_dielectric_1 perpendicular to cylindricalaxis',...
    'e_293K_alumina','relative permittivity of alumina at roomtemperature','c','speed of light (exact!)',...
    'delta_epara1','fractional increment(for determining filling factors)','delta_eperp1',...
    'fractional increment (for determining filling factors)','eperp2',...
    'relative permittivity of uniaxial_dielectric_2 perpendicular to cylindrical axis','M','azimuthal mode order',...
    'fc','constant used internally --do not modify','eperp_4K_sapph_UWA','UWA values for cryogenic HEMEX sapphire',...
    'delta_e','fractional increment (for determining filling factors)','e1','relative permittivity of isotropic_dielectric_1',...
    'epara1','relative permittivity of uniaxial_dielectric_1 parallel to cylindrical axis',...
    'e2','ditto for isotropic_dielectric_2','epara2','ditto but parallel to cylindrical axis',...
    'cMW','Magnetic-Wall-ness','alpha','penalty coefficient on Div H','eperp_293K_sapph',...
    'nominal room temperature values for same','eperp_4K_sapph_NPL','NPL values','mf',...
    'match frequency','n_AlGaAs','average refractive index of GaAs and AlGaAs layers (p. 172 of Srinivasan)',...
    'cEW','Electric-Wall-ness','mixing_angle','Electric-Magnetic Mixing Angle (in degrees)'};
fem.descr = descr;

% Multiphysics
fem=multiphysics(fem);

% Extend mesh
fem.xmesh=meshextend(fem, ...
    'linshape',[]);

% Solve problem
fem.sol=femeig(fem, ...
    'conjugate','on', ...
    'symmetric','on', ...
    'solcomp',{'Hazi','Hrad','Haxi'}, ...
    'outcomp',{'Hazi','Hrad','Haxi'}, ...
    'neigs',num_modes, ...
    'shift',f_Hz, ...
    'linsolver','spooles');

% Save current fem structure for restart purposes
fem0=fem;
[sol,sol_ind]=sort(fem.sol.lambda,'ascend');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for count=1:num_modes
    figure;clf;
    % Plot solution
    postplot(fem,'tridata',{'ElecTransSqrd','cont','internal'}, ...
        'trimap','jet(1024)', ... %'contdata',{'log10(Poynting+1e-2)','cont','internal'}, ...
        'contlevels',10, ...
        'contlabel','off', ...
        'contmap','cool(1024)', ...
        'arrowdata',{'Hrad','Haxi'}, ...
        'arrowxspacing',15, ...
        'arrowyspacing',13, ...
        'arrowscale',1.2, ...
        'arrowtype','arrow', ...
        'arrowstyle','proportional', ...
        'arrowcolor',[1.0,1.0,1.0], ...
        'maxminsub','Poynting', ...
        'solnum',sol_ind(count), ...
        'phase',(90)*pi/180, ...
        'title',['Order :' int2str(M_Guess) '; lambda=' ...
        num2str(c_m_per_sec/fem.sol.lambda(sol_ind(count))*1e9,'%8.5f') '[nm]   Surface: ElecTransSqrd'], ...
        'axis',[wnd_corner_X,wnd_sz_X+wnd_corner_X,wnd_corner_Y,wnd_sz_Y+wnd_corner_Y]);
end

y=fem0;
end