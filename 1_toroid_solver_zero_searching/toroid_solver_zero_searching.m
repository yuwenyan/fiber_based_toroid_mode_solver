% new toroid mode solver based on normal fiber modes
% zero searching
% 2013.06

clear all;
close all;

tic;

% constant
c_m_per_sec=299792458;
u0=4*pi*1e-7;
e0=8.8542e-12;

% fiber parameter
fiber_geom = toroid_structure;
val_m=fiber_geom.M_guess;
scale_fac=1e13;


toroid_comp_str={'Erad','Eaxi','Eazi','Hrad','Haxi','Hazi'};
fiber_comp_str={'Ex','Ey','Ez','Hx','Hy','Hz'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fiber_res=fiber_geom;
disp('Start computing fiber modes...');
[fiber_modes,fem_fib]=fiber_fem_solver(fiber_res);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fiber_res.mode_count=length(fiber_modes.Neff);
disp(['Total ' int2str(fiber_res.mode_count) ' modes found']);

lambda_range_175=[1560e-9; 1460e-9];
lambda_range_176=[1550e-9; 1450e-9];
lambda_range_177=[1540e-9; 1440e-9];

m_result=val_m;

if m_result==175
    lambda_range=lambda_range_175;
elseif m_result==176
    lambda_range=lambda_range_176;
elseif m_result==177
    lambda_range=lambda_range_177;
end

scan_step=0.5e-9;
lambda_try=(lambda_range(1):-scan_step:lambda_range(2));
lambda_num=length(lambda_try);
lambda_zero=zeros(lambda_num,1);
lambda_det=zeros(lambda_num,1);
compensate=3e122;


disp('Start scanning wavelength');
for lambda_count=1:lambda_num
    % lambda_det(lambda_count)=det(calc_T_det(fiber_res,fiber_modes,fem_fib,lambda_try(lambda_count),m_result));
    lambda_det(lambda_count)=calc_L_det(lambda_try(lambda_count),fiber_res,fiber_modes,fem_fib,m_result,scale_fac);
    if (mod(lambda_count,20)==1)
      %  disp(['\lambda=' num2str(lambda_try(lambda_count)*1e-9) ' [nm]']);
      %  lambda_det(1:lambda_count)
        figure(1);clf;plot(lambda_try(1:lambda_count)*1e-9,lambda_det(1:lambda_count));
        xlabel('\lambda [nm]');
        ylabel('det');
		disp(['Scanning...   ' num2str(lambda_count) '/' num2str(lambda_num)]);
    end
end

disp('Scan completed.');

figure;
hold on;
plot(lambda_try,real(lambda_det));title(['m = ', num2str(m_result)]);xlabel('\lambda');ylabel('linear');
plot(lambda_try,lambda_zero,'r');
set(gca,'XDir','reverse'); % X reverse
hold off;

%%%%%%%%%%%%%%%%%

lambda_finding_min = [];
for lambda_count=2:lambda_num
    if (lambda_det(lambda_count-1)*lambda_det(lambda_count))<=0
        lambda_finding_min(end+1)=lambda_try(lambda_count);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test(1,:)=lambda_finding_min;
test(2,:)=test(1,:)+scan_step;
test_num=length(test(1,:));
test_root=zeros(1,test_num);
for test_count=1:test_num
    test_root(test_count) = fzero(@(x) calc_L_det(x,fiber_res,fiber_modes,fem_fib,m_result,scale_fac),[test(1,test_count) test(2,test_count)]);
	disp(['Calculating...   ' num2str(test_count) '/' num2str(test_num)]);
end % test_count

res_wavelength_nm=test_root'*1e9; % resonance wavelength in nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time=toc;

% %%%%%%%%%%%%%% plot toroid fields %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% fig_id=100;
% res_lambda=test_root;
% res_m=m_result;
% root_num=length(res_lambda);
% Resonance=cell(1,root_num);
% Field_temp=zeros(6,size(fem_fib.mesh.p,2));
% for root_count=1:root_num
%     root_try=res_lambda(root_count);
%     res_root = fminsearch(@(x) real(det(calc_L_mat(fiber_res,fiber_modes,fem_fib,x,res_m,scale_fac)))^2,root_try,optimset('TolX',1e-35));
%     res_L_mat=calc_L_mat(fiber_res,fiber_modes,fem_fib,res_root,res_m,scale_fac);
%     res_V=null(res_L_mat);
%     Resonance{root_count}.wavelength=res_root;
%     Resonance{root_count}.m=res_m;
%     Resonance{root_count}.V=res_V;
%     Resonance{root_count}.Field=zeros(size(Field_temp));
%     
%     for count=1:size(res_V,1);
%         for count_field=1:length(fiber_comp_str)
%             Field_temp(count_field,:)=postinterp(fem_fib,fiber_comp_str{count_field},fem_fib.mesh.p,'solnum',fiber_modes.Neff_indx(count));
%         end %for count_field
%         Resonance{root_count}.Field=Resonance{root_count}.Field+res_V(count)*Field_temp;
%     end %for count
%     fiber_fields_bend=Resonance{root_count}.Field;
%     file_label='Method_one_toroid_';
%     fig_id=fig_id+1;
%     figid=plotfields_mode(fem_fib,fiber_fields_bend,fig_id,file_label);
% end % for root_count
    