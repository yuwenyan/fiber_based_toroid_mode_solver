% toroid mode solver based on normal fiber modes
% eigen solutions
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
num_toroid_mode=6;
scale_fac=1;

toroid_comp_str={'Erad','Eaxi','Eazi','Hrad','Haxi','Hazi'};
fiber_comp_str={'Ex','Ey','Ez','Hx','Hy','Hz'};

m_result=val_m; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res_m=m_result;
num_mode_calculation=num_toroid_mode;
res_lambda_record=zeros(num_mode_calculation,1);
iteratioin=cell(1,num_mode_calculation);
Resonance=cell(1,num_mode_calculation);


for mode_order_count=1:num_mode_calculation
    
    disp(['calculating No.' num2str(mode_order_count) ' of ' num2str(num_mode_calculation) ' toroid modes' ]);
    fiber_res=fiber_geom;
    
    iteration_num=30;
    iteration_count=0;
    iteration_wavelength_diff=0.002e-9;
    wavelength_diff=1;
    wavelength_error=zeros(iteration_num,1);
    wavelength_vec=zeros(iteration_num,1);
    
    while (iteration_count < iteration_num)&&(wavelength_diff > iteration_wavelength_diff)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [fiber_modes,fem_fib]=fiber_fem_solver(fiber_res);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fiber_res.mode_count=length(fiber_modes.Neff);
        
        res_lambda0=fiber_res.lambda;
        res_M_mat=calc_L_mat(fiber_res,fiber_modes,fem_fib,res_lambda0,res_m,scale_fac);
        res_A_mat=-diag((2*pi*fiber_modes.Neff).^2);
        res_T_mat=inv(res_A_mat)*res_M_mat;
        [V_matrix,D_mat]=eig(res_T_mat);
        % [V_mat,D_mat]=eigs(res_T_mat,1,'SR');  % 'LM' or 'SM' - Largest or Smallest Magnitude, 'LR' or 'SR' - Largest or Smallest Real part
        % D_vec=diag(D_mat);
        D_array=real(diag(D_mat));
        [sol,sol_ind]=sort(D_array,'ascend');
        D_indx=sol_ind(mode_order_count);
        D_vec=D_array(D_indx);
        V_mat=V_matrix(:,D_indx);
        
        D_dx=D_vec;
        res_lambda = abs(real(res_lambda0*(1/(D_dx*res_lambda0^2 + 1))^(1/2))); % negative solution ingored
        
        iteration_count=iteration_count+1;
        wavelength_diff=abs(res_lambda0-res_lambda);
        fiber_res.lambda=res_lambda;
        wavelength_error(iteration_count)=wavelength_diff;
        wavelength_vec(iteration_count)=res_lambda;
        
    end % while
    
    res_lambda_record(mode_order_count)=res_lambda;
    
    %     %======================================================================
    %     % plot iteration_vs_wavelength
    %     x_axis=(1:iteration_count);
    %     figure;clf;
    %     [ax,h1,h2]=plotyy(x_axis,wavelength_vec(1:iteration_count)*1e9,x_axis,wavelength_error(1:iteration_count)*1e9);
    %     set(ax(2),'YScale','log');
    %     title(['m = ' num2str(res_m) '   mode order = ' num2str(mode_order_count)]);
    %     xlabel('iteration number');
    %     ylabel(ax(1),'\lambda (nm)');
    %     ylabel(ax(2),'\Delta \lambda (nm)');
    %     %======================================================================
    
    iteratioin{mode_order_count}.lambda_record=wavelength_vec(1:iteration_count);
    iteratioin{mode_order_count}.lambda_diff_record=wavelength_error(1:iteration_count);
    
    root_num=length(D_vec);
    Field_temp=zeros(6,size(fem_fib.mesh.p,2));
    
    
    V=V_mat(:,root_num);
    Resonance{mode_order_count}.lambda=res_lambda;
    Resonance{mode_order_count}.D=D_vec(root_num);
    Resonance{mode_order_count}.m=res_m;
    Resonance{mode_order_count}.V=V;
    Resonance{mode_order_count}.Field=zeros(size(Field_temp));
    
    for count=1:size(V,1);
        for count_field=1:length(fiber_comp_str)
            Field_temp(count_field,:)=postinterp(fem_fib,fiber_comp_str{count_field},fem_fib.mesh.p,'solnum',fiber_modes.Neff_indx(count));
        end %for count_field
        Resonance{mode_order_count}.Field=Resonance{mode_order_count}.Field+V(count)*Field_temp;
    end %for count
    
    fiber_fields_bend=Resonance{mode_order_count}.Field;
    file_label='Method_one_toroid_';
    fig_id=100+mode_order_count;
    figid=plotfields_mode(fem_fib,fiber_fields_bend,fig_id,file_label);
    
    
end % mode_order_count

time=toc;



