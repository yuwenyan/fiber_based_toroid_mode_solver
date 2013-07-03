% toroid mode solver based on graded fiber modes
% graded index fiber modes
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
L_ring=2*pi*fiber_geom.toroid_R;
val_m=fiber_geom.M_guess;
num_toroid_mode=6;
scale_fac=1;

fiber_comp_str={'Ex','Ey','Ez','Hx','Hy','Hz'};
toroid_comp_str={'Erad','Eaxi','Eazi','Hrad','Haxi','Hazi'};

m_result=val_m; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode_calculation=num_toroid_mode;
Resonance=cell(1,mode_calculation);
% Resonance.Field=[];
% Resonance.Field2=[];
% Resonance.m=[];

for mode_order_count=1:mode_calculation
    
    max_iteration=30;
    count_iteration=0;
    limit_m_diff=0.00001;
    limit_lambda_diff=0.002e-9;
    m_diff=1;
    lambda_diff=1;
    res_lambda_vec=zeros(1,max_iteration);
    res_lambda_error=zeros(1,max_iteration);
    res_m_error=zeros(1,max_iteration);
    
    fiber_res=fiber_geom;
    
    while (count_iteration<max_iteration)&&(lambda_diff>limit_lambda_diff)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [fiber_modes,fem_fib]=fiber_fem_solver_graded(fiber_res);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Neff_indx=fiber_modes.Neff_indx(mode_order_count);
        Neff=fiber_modes.Neff2(Neff_indx);
        m=L_ring*Neff/fiber_res.lambda;
        res_m=m_result;
        
        res_lambda_temp = fiber_res.lambda*m/res_m;
        lambda_diff= abs(res_lambda_temp-fiber_res.lambda);
        m_diff=abs(m-res_m);
        %             disp(['deviation: ' int2str(deviation)]);
        res_lambda_vec(count_iteration+1)= res_lambda_temp;
        res_lambda_error(count_iteration+1)=lambda_diff;
        res_m_error(count_iteration+1)=m_diff;
        
        fiber_res.lambda=res_lambda_temp;
        count_iteration=count_iteration+1;
        
    end % while
    
    if (count_iteration==max_iteration);
        disp('Maxium iterations!');
    end
    
%     x_axis=(1:count_iteration);
%     figure;
%     [ax,h1,h2]=plotyy(x_axis,res_lambda_vec(1:count_iteration),x_axis,res_lambda_error(1:count_iteration));
%     set(ax(2),'YScale','log');
%     ylabel(ax(2),'lambda error');
%     figure;
%     [m_ax,m_h1,m_h2]=plotyy(x_axis,res_lambda_vec(1:count_iteration),x_axis,res_m_error(1:count_iteration));
%     set(m_ax(2),'YScale','log');
%     ylabel(m_ax(2),'m error');
    
    Resonance{mode_order_count}.m=res_m;
    Field_temp=zeros(6,size(fem_fib.mesh.p,2));    
    Resonance{mode_order_count}.Field=zeros(size(Field_temp));
    
    for count_field=1:length(fiber_comp_str)
        Field_temp(count_field,:)=postinterp(fem_fib,fiber_comp_str{count_field},fem_fib.mesh.p,'solnum',Neff_indx);
    end
    
    Resonance{mode_order_count}.Field=Field_temp;
    fiber_fields_bend=Resonance{mode_order_count}.Field;

    file_label='toroid';
    fig_id=100+mode_order_count;
    figid=plotfields_mode(fem_fib,fiber_fields_bend,fig_id,file_label);
    
end % mode_order_count

time=toc;
