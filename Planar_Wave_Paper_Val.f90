Module Planar_Wave_Paper_Val
    use Sim_parameters
    use Swapper
    use graphene
    use Plot_Pgplot
    use Data_Loader
    implicit none    
    
    contains
    subroutine jiang_case
        use Modified_Otto
        implicit none 
        type(dots), allocatable :: p_x(:), p_y(:)
        procedure(fun_config), pointer :: fun_p
        fun_p => mod_otto_config
        
        ! graphene parameters
        muc = 1.2
        tau = 0.5e-12
        b0 = 0.0
        
        allocate(data_1(1001)) 
        allocate(data_2(3,1001))
        allocate(data_3(3,1001))
        
        allocate(p_x(3))
        allocate(p_y(3))
        
        call theta_swap(fun_p,5.1742e12_wp,0.0_wp, 89.999_wp, 1001)
        call csv_data_1d('jiang_fig2_t.csv',p_x(1),p_y(1))
        call csv_data_1d('jiang_fig2_r.csv',p_x(2),p_y(2))
        call csv_data_1d('jiang_fig2_a.csv',p_x(3),p_y(3))
        
        call plot_1d(data_1,data_2, dots_x = p_x,dots_y = p_y, x_label = '\(0685) (Degrees)', y_label = 'Amplitude (A.U.)', title = '', color = (/1,2,3/), style = (/1,2,3/),dev = 'case1.ps/CPS',legend=(/'Transmittance','Reflectance','Absorption'/),dots_legend=(/'Transmitance (Ref.)','Reflectance (Ref.)','Absorption (Ref.)'/),dots_style=(/4,5,6/),yrange=(/0.0,1.0/))

    end subroutine jiang_case
    
    subroutine da_case
        use Graphene_PC_OM
        use data_global
        implicit none 
        type(dots), allocatable :: p_x(:), p_y(:), p_x_2(:), p_y_2(:)
        procedure(fun_config), pointer :: fun_p
        fun_p => graphene_pc_config
        
        counter = 1
        allocate(data_1(51)) 
        allocate(data_2(3,51))
        allocate(data_3(3,51))
        allocate(p_x(3))
        allocate(p_y(3))
        allocate(p_x_2(1))
        allocate(p_y_2(1))
        ! graphene parameters
        muc = 8.621738e-5 * 50
        tau = 1.6-12
        b0 = 1.0
        
        call freq_swap(fun_p,23e12_wp,26e12_wp, 51)
        
        ! graphene parameters
        muc = 8.621738e-5 * 510
        call freq_swap(fun_p,23e12_wp,26e12_wp, 51)
        
        ! graphene parameters
        muc = 8.621738e-5 * 660
        call freq_swap(fun_p,23e12_wp,26e12_wp, 51)
        
        
        call csv_data_1d('da_fig3_50k.csv',p_x(1),p_y(1))
        call csv_data_1d('da_fig3_510k.csv',p_x(2),p_y(2))
        call csv_data_1d('da_fig3_660k.csv',p_x(3),p_y(3))
        
        call csv_data_1d('da_fig3_t.csv',p_x_2(1),p_y_2(1))
        
        call plot_1d(data_1, data_2, dots_x =  p_x, dots_y = p_y, x_label = '\(2156) (THz)', y_label = 'Angle (Degrees) ', title = '', dev='jiang_case_2_kerrangle.ps/CPS', color = (/1,2,3/), style = (/1,2,3/), legend = (/'Kerr Angle, \(2138) = 50K','Kerr Angle, \(2138) = 510K','Kerr Angle, \(2138) = 660K'/),dots_legend = (/'Kerr Angle, \(2138) = 50K, Ref.','Kerr Angle, \(2138) = 510K, Ref.','Kerr Angle, \(2138) = 660K, Ref.'/),dots_style=(/4,5,6/),yrange=(/-0.6,0.3/))
        call plot_1d(data_1, data_3, dots_x = p_x_2, dots_y = p_y_2, x_label = '\(2156) (THz)', y_label = 'Amplitude (A.U.) ', title = '', dev='jiang_case_2_tc.ps/CPS', color = (/1,2,3/), style = (/1,2,3/), legend = (/'Transmittance, \(2138) = 50K','Transmittance, \(2138) = 510K','Transmittance, \(2138) = 660K'/),dots_legend = (/'Transmittance, Ref.'/),dots_style=(/4/),yrange=(/0.0,1.0/))
    end subroutine  da_case
         
end Module Planar_Wave_Paper_Val