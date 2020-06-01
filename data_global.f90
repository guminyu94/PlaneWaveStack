module data_global
    use Sim_parameters, only : wp
    implicit none
    real, allocatable :: data_1(:)
    real, allocatable :: data_2(:,:) 
    real, allocatable :: data_3(:,:) 
    integer :: counter
end module data_global
    