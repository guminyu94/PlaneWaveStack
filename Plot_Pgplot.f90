Module Plot_Pgplot
    contains
    subroutine plot_1d(x,y,x_label,y_label,title)
        real, allocatable, intent(in) :: x(:),y(:) 
        character(len=*), intent(in), optional  :: x_label, y_label, title
        integer :: n_p
        n_p = size(x)
        !   initialize the pgplot, and select windows as plotting device     
        IER = pgbeg(0,'/WZ',1,1)
        
        ! adjust x and y lim automatically
        !   start a plot, xlim, ylim, box style, and x y scale    
        CALL pgenv(minval(x),maxval(x),minval(y),maxval(y),0,1)
        !   labeling
        CALL pglab(x_label, y_label, title)
        !   draw point mark
        !     CALL PGPT(5,XS,YS,9)

        ! plot line, num of data points, data
        CALL PGLINE(n_p,x,y)
        CALL pgend
    end subroutine plot_1d
end module Plot_Pgplot