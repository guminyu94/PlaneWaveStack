!****************************************************************************
!
!   Module: Data_Loader
!
!   PURPOSE: A subroutine wrapper for loading the data
!
!****************************************************************************
Module Data_Loader
    
    contains

    subroutine csv_data_1d(file_name,x,y)
        use Plot_Pgplot, only : dots
    !*****************************************************************************
    ! load the 1d x, y data from cvs
        character(len=*), intent(in), optional :: file_name
        type(dots) :: x, y
        integer :: file_len, file_index, i
        
        ! open the file
        call csv_file_open_read ( file_name, file_index )
        ! count the line
        call csv_file_line_count(file_name,file_len)
        
        if (allocated(x%dots_data)) then
            deallocate(x%dots_data)
            allocate(x%dots_data(file_len))
        else
            allocate(x%dots_data(file_len))
        end if
        
        if (allocated(y%dots_data)) then
            deallocate(y%dots_data)
            allocate(y%dots_data(file_len))
        else
            allocate(y%dots_data(file_len))
        end if
        
        do i=1,file_len
            read (file_index,*) x%dots_data(i),y%dots_data(i)
        end do
        
        ! close 
        call csv_file_close_read(file_index)
        
    end subroutine csv_data_1d
    
    subroutine csv_value_count ( csv_record, csv_record_status, value_count )

!*****************************************************************************80
!
!! CSV_COUNT counts the number of values in a CSV record.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  character csv_char
  character csv_char_old
  integer ( kind = 4 ) csv_len
  integer ( kind = 4 ) csv_loc
  character ( len = * ) csv_record
  integer ( kind = 4 ) csv_record_status
  character :: TAB = achar ( 9 )
  integer ( kind = 4 ) value_count
  integer ( kind = 4 ) word_length
!
!  No values so far.
!
  value_count = 0
!
!  We begin in "unquoted" status.
!
  csv_record_status = 0
!
!  How many characters in the record?
!
  csv_len = len_trim ( csv_record )
!
!  Count number of characters in each word.
!
  word_length = 0
!
!  Consider each character.
!
  csv_char_old = ','

  do csv_loc = 1, csv_len

    csv_char = csv_record(csv_loc:csv_loc)
!
!  Each comma divides one value from another.
!
    if ( csv_char_old == ',' ) then

      value_count = value_count + 1
      word_length = 0
!
!  For quotes, try using CSV_RECORD_STATUS to count the number of
!  quoted characters.
!
    else if ( csv_char == '"' ) then

      if ( 0 < csv_record_status ) then
        csv_record_status = 0
      else
        csv_record_status = csv_record_status + 1
      end if
!
!  Ignore blanks
!
    else if ( csv_char == ' ' .or. csv_char == TAB ) then
!
!  Add character to length of word.
!
    else

      word_length = word_length + 1

      if ( value_count == 0 ) then
        value_count = 1
      end if

    end if

    csv_char_old = csv_char

  end do

  return
    end
    
subroutine csv_file_close_read( csv_file_unit )

!*****************************************************************************80
!
!! CSV_FILE_CLOSE_READ closes a CSV file for reading.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
!
  implicit none

  integer ( kind = 4 ) csv_file_unit

  close ( unit = csv_file_unit )

  return
end

subroutine csv_file_close_write ( csv_file_unit )

!*****************************************************************************80
!
!! CSV_FILE_CLOSE_WRITE closes a CSV file for writing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
!
  implicit none

  integer ( kind = 4 ) csv_file_unit

  close ( unit = csv_file_unit )

  return
end

subroutine csv_file_header_write ( csv_file_unit, header )

!*****************************************************************************80
!
!! CSV_FILE_HEADER_WRITE writes a header to a CSV file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
!
!    Input, character ( len = * ) HEADER, the header.
!
  implicit none

  integer ( kind = 4 ) csv_file_unit
  character ( len = * ) header

  write ( csv_file_unit, '(a)' ) trim ( header )

  return
end

subroutine csv_file_line_count ( csv_file_name, line_num )

!*****************************************************************************80
!
!! CSV_FILE_LINE_COUNT counts the number of lines in a CSV file.
!
!  Discussion:
!
!    This routine does not try to distinguish the possible header line,
!    blank lines, or cases where a single CSV record extends over multiple 
!    lines.  It simply counts the number of lines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CSV_FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of lines.
!
  implicit none

  character ( len = * ) csv_file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  character ( len = 1023 ) line
  integer ( kind = 4 ) line_num

  line_num = -1

  call get_unit ( input_unit )

  open ( unit = input_unit, file = csv_file_name, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSV_FILE_LINE_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open "' // trim ( csv_file_name ) // '".'
    stop
  end if

  line_num = 0

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = line_num
      exit
    end if

    line_num = line_num + 1

  end do

  close ( unit = input_unit )

  return
end

subroutine csv_file_record_write ( csv_file_unit, record )

!*****************************************************************************80
!
!! CSV_FILE_RECORD_WRITE writes a record to a CSV file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
!
!    Input, character ( len = * ) RECORD, the record.
!
  implicit none

  integer ( kind = 4 ) csv_file_unit
  character ( len = * ) record

  write ( csv_file_unit, '(a)' ) trim ( record )

  return
end

subroutine csv_file_open_read ( csv_file_name, csv_file_unit )

!*****************************************************************************80
!
!! CSV_FILE_OPEN_READ opens a CSV file for reading.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CSV_FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
!
  implicit none

  character ( len = * ) csv_file_name
  integer ( kind = 4 ) csv_file_status
  integer ( kind = 4 ) csv_file_unit

  call get_unit ( csv_file_unit )

  open ( unit = csv_file_unit, file = csv_file_name, status = 'old', &
    iostat = csv_file_status )

  if ( csv_file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSV_FILE_OPEN_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open "' // trim ( csv_file_name ) // '".'
    csv_file_unit = - 1
    stop
  end if

  return
end

subroutine csv_file_open_write ( csv_file_name, csv_file_unit )

!*****************************************************************************80
!
!! CSV_FILE_OPEN_WRITE opens a CSV file for writing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CSV_FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) CSV_FILE_UNIT, the unit number
!
  implicit none

  character ( len = * ) csv_file_name
  integer ( kind = 4 ) csv_file_status
  integer ( kind = 4 ) csv_file_unit

  call get_unit ( csv_file_unit )

  open ( unit = csv_file_unit, file = csv_file_name, status = 'replace', &
    iostat = csv_file_status )

  if ( csv_file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSV_FILE_OPEN_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open "' // trim ( csv_file_name ) // '".'
    stop
  end if

  return
end



subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
function i4_log_10 ( i )

!*****************************************************************************80
!
!! I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I  I4_LOG_10
!    -----  --------
!        0    0
!        1    0
!        2    0
!        9    0
!       10    1
!       11    1
!       99    1
!      100    2
!      101    2
!      999    2
!     1000    3
!     1001    3
!     9999    3
!    10000    4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 10
!    is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_10, the integer part of the
!    logarithm base 10 of the absolute value of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) ten_pow

  if ( i == 0 ) then

    i4_log_10 = 0

  else

    i4_log_10 = 0
    ten_pow = 10

    i_abs = abs ( i )

    do while ( ten_pow <= i_abs )
      i4_log_10 = i4_log_10 + 1
      ten_pow = ten_pow * 10
    end do

  end if

  return
end


end module Data_Loader