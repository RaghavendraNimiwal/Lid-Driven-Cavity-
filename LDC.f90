program LDC

  implicit none

  integer, parameter :: NROWS = 32
  integer, parameter :: NCOLS = 32

  real, parameter :: Re = 100.0
  real, parameter :: tolerance= 0.0000001
  real, parameter :: delt = 0.003125

  integer :: iter = 0
  integer :: iter_p = 0
  integer :: r,c

  real :: delx = 1/real(NCOLS)
  real :: dely = 1/real(NROWS)

  real :: avg_L2_u = 1000.0
  real :: avg_L2_v = 1000.0

  real :: residue(NROWS+2,NCOLS+2)
  real :: avg_L2_p, ae, aw, an, as, ap, Source

  real :: Adv, Dif
  real :: ue,uw,un,us,ve,vw,vn,vs

  real :: relaxation = 1.0

  real :: u_vel(NROWS+2,NCOLS+1),u_vel_star(NROWS+2,NCOLS+1),u_vel_old(NROWS+2,NCOLS+1),&
          v_vel(NROWS+1,NCOLS+2),v_vel_old(NROWS+1,NCOLS+2),v_vel_star(NROWS+1,NCOLS+2),pres(NROWS+2,NCOLS+2)
  real :: u_vel_p(NROWS,NROWS), v_vel_p(NROWS,NROWS)

  do r=1,NROWS+2
    do c=1,NCOLS+2
      pres(r,c) = 0
    enddo
  enddo
  do r=1,NROWS+2
    do c=1,NCOLS+1
      u_vel(r,c) = 0
      u_vel_old(r,c) = 0
      u_vel_star(r,c) = 0
    enddo
  enddo
  do r=1,NROWS+1
   do c=1,NCOLS+2
      v_vel(r,c) = 0
      v_vel_old(r,c) = 0
      v_vel_star(r,c) = 0
    enddo
  enddo

  open(4, file='AVG_L2_a.dat', status='unknown')

  do while (avg_L2_u .gt. tolerance .or. avg_L2_v .gt. tolerance)

    !UPDATING BOUNDARY VALUES
    !boundary condition for v velocity
    do c=1,NCOLS+2
      v_vel(NROWS+1,c) = 0
      v_vel(1,c) = 0
    enddo
    do r=1,NROWS+1
      v_vel(r,1) = 2*0 - v_vel(r,2)
      v_vel(r,NCOLS+2) = 2*0 - v_vel(r,NCOLS+1)
    enddo
    !boundary condition for u velocity
    do r=1,NROWS+2
      u_vel(r,1) = 0
      u_vel(r,NCOLS+1) = 0
    enddo
    do c=1,NCOLS+1
      u_vel(1,c) = 2*0 - u_vel(2,c)
      u_vel(NROWS+2,c) = 2*1 - u_vel(NROWS+1,c)
    enddo

    ! CALCULATING THE VELOCITIES AT HALF STEP
    if(iter == 0) then
      do r=2,NROWS+1
        do c=2,NCOLS
          ue = (u_vel(r,c) + u_vel(r,c+1))/2
          uw = (u_vel(r,c) + u_vel(r,c-1))/2
          us = (u_vel(r,c) + u_vel(r-1,c))/2
          un = (u_vel(r,c) + u_vel(r+1,c))/2
          vn = (v_vel(r,c) + v_vel(r,c+1))/2
          vs = (v_vel(r-1,c) + v_vel(r-1,c+1))/2
          Adv = -(ue*ue - uw*uw)/delx - (un*vn - us*vs)/dely
          Dif = (u_vel(r+1,c)-2*u_vel(r,c)+u_vel(r-1,c))/(dely*dely*Re) +&
                (u_vel(r,c+1)-2*u_vel(r,c)+u_vel(r,c-1))/(delx*delx*Re)
          u_vel_star(r,c) = u_vel(r,c) + delt*Adv + delt*Dif
        enddo
      enddo
      do r=2,NROWS
        do c=2,NCOLS+1
          ve = (v_vel(r,c) + v_vel(r,c+1))/2
          vw = (v_vel(r,c) + v_vel(r,c-1))/2
          ue = (u_vel(r+1,c) + u_vel(r,c))/2
          uw = (u_vel(r+1,c-1) + u_vel(r,c-1))/2
          vn = (v_vel(r+1,c) + v_vel(r,c))/2
          vs = (v_vel(r-1,c) + v_vel(r,c))/2
          Adv = -(vn*vn - vs*vs)/dely - (ve*ue - vw*uw)/delx
          Dif = (v_vel(r+1,c)-2*v_vel(r,c)+v_vel(r-1,c))/(dely*dely*Re) + &
                (v_vel(r,c+1)-2*v_vel(r,c)+v_vel(r,c-1))/(delx*delx*Re)
          v_vel_star(r,c) = v_vel(r,c) + delt*Adv + delt*Dif
        enddo
      enddo

    else

      do r=2,NROWS+1
        do c=2,NCOLS
          ue = (u_vel(r,c) + u_vel(r,c+1))/2
          uw = (u_vel(r,c) + u_vel(r,c-1))/2
          us = (u_vel(r,c) + u_vel(r-1,c))/2
          un = (u_vel(r,c) + u_vel(r+1,c))/2
          vn = (v_vel(r,c) + v_vel(r,c+1))/2
          vs = (v_vel(r-1,c) + v_vel(r-1,c+1))/2
          Adv = -(ue*ue - uw*uw)/delx - (un*vn - us*vs)/dely
          Dif = (u_vel(r+1,c)-2*u_vel(r,c)+u_vel(r-1,c))/(dely*dely*Re) +&
                (u_vel(r,c+1)-2*u_vel(r,c)+u_vel(r,c-1))/(delx*delx*Re)
          u_vel_star(r,c) = u_vel(r,c) + (3/2)*delt*Adv + (3/2)*delt*Dif
        enddo
      enddo
      do r=2,NROWS+1
        do c=2,NCOLS
          ue = (u_vel_old(r,c) + u_vel_old(r,c+1))/2
          uw = (u_vel_old(r,c) + u_vel_old(r,c-1))/2
          us = (u_vel_old(r,c) + u_vel_old(r-1,c))/2
          un = (u_vel_old(r,c) + u_vel_old(r+1,c))/2
          vn = (v_vel_old(r,c) + v_vel_old(r,c+1))/2
          vs = (v_vel_old(r-1,c) + v_vel_old(r-1,c+1))/2
          Adv = -(ue*ue - uw*uw)/delx - (un*vn - us*vs)/dely
          Dif = (u_vel_old(r+1,c)-2*u_vel_old(r,c)+u_vel_old(r-1,c))/(dely*dely*Re) +&
                (u_vel_old(r,c+1)-2*u_vel_old(r,c)+u_vel_old(r,c-1))/(delx*delx*Re)
          u_vel_star(r,c) = u_vel_star(r,c) - (1/2)*delt*Adv - (1/2)*delt*Dif
        enddo
      enddo

      do r=2,NROWS
        do c=2,NCOLS+1
          ve = (v_vel(r,c) + v_vel(r,c+1))/2
          vw = (v_vel(r,c) + v_vel(r,c-1))/2
          ue = (u_vel(r+1,c) + u_vel(r,c))/2
          uw = (u_vel(r+1,c-1) + u_vel(r,c-1))/2
          vn = (v_vel(r+1,c) + v_vel(r,c))/2
          vs = (v_vel(r-1,c) + v_vel(r,c))/2
          Adv = -(vn*vn - vs*vs)/dely - (ve*ue - vw*uw)/delx
          Dif = (v_vel(r+1,c)-2*v_vel(r,c)+v_vel(r-1,c))/(dely*dely*Re) + &
                (v_vel(r,c+1)-2*v_vel(r,c)+v_vel(r,c-1))/(delx*delx*Re)
          v_vel_star(r,c) = v_vel(r,c) + (3/2)*delt*Adv + (3/2)*delt*Dif
        enddo
      enddo
      do r=2,NROWS
        do c=2,NCOLS+1
          ve = (v_vel_old(r,c) + v_vel_old(r,c+1))/2
          vw = (v_vel_old(r,c) + v_vel_old(r,c-1))/2
          ue = (u_vel_old(r+1,c) + u_vel_old(r,c))/2
          uw = (u_vel_old(r+1,c-1) + u_vel_old(r,c-1))/2
          vn = (v_vel_old(r+1,c) + v_vel_old(r,c))/2
          vs = (v_vel_old(r-1,c) + v_vel_old(r,c))/2
          Adv = -(vn*vn - vs*vs)/dely - (ve*ue - vw*uw)/delx
          Dif = (v_vel_old(r+1,c)-2*v_vel_old(r,c)+v_vel_old(r-1,c))/(dely*dely*Re) + &
                (v_vel_old(r,c+1)-2*v_vel_old(r,c)+v_vel_old(r,c-1))/(delx*delx*Re)
          v_vel_star(r,c) = v_vel_star(r,c) - (1/2)*delt*Adv - (1/2)*delt*Dif
        enddo
      enddo

    endif


    !PRESSURE POISSON EQUATION

    avg_L2_p = 1000.0
    iter_p = 0

    do while (avg_L2_p .gt. 0.00001)
      do r=2,NROWS+1
        do c=2,NCOLS+1
          ae = -delt*dely/delx
          aw = -delt*dely/delx
          an = -delt*delx/dely
          as = -delt*delx/dely
          ap = -(ae + aw + an + as)
          Source = (u_vel_star(r,c-1) - u_vel_star(r,c))*dely + (v_vel_star(r-1,c) - v_vel_star(r,c))*delx
          pres(r,c) = relaxation*((-ae*pres(r,c+1) -aw*pres(r,c-1) -an*pres(r+1,c) -as*pres(r-1,c) + Source)/ap)&
                      + (1-relaxation)*pres(r,c)
        enddo
      enddo

      !updating pressure boundary conditions
      do r=2,NROWS+1
        pres(r,1) = pres(r,2)
        pres(r,NCOLS+2) = pres(r,NCOLS+1)
      enddo
      do c=2,NCOLS+1
        pres(1,c) = pres(2,c)
        pres(NROWS+2,c) = pres(NROWS+1,c)
      enddo
        !calculating average L2 norm
      avg_L2_p = 0
      do r=2,NROWS+1
        do c=2,NCOLS+1
          ae = -delt*dely/delx
          aw = -delt*dely/delx
          an = -delt*delx/dely
          as = -delt*delx/dely
          ap = -(ae + aw + an + as)
          Source = (u_vel_star(r,c-1) - u_vel_star(r,c))*dely + (v_vel_star(r-1,c) - v_vel_star(r,c))*delx
          residue(r,c) = Source -ap*pres(r,c) -ae*pres(r,c+1) -aw*pres(r,c-1) -an*pres(r+1,c) -as*pres(r-1,c)
            avg_L2_p = avg_L2_p + residue(r,c)*residue(r,c)
        enddo
      enddo
        avg_L2_p = sqrt(avg_L2_p/(NROWS*NCOLS))
        iter_p = iter_p + 1
        if (iter_p .gt. 10000000) then
        print*, "*****************************************"
        print*, "ERROR IN PRESSURE POISSON SOLVER"
        print*, "*****************************************"
      endif
    enddo !END OF PRESSURE LOOP

    !CORRECTED VELOCITY
    do r=2,NROWS+1
      do c=2,NCOLS
        u_vel(r,c) = u_vel_star(r,c) - (pres(r,c+1) - pres(r,c))*delt/delx
      enddo
    enddo
    do r=2,NROWS
      do c=2,NCOLS+1
        v_vel(r,c) = v_vel_star(r,c) - (pres(r+1,c) - pres(r,c))*delt/dely
      enddo
    enddo

    !Average L2 norm for velocities
    avg_L2_u = 0
    do r=2,NROWS+1
      do c=1,NCOLS+1
        avg_L2_u = avg_L2_u + (u_vel(r,c) - u_vel_old(r,c))*(u_vel(r,c) - u_vel_old(r,c))
      enddo
    enddo
    avg_L2_u = sqrt(avg_L2_u/(NROWS*(NCOLS+1)))

    avg_L2_v = 0
    do r=1,NROWS+1
      do c=2,NCOLS+1
        avg_L2_v = avg_L2_v + (v_vel(r,c) - v_vel_old(r,c))*(v_vel(r,c) - v_vel_old(r,c))
      enddo
    enddo
    avg_L2_v = sqrt(avg_L2_v/((NROWS+1)*NCOLS))

    !UPDATING OLD VELOCITIES
    do r=1,NROWS+2
      do c=1,NCOLS+1
        u_vel_old(r,c) = u_vel(r,c)
      enddo
    enddo
    do r=1,NROWS+1
      do c=1,NCOLS+2
        v_vel_old(r,c) = v_vel(r,c)
      enddo
    enddo

    iter = iter+1

    write(4 ,* ) avg_L2_u, avg_L2_v, delt*iter

    if(iter .gt. 10000000) then
      print*, "*****************************************"
      print*, "ERROR IN OUTER ITERATION"
      print*, "*****************************************"
    endif

    print*,"ITER",iter,"L2u",avg_L2_u,"L2v",avg_L2_v

  enddo !END OF OUTER LOOP

  close(4)

  do r = 1,NROWS
    do c = 1,NCOLS
      u_vel_p(r,c) = (u_vel(r+1,c) + u_vel(r+1,c+1))/2
    enddo
  end do
  do r = 1,NROWS
    do c = 1,NCOLS
      v_vel_p(r,c) = (v_vel(r,c+1) + v_vel(r+1,c+1))/2
    enddo
  end do


  open(5, file ='Mid_u64.dat', status ='unknown')
  do r = 2,NROWS+1
    write(5 ,* ) u_vel(r,(NROWS/2)+1)
  end do
  close(5)

  open(6, file ='Mid_v64.dat', status ='unknown')
  do c = 2,NCOLS+1
      write(6 ,* ) v_vel((NCOLS/2)+1,c)
  end do
  close(6)

  open(1, file='Velocity_x.dat', status='unknown')
  do r = 1,NROWS
    do c = 1,NCOLS
      write(1 ,* ) u_vel_p(r,c)
    end do
  end do
  close(1)

  open(2, file='Velocity_y.dat', status='unknown')
  do r = 1,NROWS
    do c = 1,NCOLS
      write(2 ,* ) v_vel_p(r,c)
    end do
  end do
  close(2)

  open(3, file ='Pressure.dat', status ='unknown')
  do r = 2,NROWS+1
    do c = 2,NCOLS+1
      write(3 ,* ) pres(r,c)
    end do
  end do
  close(3)


end program LDC
