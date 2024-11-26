      program writeParameter

      open(unit=20,file="lake-params.txt",status="old")  ! latin hypercube scalings
      open(unit=22,file="row.txt",status="old")     ! row of LHS file to use, from shell script

      read(22,*) inum

      do j=1,inum

        read(20,*) v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18  ! match number of parameters in lake-params.txt

      end do

! CHANGE following lines to specify ranges for lake.inc parameters you wish to vary
      v1 = 2.e-3*v1 + 1.e-3   ! cdrn from 1 to 3 e-3
      v2 = 0.5*v2 + 0.2       ! eta from 0.2 to 0.7
      v3 = 0.09*v3 + 0.7       ! alb_snow from 0.7 to 0.79
      v4 = 0.15*v4 + 0.4       ! alb_slush from 0.4 to 0.55
      v5 = 2.e6*v5 + 2.e6     ! csed from 2e6 to 4e6
      v6 = 2.*v6 + 0.5       ! condsed from 0.5 to 2.5
      v7 = 0.15*v7 + 0.05     ! albsed from 0.05 to 0.2
      v8 = 24*v8 - 42.1     ! d18Oa from -42.1 to -18.1 -- use ranges from Hannah Bailey dataset
      v9 = 188.9*v9 - 322.8   ! d2Ha from -322.8 to -133.9
      v10 = 0.9* v10 + 0.0      ! f from 0 to 1  (changed to 0.9 to avoid non-realistic isotope behavior)

! CHANGE following lines to specify ranges for met-generate.py parameters you wish to vary
      v11 = 0.95*v11 + 0.05     ! melt_ratio from 0.05? to 1
      v12 = 0.95*v12 + 0.05     ! rp_ratio_summer from 0.05 to 1
      v13 = 0.95*v13 + 0.05     ! rp_ratio_winter from 0.05 to 1
      v14 = 0.95*v14 + 0.05     ! rsm_ratio from 0.05 to 1
      v15 = 100*v15 + 0         ! p
      v16 = 10*v16 + 0          ! S
      v17 = 20 * v17            ! spring threshold of days above freezing from 0 to 20
      v18 = 20 * v18            ! fall threshold of days above freezing from 0 to 20    

! CHANGE following lines to write lake.inc parameter statements for parameters you wish to vary
      write(23,*) "      parameter (cdrn = ", v1,")"
      write(23,*) "      parameter (eta = ", v2,")"
      write(23,*) "      parameter (alb_snow = ", v3,")"
      write(23,*) "      parameter (alb_slush = ", v4,")"
      write(23,*) "      parameter (csed = ", v5,")"
      write(23,*) "      parameter (condsed = ", v6,")"
      write(23,*) "      parameter (alb_sed = ", v7,")"
      write(23,*) "      parameter (d18Oa = ", v8,")"
      write(23,*) "      parameter (d2Ha = ", v9,")"
      write(23,*) "      parameter (f = ", v10,")"
      close(23)

! CHANGE following lines to write met-generate parameter statements for parameters you wish to vary
      write(24,"(a,f16.8)") "melt_ratio = ", v11
      write(24,"(a,f16.8)") "rp_ratio_summer = ", v12
      write(24,"(a,f16.8)") "rp_ratio_winter = ", v13
      write(24,"(a,f16.8)") "rp_ratio_winter = ", v13
      write(24,"(a,f16.8)") "rsm_ratio = ", v14
      write(24,"(a,f16.8)") "p = ", v15
      write(24,"(a,f16.8)") "s = ", v16
      write(24,"(a,f16.8)") "thresh_spring = ", v17
      write(24,"(a,f16.8)") "thresh_fall = ", v18
      close(24)


      end
