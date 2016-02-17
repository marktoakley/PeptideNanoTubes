! Calculate properties of molecules from PATHSAMPLE/AMBER runs
! Works for any cyclic peptide where all residues are the same type
program metric
   use hb
   implicit none
   integer::  natoms
   integer:: res_size ! Number of atoms/residue
   integer:: res_coor !Number of coordinates per residue
   double precision, allocatable :: coords(:)
   integer recno
   integer i,j
   character*50 arg,filename
! Set up stuff (assuming only Ala for now)
   if (COMMAND_ARGUMENT_COUNT().lt.3) then
      write(*,*) "usage: PSHBond <int residues> <int rings> <int res_size>"
      stop
   endif
   call GET_COMMAND_ARGUMENT(1,arg)
   read(arg,*)residues
   call GET_COMMAND_ARGUMENT(2,arg)
   read(arg,*)rings
   call GET_COMMAND_ARGUMENT(3,arg)
   read(arg,*)res_size
   natoms=rings*residues*res_size
   res_coor=3*res_size
   allocate(coords(3*natoms))

! Open output file
   open(30,file="hbonds.csv")
! Open PATHSAMPLE minima file
   open(10,file="points.min",ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=8*3*NATOMS)
   recno=0
! Start reading pathsample file here
150 recno=recno+1
   read(10,rec=recno,err=200) (coords(j),j=1,3*natoms)
! Copy coords to hbond array
! Atom indices for Ala N 1, H 2, C 9 O 10
   do i=1,rings*residues
      !N
      peptides(i,3,1)=coords((res_coor*(i-1))+1)
      peptides(i,3,2)=coords((res_coor*(i-1))+2)
      peptides(i,3,3)=coords((res_coor*(i-1))+3)
      !H
      peptides(i,4,1)=coords((res_coor*(i-1))+4)
      peptides(i,4,2)=coords((res_coor*(i-1))+5)
      peptides(i,4,3)=coords((res_coor*(i-1))+6)
      if (mod(i,residues).eq.0) then
         !C
         peptides(i-residues+1,2,1)=coords((res_coor*(i))-5)
         peptides(i-residues+1,2,2)=coords((res_coor*(i))-4)
         peptides(i-residues+1,2,3)=coords((res_coor*(i))-3)
         !O
         peptides(i-residues+1,1,1)=coords((res_coor*(i))-2)
         peptides(i-residues+1,1,2)=coords((res_coor*(i))-1)
         peptides(i-residues+1,1,3)=coords((res_coor*(i)))
      else
         !C
         peptides(i+1,2,1)=coords((res_coor*(i))-5)
         peptides(i+1,2,2)=coords((res_coor*(i))-4)
         peptides(i+1,2,3)=coords((res_coor*(i))-3)
         !O
         peptides(i+1,1,1)=coords((res_coor*(i))-2)
         peptides(i+1,1,2)=coords((res_coor*(i))-1)
         peptides(i+1,1,3)=coords((res_coor*(i)))
      endif
   enddo

   call count_hb()
   write(30,*)dble(max_hb+num_par-num_anti)/dble(2*max_hb)
   write(*,*) "Stucture:",recno,"H Bonds:",num_hb,"Parallel:",num_par,"Antiparallel:",num_anti
!Reset counts
   num_hb=0
   num_par=0
   num_anti=0
   goto 150
200 close(10)
   close(30)
   recno=recno-1
   stop
end
