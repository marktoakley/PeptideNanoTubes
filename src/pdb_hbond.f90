program hbond
use hb
! Analyse hydrogen-bonding pattern in cyclic peptide nanotubes
! Get input from a pdb file
implicit none
character*50 arg,filename
character*4 atname
integer resnum
double precision x, y, z
integer :: ires=0, iring=0
integer i,j,k
logical :: finished=.false.

if (COMMAND_ARGUMENT_COUNT().lt.3) then
   write(*,*) "usage: TubeHbond <string file>  <int residues> <int rings>"
   stop
endif


call GET_COMMAND_ARGUMENT(1,arg)
read(arg,*)filename
call GET_COMMAND_ARGUMENT(2,arg)
read(arg,*)residues
call GET_COMMAND_ARGUMENT(3,arg)
read(arg,*)rings

! Atom indices in peptides array
! 1 O, 2 C, 3 N, 4 H

open(10,file=filename)
do while (.not.finished)
   read(10,'(12X,A4,6X,I4,4X,3(F8.3))')atname,resnum,x,y,z
   if (atname=="N   ") then
!      write(*,*) atname ,resnum,x,y,z
      peptides(resnum,3,1)=x
      peptides(resnum,3,2)=y
      peptides(resnum,3,3)=z
   else if (atname=="H   ") then
!      write(*,*) atname ,resnum,x,y,z
      peptides(resnum,4,1)=x
      peptides(resnum,4,2)=y
      peptides(resnum,4,3)=z
   else if (atname=="C   ") then
!      write(*,*) atname ,resnum,x,y,z
      if (mod(resnum,residues).eq.0) then
         peptides(resnum-residues+1,2,1)=x
         peptides(resnum-residues+1,2,2)=y
         peptides(resnum-residues+1,2,3)=z
      else
         peptides(resnum+1,2,1)=x
         peptides(resnum+1,2,2)=y
         peptides(resnum+1,2,3)=z
      endif
   else if (atname=="O   ") then
!      write(*,*) atname ,resnum,x,y,z
      if (mod(resnum,residues).eq.0) then
         peptides(resnum-residues+1,1,1)=x
         peptides(resnum-residues+1,1,2)=y
         peptides(resnum-residues+1,1,3)=z
      else
         peptides(resnum+1,1,1)=x
         peptides(resnum+1,1,2)=y
         peptides(resnum+1,1,3)=z
      endif
      if (resnum.eq.rings*residues) then
         finished=.true.
   endif
 endif
enddo

call count_hb()

write(*,*)"HBonds: max",max_hb,"total",num_hb,"parallel",num_par,"antiparallel",num_anti

stop
end
