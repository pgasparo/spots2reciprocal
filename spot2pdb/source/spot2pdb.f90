! Kay Diederichs 5.3.2018 for use with coot, or with coot + dials.rs_mapper
!
! read SPOT.XDS and XPARM.XDS; write pseudo-PDB files visualizing reciprocal space
!
! "coot"ing the pseudo-PDB files needs the file ~/.coot with this line:
! (allow-duplicate-sequence-numbers)
!
! TODO: write out the intensity value in the PDB file's B-factor field.
! TODO: understand why -180 offset in phi for compatibility with dials.rs_mapper
!
! known limitation: detector segments are not implemented
! undesired feature: coot connects spots by "bonds" if they a) are close b) are in same chain c) have residue number difference <2
!                    the program tries to avoid this by cycling through chain ids, but this does not work perfectly.
!
! on Mac:
! gfortran -C -O -static-libgfortran -static-libgcc -Wno-argument-mismatch spot2pdb.f90 -o spot2pdb
! on Linux -static-libgcc is not needed
! 2020-02-21: remove existing SPOT-indexed.pdb if it exists

      IMPLICIT NONE
      INTEGER :: ih,ik,il,starting_frame,spacegroup,nseg,nx,ny,i,j,k 
      CHARACTER :: option*80,line*80,atomname*4,residuename*3,chain*1,spots*80='SPOT.XDS',xparm*80='XPARM.XDS'
      REAL :: rot_ax(3),starting_angle,oscillation_range,wavelength, &
           incident_beam(3),unit_cell(6),unitcell_axes(3,3),qx,qy,orgx,orgy,f, &
           det_x(3),det_y(3),det_z(3),x,y,phi,intensity,rot_r(3),r(3), &
           det,rr(3),resolmin=999,resolmax=6,xxx  ! resolmax=6A is dials.rs_mapper default
      LOGICAL :: idxref=.TRUE.
      
      PRINT '(a)','SPOT2PDB to create pseudo-PDB files from SPOT.XDS . KD 21/6/2019'
      
! read command line
      i=COMMAND_ARGUMENT_COUNT()
      IF (MOD(i,2)/=0) STOP 'number of arguments must be even; default: -q 999. -r 6. -s SPOT.XDS -x XPARM.XDS'
      DO j=1,i,2
        CALL GET_COMMAND_ARGUMENT(j,option)
        CALL GET_COMMAND_ARGUMENT(j+1,line)
        IF      (option(1:2)=='-q') THEN  
          READ(line,*) resolmin 
        ELSE IF (option(1:2)=='-r') THEN  
          READ(line,*) resolmax
        ELSE IF (option(1:2)=='-s') THEN  
          spots=TRIM(line) 
        ELSE IF (option(1:2)=='-x') THEN  
          xparm=TRIM(line) 
        ELSE
          WRITE(*,'(a,a)') 'unknown option: ',TRIM(option)
          STOP 'error reading command line'
        END IF
      END DO  
      PRINT '(a,f8.3)','min resolution (change with -q commandline option):  ',resolmin
      PRINT '(a,f8.3)','max resolution (change with -r commandline option):  ',resolmax
      PRINT '(a,a)'   ,'    spots file (change with -s commandline option):  ',TRIM(spots)
      PRINT '(a,a)'   ,'    xparm file (change with -x commandline option): ' ,TRIM(xparm)
      
      PRINT *,''
      
! read XPARM.XDS
      OPEN(1,FILE=xparm,ACTION='READ',IOSTAT=i)
      IF (i==0) THEN
        READ(1,'(a)') line 
        IF (INDEX(line,'XPARM.XDS')==0) STOP '1st line of xparm file is not XPARM.XDS'
        READ(1,*) starting_frame,starting_angle,oscillation_range,rot_ax
        READ(1,*) wavelength,incident_beam
        READ(1,*) spacegroup,unit_cell
        READ(1,*) (unitcell_axes(1,i),i=1,3) 
        READ(1,*) (unitcell_axes(2,i),i=1,3) 
        READ(1,*) (unitcell_axes(3,i),i=1,3) 
        READ(1,*) nseg,nx,ny,qx,qy   
        IF (nseg>1) STOP 'number of detector segments > 1 not implemented'
        READ(1,*) orgx,orgy,f
        READ(1,*) det_x                 ! detector x-axis
        READ(1,*) det_y                 ! detector y-axis
! no further lines from XPARM.XDS are read
      ELSE
        PRINT*,'reading geometry from XDS.INP'
        OPEN(1,FILE='XDS.INP',ACTION='READ',IOSTAT=i)
        IF (i/=0) STOP 'could open neither XPARM.XDS nor XDS.INP'
        CALL getfromxdsinp('STARTING_ANGLE',1,starting_angle)
        IF (starting_angle==HUGE(starting_angle)) starting_angle=0.
        CALL getfromxdsinp('OSCILLATION_RANGE',1,oscillation_range)
        CALL getfromxdsinp('ROTATION_AXIS',3,rot_ax)  
        CALL getfromxdsinp('X-RAY_WAVELENGTH',1,wavelength)
        CALL getfromxdsinp('INCIDENT_BEAM_DIRECTION',3,incident_beam)
        CALL getfromxdsinp('NX',1,xxx) ; nx=xxx
        CALL getfromxdsinp('NY',1,xxx) ; ny=xxx
        CALL getfromxdsinp('QX',1,qx) 
        CALL getfromxdsinp('QY',1,qy) 
        CALL getfromxdsinp('ORGX',1,orgx) 
        CALL getfromxdsinp('ORGY',1,orgy) 
        CALL getfromxdsinp('DETECTOR_DISTANCE',1,f) 
        CALL getfromxdsinp('DIRECTION_OF_DETECTOR_X-AXIS',3,det_x) 
        CALL getfromxdsinp('DIRECTION_OF_DETECTOR_Y-AXIS',3,det_y)    
      END IF
      CLOSE(1)
      det_z(1)=det_x(2)*det_y(3)-det_x(3)*det_y(2)   ! calculate detector normal -
      det_z(2)=det_x(3)*det_y(1)-det_x(1)*det_y(3)   ! XDS.INP does not have
      det_z(3)=det_x(1)*det_y(2)-det_x(2)*det_y(1)   ! this item.
      det_z = det_z/SQRT(DOT_PRODUCT(det_z,det_z))   ! normalize (usually not req'd)
      
      det=200./resolmax   ! unit cell axes; typically 33.333 if resolmax=6
! read SPOT.XDS, and write SPOT-indexed.pdb and SPOT-notindexed.pdb
      OPEN(1,FILE=spots,ACTION='READ',IOSTAT=i)
      IF (i/=0) STOP 'could not open spots file'
      
! find out whether SPOT.XDS is from COLSPOT or IDXREF
      READ(1,*,IOSTAT=i) x,y,phi,intensity,ih,ik,il
      REWIND(1)
      OPEN(2,FILE='SPOT-indexed.pdb',ACTION='WRITE',STATUS='REPLACE') ! remove old file
      IF (i/=0) THEN         ! is from COLSPOT
        idxref=.FALSE.       ! lines don't have ih,ik,il
        ih=0
        ik=0
        il=0
      ELSE                   ! is from IDXREF
        WRITE(2,'(a,3f9.3,a)') 'CRYST1',det,det,det,'  90.00  90.00  90.00 P 1'
      END IF
      OPEN(3,FILE='SPOT-notindexed.pdb',ACTION='WRITE')
      WRITE(3,'(a,3f9.3,a)') 'CRYST1',det,det,det,'  90.00  90.00  90.00 P 1'
      i=0   ! reflection counter
! sequence number for not-indexed reflections; increment is 2 to prevent bonds
      j=0  
      k=0   ! for cycling through dummy chains, using ASCII characters from " " to "~"
      DO
        IF (idxref) THEN       ! SPOT.XDS re-written by IDXREF
          READ(1,*,end=99) x,y,phi,intensity,ih,ik,il 
        ELSE                   ! SPOT.XDS written by COLSPOT
          READ(1,*,end=99) x,y,phi,intensity
        END IF

! convert detector coordinates to local coordinate system
        r(1)=(x-orgx)*qx*det_x(1) + (y-orgy)*qy*det_y(1) +f*det_z(1)
        r(2)=(x-orgx)*qx*det_x(2) + (y-orgy)*qy*det_y(2) +f*det_z(2)
        r(3)=(x-orgx)*qx*det_x(3) + (y-orgy)*qy*det_y(3) +f*det_z(3)

! normalize scattered vector to obtain S1
        r=r/(wavelength*SQRT(DOT_PRODUCT(r,r)))

! obtain reciprocal space vector S = S1-S0
        r=r-incident_beam
        IF (SQRT(DOT_PRODUCT(r,r))>1./resolmax) CYCLE  ! outer resolution limit
        IF (SQRT(DOT_PRODUCT(r,r))<1./resolmin) CYCLE  ! inner resolution limit

! rotate  
! NB: the term "-180." (found by trial&error) seems to make it match dials.rs_mapper
        phi=(starting_angle+oscillation_range*phi -180.)/180.*3.141592654   
        CALL rodrigues(r,phi,rot_ax,rot_r)
        rot_r=100.*rot_r + 100./resolmax  ! transform to match dials.rs_mapper
        i=i+1
        k=k+1
        WRITE(chain,'(a)') char(k+31)  ! ASCII code 32 is space (blank); alloweded as chain id
        IF (k+31==126) k=0             ! ASCII code 126 is ~    (tilde); also allowed         
        IF (ih==0.AND.ik==0.AND.il==0) THEN   ! not indexed
          j=j+2               ! 10000 would overflow the sequence number field
          IF (j>9999) j=-999  ! NB we may have repeated sequence numbers
          atomname=  '   X' 
          residuename='  Y'
          WRITE(3,101) i,atomname,' ',residuename,chain,j,rot_r,1.,10.,'O' ! red atom
        ELSE 
          WRITE(atomname,'(i4)')ih
          WRITE(residuename,'(i3)')ik
          WRITE(2,101) i,atomname,' ',residuename,chain,il,rot_r,1.,10.,'C' ! yellow atom
        END IF
      END DO
      
  99  CONTINUE
 
      PRINT '(i0,a)',i,' reflections from min to max resolution written to:'
      IF (idxref) THEN
        PRINT '(a)','SPOT-indexed.pdb with (atomname residuename sequence#) = (H K L)'
        PRINT '(a)','SPOT-notindexed.pdb with dummy atomname residuename sequence#'
        PRINT *,''
        PRINT '(a)','If ~/.coot has the line: (allow-duplicate-sequence-numbers) , you can:'
        PRINT '(a)','coot SPOT-*.pdb'
        PRINT '(a)','clicking on atoms of SPOT-indexed.pdb shows their H K L (ignore chain)'
        PRINT *,''
        PRINT '(a)','After e.g. "dials.rs_mapper /raw/data/of/SPOT_RANGE/*.cbf"  , you can: '
        PRINT '(a)','coot SPOT-*.pdb --map rs_mapper_output.ccp4'
! write coordinates for rotation axis and origin
        CALL axisandorigin(2,resolmax,rot_ax)
      ELSE
        PRINT '(a)','SPOT-notindexed.pdb with dummy atomname residuename sequence#'
        PRINT *,''
        PRINT '(a)','If ~/.coot has the line: (allow-duplicate-sequence-numbers) , you can:'
        PRINT '(a)','coot SPOT-notindexed.pdb'
        PRINT *,''
        PRINT '(a)','After e.g. "dials.rs_mapper /raw/data/of/SPOT_RANGE/*.cbf"  , you can: '
        PRINT '(a)','coot SPOT-notindexed.pdb ROTATION_AXIS.pdb --map rs_mapper_output.ccp4'
! write coordinates for rotation axis and origin
        CALL axisandorigin(3,resolmax,rot_ax)
      END IF
      PRINT '(a)','In coot, adjust the map level and map radius!'
 101  FORMAT('HETATM',i5,t13,a4,a1,a3,t22,a1,i4,t31,3f8.3,t55,f6.2,t61,f6.2,t78,a1)
!  http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#HETATM  
      END

!   
      SUBROUTINE rodrigues(h,phi,rot_ax,rot_h)
 ! this implements Rodrigues' rotation formula (see Wikipedia)
! {v}_rot = {v} cos(theta) + {k} x {v} sin(theta) + {k} {k}.{v}) (1 - cos(theta))
! below is from http://mathworld.wolfram.com/RodriguesRotationFormula.html
!
! -phi is used because we are rotating back to to unrotated crystal orientation
      IMPLICIT NONE
      REAL rot_h(3),rot_ax(3),h(3),phi,omcp,cp,sp
      cp=cos(-phi)
      sp=sin(-phi)
      omcp=1.-cp
      rot_h(1)= (                  cp+rot_ax(1)**2*omcp)*h(1)+ &
                (-rot_ax(3)*sp+rot_ax(1)*rot_ax(2)*omcp)*h(2)+ &
                ( rot_ax(2)*sp+rot_ax(1)*rot_ax(3)*omcp)*h(3)
      rot_h(2)= ( rot_ax(3)*sp+rot_ax(1)*rot_ax(2)*omcp)*h(1)+ &
                (                  cp+rot_ax(2)**2*omcp)*h(2)+ &
                (-rot_ax(1)*sp+rot_ax(2)*rot_ax(3)*omcp)*h(3)
      rot_h(3)= (-rot_ax(2)*sp+rot_ax(1)*rot_ax(3)*omcp)*h(1)+ &
                ( rot_ax(1)*sp+rot_ax(2)*rot_ax(3)*omcp)*h(2)+ &
                (                  cp+rot_ax(3)**2*omcp)*h(3)
      END

      SUBROUTINE axisandorigin(iunit,resolmax,rot_ax)
! place dummy N atoms on ROTATION_AXIS from 0 to 200/resolmax, and mark origin 
      IMPLICIT NONE
      REAL resolmax,rot_ax(3),ori
      INTEGER iunit,i
      
      ori=100./resolmax
      DO i=-INT(ori/1.5),INT(ori/1.5)
        WRITE(iunit,101) 0,'   N',' ','  X',' ',0,ori+rot_ax*i*1.5,1.,1.,'N' 
      END DO
! mark the origin with a 3D cross
      WRITE(iunit,101) 0,'   N',' ','  X',' ',0,ori+1.5,ori,ori,1.,1.,'N' ! blue
      WRITE(iunit,101) 0,'   N',' ','  X',' ',0,ori-1.5,ori,ori,1.,1.,'N'
      WRITE(iunit,101) 0,'   N',' ','  X',' ',0,ori,ori+1.5,ori,1.,1.,'N'
      WRITE(iunit,101) 0,'   N',' ','  X',' ',0,ori,ori-1.5,ori,1.,1.,'N'
      WRITE(iunit,101) 0,'   N',' ','  X',' ',0,ori,ori,ori+1.5,1.,1.,'N'
      WRITE(iunit,101) 0,'   N',' ','  X',' ',0,ori,ori,ori-1.5,1.,1.,'N'
 101  FORMAT('ATOM  ',i5,t13,a4,a1,a3,t22,a1,i4,t31,3f8.3,t55,f6.2,t61,f6.2,t78,a1)
      END
      
      SUBROUTINE getfromxdsinp(string,length,array)
      IMPLICIT NONE
      CHARACTER*(*), INTENT(IN) :: string
      INTEGER, INTENT(IN) :: length
      REAL, INTENT(OUT) :: array(length)
      
      CHARACTER*132 line
      INTEGER i

      DO 
        READ(1,'(a)',END=99,ERR=99) line
        IF (INDEX(line,'!')>0) line(INDEX(line,'!'):)=''
        i=INDEX(line,string)
        IF (i>0) THEN
          READ(line(i+LEN(string)+1:),*,ERR=99) array
          REWIND(1)
          RETURN
        END IF
      END DO
 99   CONTINUE
      IF (string=='STARTING_ANGLE') THEN ! special case - not specified in XDS.INP
        array(1)=HUGE(array) 
      ELSE
        PRINT*,'could not find ',TRIM(string),' or read its value'
        STOP
      END IF
      REWIND(1)
      END SUBROUTINE getfromxdsinp

