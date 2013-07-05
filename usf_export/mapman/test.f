      program mapman

      include 'mapman.incl'
c
      integer maxpnt,maxmap
      parameter (maxmap = maxgk4)
      parameter (maxpnt = 100000)
      real map(maxpnt,maxmap),buffer(maxpnt)
c
      real xdum,ave,sdv
c
      integer imap,ifmt,ierr,i,jfmt,level,iunit,imax,jmax,kmax, IOS
c
code ...
c
      ifmt = 4
      ierr = -1
      jfmt = -ifmt
      level = 0
      iunit = 1
c      call prompt (' Read header')
      incore (imap) = .false.

 OPEN(UNIT=iunit, FILE='mapman_test.xE', STATUS='UNKNOWN', ACCESS='SEQUENTIAL'
     +           ,FORM='UNFORMATTED', ERR=5, IOSTAT=IOS)

     print *, 'MRH SUCCESS ', IOS
      return
 5    print *, 'MRH ERROR ', IOS
      return

      call maphdr_test ('mapman_test.xE', iunit, jfmt,
     +  origin(1,imap), extent(1,imap), grid(1,imap),
     +  uvw(1,imap), cell(1,imap), spaceg(imap),
     +  rho, maxrho, map(1,imap), maxpnt, incore(imap))
c

      if (spaceg(imap).le.0) then
        call errcon ('While opening map file')
        ierr = -2
        return
      end if
c
      call prompt (' Header done')
      return
      end

C============================================================================

      subroutine maphdr_test (file, in, flag, 
     $                   ioxyz, imxyz, mygrid, iuvw, cell, spgrp,
     $                   rho, size, savrho, sizsav, incore)
c
c ---	A routine to read in 'header' information from the
c	density maps
c ---	For CCP4 maps, if order is wrong, it gets stored away in-core
c ---	Alwyn Jones 3-Aug-91
c
c ---	Alterations
c
c ---	Phil Evans   MRC LMB, Cambridge                August 1991
c	Added option to read formatted CCP4 maps and 
c	added common block /cmscal/
c
c ---   Gerard Kleywegt March/April 1993
c       - added MASK format
c       - improved error trapping
c       - added NEWEZD format
c       - added binary XPLOR format
c
c --- gjk @ 940212 - interface with new MAPCLO routine
c                    to properly close read maps
c
      implicit none
c
      character file*(*)
      integer size, sizsav
      integer nfast,nmed,nslow
      integer cmnd, imxyz(3), in, ioxyz(3), iuvw(3), mygrid(3), 
     $  spgrp, flag
      real cell(6), rho(size), savrho(sizsav)
      real*8 drcell(6)
      logical incore,xinter
c
c --- For CCP4 maps
c
      integer nu, nv
      logical lform
      real cscale,cplus
      common /cmscal/ cscale,cplus,nu,nv,lform
      save /cmscal/
c
c ---	For Protein maps
c
      integer title(20),nrho
      character*4 title_c(20) ! mrh   

      real rhorng(6)
      integer ibp0, iperm(3), igrid(3), mifd(9)
c
      character tit*80, fmt*80, line*80
      character qxyz(3)
      integer titnum
      integer lmode, i, ierror, ipt, ip1, ip2, ir1, ir2, isec, 
     $  j, nsym, numbay
      integer ct, ct1,ierr
      real scale, avea,sdva,vara,xmin,xmax,xtot
c
      integer*2 jj1,jj2,jj3,jj4,jj5
c
ccc      data  title /20*'    '/
      data  title_c /20*'    '/ ! mrh   
      equivalence (title, title_c) ! mrh

      data qxyz /'X','Y','Z'/
c
code ...
c
      call textut (' Input map :',file)
c
      cmnd = iabs(flag)
      incore = .false.
      spgrp = 1
c
      call setmgt (in,cmnd)
c
c --- "PROTEIN" style map
c
      if (cmnd .eq. 1) then
c ---	Rewind the electron density file.
        close (in)
c
        call xopxob (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
        rhorng(1) = 10000000.0
        rhorng(2) = -10000000.0
        rhorng(3) = 0.0
        rhorng(4) = 0.0
        rhorng(5) = 0.0
        rhorng(6) = 0.0
        nrho	= 0
        ibp0	= 0
c ---	  Read the type 0 header.
        call getr0 (in, title, ierror)
        if (ierror .ne. 0) goto 8900
        write(6,10) title
10      format (/,' Map title is :',20A4,/)
c ---	  Read the type 1 header and symmetry card.
        call getr1 (in, cell, nsym, ierror)
        if (ierror .ne. 0) goto 8900
c ---	  Read the type 2 map index parameter record
        call getr2 (in, iperm, igrid, mifd, ierror)
        if (ierror .ne. 0) goto 8900
c ---	  Load common blocks
        do 100 i=1,3
          ipt = (i-1)*3+1
          ioxyz(iperm(i)) = mifd(ipt)
          imxyz(iperm(i)) = mifd(ipt+1)- mifd(ipt)+1
100       mygrid(iperm(i)) = igrid(i)
        iuvw(1)=iperm(2)
        iuvw(2)=iperm(1)
        iuvw(3)=iperm(3)
c
	else if (cmnd .eq. 2 .or. cmnd .eq. 3) then
c
c ---   Ten-eycke style maps: cmnd = 2 for FFT-Y, 3 for TENEYCK2
c
        close(in)
c
        call xopxob (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
c ---	  read title and as much information as possible
        iuvw(1) = 2
        iuvw(2) = 1
        iuvw(3) = 3
        call ivalin (' Fast, medium, slow axes ?',3,iuvw)
cc        write (6,20)
cc        read (5,*) iuvw
        read (in,err=8900) title
        if (cmnd .eq. 2) then
          read (in,err=8900) isec, ip1, ip2, ir1, ir2	! SECN, FAST, MEDIUM AXES
        else					!TENEYCK2
          read (in,err=8900) jj1, jj2, jj3, jj4, jj5
          isec = jj1
          ir1 = jj2
          ir2 = jj3
          ip1 = jj4
          ip2 = jj5
          rewind (in)
c ---       Skip over title
          read (in,err=8900)
        end if
        mygrid (1) = 10
        mygrid (2) = 10
        mygrid (3) = 10
        call ivalin (' Grid ?',3,mygrid)
cc        write (6,30)
cc        read (5,*) mygrid
        numbay = 10
        call ivalin (' Number of Y-sections ?',1,numbay)
cc	  write (6,40)
cc        read (5,*) numbay
        cell (1) = 100.0
        cell (2) = 100.0
        cell (3) = 100.0
        cell (4) =  90.0
        cell (5) =  90.0
        cell (6) =  90.0
        call fvalin (' Unit cell constants ?',6,cell)
cc        write (6,80)
cc        read (5, *) cell
c
        ioxyz(iuvw(1))=ip1
        ioxyz(iuvw(2))=ir1
        ioxyz(iuvw(3))=isec
        imxyz(iuvw(1))=ip2-ip1+1
        imxyz(iuvw(2))=ir2-ir1+1
        imxyz(iuvw(3))=numbay
c
      else if (cmnd .eq. 4) then
c
c ---    CCP4 maps, formatted or binary
c
c Is it formatted or binary? lform = .true. if formatted
c
       call prompt ('MGULPR approach 1')

        call opnmfl_test(file,lform,in,ierr)

       call prompt ('MGULPR leave')

        if (ierr .ne. 0) goto 9000
c
        if (lform) then
c Formatted, so read formatted header
          call rdfhdr (in, file, tit, numbay, iuvw, mygrid,
     $        isec,ip1,ip2,ir1,ir2,cell,spgrp,lmode,
     $        rhorng(1),rhorng(2),rhorng(3),rhorng(4),
     $        cscale, cplus)
c
c ... flag file type as FORMATTED CCP4
c
          call setmgt (in,-4)
c
        else
c Binary
          call mrdhdr (in, file, tit, numbay, iuvw, mygrid,
     $        isec,ip1,ip2,ir1,ir2,cell,spgrp,lmode,
     $        rhorng(1),rhorng(2),rhorng(3),rhorng(4))
        endif
cxyz
c        call prompt ('MRDHDR done')
c
c ---	 Switch these arrays around, so ordered according to x,y,z.
c
        ioxyz(iuvw(1)) = ip1
        ioxyz(iuvw(2)) = ir1
        ioxyz(iuvw(3)) = isec
        imxyz(iuvw(1)) = ip2-ip1+1
        imxyz(iuvw(2)) = ir2-ir1+1
        imxyz(iuvw(3)) = numbay
        nu = ip2-ip1+1
        nv = ir2-ir1+1
c
        if (iuvw(1) .eq. 2 .and. iuvw(2) .eq. 1 .and. iuvw(3) .eq. 3)
     $    goto 130
        if (imxyz(1)*imxyz(2)*imxyz(3) .gt. sizsav) goto 130
        do 110 i=1,numbay
          if (lform) then
            call rdfsec(in, ierror, rho, nu, nv, cscale, cplus)
          else


c ... 940225 - replace MGULP by MGULPR so it can read CCP4 masks
c
            call mgulpr (in, rho, ierror) ! CCP4 routine
cxyz
c        call prompt ('MGULPR done')
c
            if (ierror .ne. 0) goto 8900
          endif
cxyz
c        call prompt ('PCKRHO call next')
          call pckrho (savrho, imxyz(1), imxyz(2), imxyz(3), i,
     $          rho, imxyz(iuvw(1)), imxyz(iuvw(2)), iuvw)
cxyz
c        call prompt ('PCKRHO done')
c
110     continue
        incore = .true.
c
c ---	X-plor style giant maps. Assume sorted ZYX
c     010513 - also AMBER-style xplor maps (type=23)
c
      else if (cmnd .eq. 5 .or. cmnd .eq. 23) then
c
        call xopxoa (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
        read (in, 1020,err=8900) i
        do 140 j=1,i
          read (in,1030,err=8900) tit
          call textut (' Title :',tit)
140     continue
        read (in, 1050,err=8900) (mygrid(i), ioxyz(i), imxyz(i), i=1,3)
        do 150 i=1,3
150       imxyz(i) = imxyz(i)- ioxyz(i) + 1
        read (in, 1060,err=8900) cell
        read (in, 1030,err=8900) tit
        if (tit(1:3) .ne. 'ZYX') then
          call errcon (' Map should be sorted ZYX !')
          goto 9000
        end if
        iuvw(1) = 2
        iuvw(2) = 1
        iuvw(3) = 3
c
c ---	 OLDEZDensity style giant maps
c
      else if (cmnd .eq. 6) then
c
        call xopxoa (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
        read (in, 2010,err=8900) ioxyz
        read (in, 2010,err=8900) imxyz
        read (in, 2010,err=8900) mygrid
        read (in, 2020,err=8900) cell
        read (in, 2030,err=8900) line
        if (imxyz(1)*imxyz(2)*imxyz(3) .gt. sizsav) then
          call errcon (' Map too big !')
          call jvalut (' Requested size :',1,
     +                 (imxyz(1)*imxyz(2)*imxyz(3)))
          call jvalut (' Available size :',1,sizsav)
          goto 9000
        end if
c
        if (line(1:4) .ne. 'MAP ') then
          call errcon ('Fifth line of OLDEZD header corrupted')
          goto 9000
        end if
        ct = 4
2100    if (line(ct:ct) .eq. ' ') then
          ct = ct+1
          if (ct .gt. 80) then
            call errcon ('Fifth line of OLDEZD header corrupted')
            goto 9000
          end if
          goto 2100
        end if
        ct1 = ct
2110    if (line(ct:ct) .ne. ')') then
          ct = ct+1
          if (ct .gt. 80) then
            call errcon ('Fifth line of OLDEZD header corrupted')
            goto 9000
          end if
          goto 2110
        end if
        fmt = line(ct1:ct)
        line(1:ct) = ' '
        call remspa (line)
cc        call alignl (line)
        read (line, 2060,err=8900) scale
        call rvalut (' Scale constant :',1,scale)
cc        type*,'scale constant=', scale
        ct = imxyz(1)*imxyz(2)*imxyz(3)
        read (in, fmt, end=8900, err=8900) (savrho(i), i=1,ct)
        do 2120 i=1, ct
2120      savrho(i) = savrho(i)/scale
        iuvw(1) = 2
        iuvw(2) = 1
        iuvw(3) = 3
        incore = .true.
        call xstats (savrho,ct,avea,sdva,xmin,xmax,xtot)
        vara = sdva*sdva
        write (*,6010) ct,xmin,xmax,xtot,avea,vara,sdva
c
        goto 130
c
c ---	 MASKs
c
      else if (cmnd .eq. 7) then
c
        call xopxoa (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
        call maskin (in,savrho,ioxyz,imxyz,mygrid,cell,sizsav,ierr)
        if (ierr .ne. 0) goto 9000
        ct = imxyz(1)*imxyz(2)*imxyz(3)
c
c ... MASKIN has written integers into array SAVRHO
c     convert them back into "real integers"
c
        do i=1,ct
          call r2r (savrho(i),savrho(i))
        end do
        iuvw(1) = 2
        iuvw(2) = 1
        iuvw(3) = 3
        incore = .true.
        goto 130
c
c ---	 (NEW)EZDensity style giant maps
c
      else if (cmnd .eq. 8) then
c
        call xopxoa (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
        do i=1,3
          ioxyz(i) = 0
          imxyz(i) = 0
          mygrid(i) = 0
          cell(i) = 0.0
          cell(i+3) = 0.0
        end do
        scale = 0.0
c
        read (in,'(a)',err=8900) line
        call upcase (line)
        if (line (1:7) .ne. 'EZD_MAP') then
          call errcon ('Not a (NEW)EZD file !')
          goto 8900
        end if
c
 2873   continue
        read (in,'(a)',err=8900) line
c
        if (line(1:1) .eq. '!') then
          call textut (' >',line)
          goto 2873
        end if
c
        call upcase (line)
        if (line(1:4).eq.'CELL') then
          read (line(5:),*,err=8900) cell
          goto 2873
        end if
c
        if (line(1:6).eq.'ORIGIN') then
          read (line(7:),*,err=8900) ioxyz
          goto 2873
        end if
c
        if (line(1:6).eq.'EXTENT') then
          read (line(7:),*,err=8900) imxyz
          goto 2873
        end if
c
        if (line(1:4).eq.'GRID') then
          read (line(5:),*,err=8900) mygrid
          goto 2873
        end if
c
        if (line(1:5).eq.'SCALE') then
          read (line(6:),*,err=8900) scale
          call rvalut (' Scale constant :',1,scale)
          goto 2873
        end if
c
        if (line(1:3) .eq. 'MAP') goto 2876
c
        call errcon ('Invalid line in (NEW)EZD file')
        call textut (' >',line)
        goto 2873
c
 2876   continue
        if (scale.eq.0.0 .or.
     +      imxyz(3).eq.0 .or. cell(6).eq.0.0) then
          call errcon ('Not all data in header')
          goto 8900
        end if
c
        if (imxyz(1)*imxyz(2)*imxyz(3) .gt. sizsav) then
          call errcon (' Map too big !')
          call jvalut (' Requested size :',1,
     +                 (imxyz(1)*imxyz(2)*imxyz(3)))
          call jvalut (' Available size :',1,sizsav)
          goto 9000
        end if
c
        ct = imxyz(1)*imxyz(2)*imxyz(3)
        read (in, *, end=8900, err=8900) (savrho(i), i=1,ct)
        do 2920 i=1, ct
2920      savrho(i) = savrho(i)/scale
        iuvw(1) = 2
        iuvw(2) = 1
        iuvw(3) = 3
        incore = .true.
        call xstats (savrho,ct,avea,sdva,xmin,xmax,xtot)
        vara = sdva*sdva
        write (*,6010) ct,xmin,xmax,xtot,avea,vara,sdva
c
        goto 130
c
c ---	BINARY X-plor style giant maps. Assume sorted ZYX
c
      else if (cmnd .eq. 9) then
c
        call xopxob (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
        read (in,err=8900) i,(tit(1:80),j=1,i)
        call textut (' Title :',tit)
c
        read (in,err=8900)
     +    mygrid(1),ioxyz(1),imxyz(1),
     +    mygrid(2),ioxyz(2),imxyz(2),
     +    mygrid(3),ioxyz(3),imxyz(3)
cc      print *,(mygrid(i), ioxyz(i), imxyz(i), i=1,3)
        do 139 i=1,3
139       imxyz(i) = imxyz(i) - ioxyz(i) + 1
cc        read (in) i
cc      print *,' I = ',i
        read (in,err=8900) drcell(1),drcell(2),drcell(3),
     +    drcell(4),drcell(5),drcell(6)
cc      print *,'... CELL ...',(drcell(i),i=1,6)
        do i=1,6
          cell(i)=drcell(i)
        end do
        read (in,err=8900) (tit(i:i),i=1,3)
cc      print *,tit(1:3)
        if (tit(1:3) .ne. 'ZYX') then
          call errcon (' Map should be sorted ZYX !')
          goto 9000
        end if
        iuvw(1) = 2
        iuvw(2) = 1
        iuvw(3) = 3
c
c ... check if not created on a 64-bit machine
c
        do i=1,3
          if (mygrid(i).le.1 .or. imxyz(i).le.1 .or.
     +        cell(i).le.1.0 .or. cell(i+3).le.1.0) then
            write (6,70) ioxyz,imxyz,mygrid,
     +        cell,(qxyz(iuvw(j)),j=1,3)
            call errcon ('Map probably created on a 64-bit machine')
            call errcon ('Cannot handle this !')
            goto 9000
          end if
        end do
c
        else if (cmnd .eq. 13) then
c
c ---   TNT style maps: cmnd = 13
c
        close(in)
c
        call xopxob (in,file,xinter(),ierr)
        if (ierr .ne. 0) goto 8800
c
c ---     read title and as much information as possible
c         Title records start with a '*'
c         except the last one...
c      
        titnum = 0
141     read(in) tit
        titnum = titnum +1
        if (tit(1:1).eq.'*') goto 141
        read(in) isec, ir1, ir2, ip1, ip2, nfast, nmed, nslow
c
c  Now we reposition the map file at the first record that contains
c  electron density information
c
        rewind(in)
        do 142 i=1,titnum
           read(in) tit
           call textut (' Title :',tit)
142     continue
        iuvw(1) = 1
        iuvw(2) = 2
        iuvw(3) = 3
        call ivalin (' Fast, medium, slow axes ?',3,iuvw)
        numbay = 10
        call ivalin (' Number of sections ?',1,numbay)
        cell (1) = 100.0
        cell (2) = 100.0
        cell (3) = 100.0
        cell (4) =  90.0
        cell (5) =  90.0
        cell (6) =  90.0
        call fvalin (' Unit cell constants ?',6,cell)
c
        mygrid(1) = nfast
        mygrid(2) = nmed
        mygrid(3) = nslow
        ioxyz(iuvw(1))=ir1
        ioxyz(iuvw(2))=ip1
        ioxyz(iuvw(3))=isec
        imxyz(iuvw(1))=ir2-ir1+1
        imxyz(iuvw(2))=ip2-ip1+1
        imxyz(iuvw(3))=numbay
c
c --- Command not found...
c
      else
        call errcon ('Unknown MAP format')
        call ivalut (' Format :',1,cmnd)
cc        write(6,60) cmnd
        goto 9000
      endif
c
c --- Come here to output some facts about map...
c
 130  continue
      write (6,70)ioxyz,imxyz,mygrid,cell,(qxyz(iuvw(i)),i=1,3)
c
      if (imxyz(iuvw(1))*imxyz(iuvw(2)) .gt. size) then
        call errcon (' Levels contain too many points')
        call jvalut (' Available :',1,size)
        call jvalut (' Requested :',1,
     +               (imxyz(iuvw(1))*imxyz(iuvw(2))))
        goto 9000
      end if
c
      if (imxyz(1)*imxyz(2)*imxyz(3) .gt. sizsav) then
        call errcon (' Map too big !')
        call jvalut (' Requested size :',1,
     +               (imxyz(1)*imxyz(2)*imxyz(3)))
        call jvalut (' Available size :',1,sizsav)
        goto 9000
      end if
c
      return
c
c ... error traps
c
 8800 call errcon ('While opening input file')
      call mapclo (in)
      goto 9000
c
 8900 call errcon ('While reading input file')
      call mapclo (in)
      goto 9000
c
 9000 continue
      if (flag .gt. 0) then
        call errstp ('MAPHDR - Sorry !')
      end if
      call errcon ('While reading map header. Sorry !')
      spgrp = -1
c
      return
c
70    format(/,' Parameters as read from the map file',/,
     $ ' Origin ...................... ',3I10,/,
     $ ' Extent ...................... ',3I10,/,
     $ ' Grid ........................ ',3I10,/,
     + ' Cell axes ................... ',3f10.2/
     + ' Cell angles ................. ',3f10.2/
     $ ' UVW (fast, medium, slow) .... ',3a10,/)
1020  format (/,i8)
1030  format (a)
1040  format (1x,a)
1050  format (9i8)
1060  format (6e12.5)
2010  format(3i5)
2020  format(6f10.0)
2030  format (a)
2040  format (' No space to store EZD map')
2050  format (' Fifth line wrong')
2060  format (e20.6)
2070  format (' Error reading map file')
c
 6010 format (' Size ',i10/
     +        ' ED min, max, total ',1p,3e12.4/
     +        ' ED ave, var, stdev ',3e12.4/)
c
      end
c
      subroutine opnmfl_test(fname,lform,iunit,ierr)
C     ====================================
C Find out if input file FNAME is formatted or not
C Leave file closed
c
      implicit none
c
      character*(*) fname
C LFORM is true if input file in formatted
      logical lform
c
      integer index, kdummy, kfail, iunit, ierr
      character*80 line
c
code ...
c
      lform = .false.
      kdummy = 0
      kfail = 1
      print *,'Approaching ccpdpn fname : ',iunit, fname

      call ccpdpn_test(iunit,fname,'READONLY','F',kdummy,kfail)
      if (kfail .eq. -1) then
        ierr = -1
      print *,'failing ccpdpn'
        return
      end if

c
C Read 1st line as is formatted, should contain string MAPEXCHANGE if so
      read (iunit,6001,err=100) line
 6001 format(A)
c
      if (index(line,'MAPEXCHANGE') .gt. 0) then
         lform = .true.
      else
         lform = .false.
      endif
c
C Close file
 100  close (unit=iunit)
      ierr = 0
c
      return 
      end
c
c
      SUBROUTINE CCPDPN_test(IUN,LOGNAM,STATUS,TYPE,LREC,IFAIL)
C     ====================================================
C
C---- Calls CCPOPN to open a file, but with mnemonic arguments
C
C Arguments:
C ==========
C
C         IUN (I)   INTEGER: UNIT NUMBER
C
C      LOGNAM (I)   CHARACTER*(*): LOGICAL FILE NAME
C
C      STATUS (I)   CHARACTER*(*): FILE STATUS FLAG:
C                                     'UNKNOWN'
C                                     'SCRATCH'
C                                     'OLD'
C                                     'NEW'
C                                     'READONLY'
C                                     'PRINTER'
C
C        TYPE (I)   CHARACTER*(*): FILE TYPE FLAG:
C                                  ='F', 'SEQUENTIAL' 'FORMATTED'
C                                  ='U', 'SEQUENTIAL' 'UNFORMATTED'
C                                  ='DF', 'DIRECT'     'FORMATTED'
C                                  ='DU', 'DIRECT'     'UNFORMATTED'
C     [STATUS and TYPE are case-insensitive]
C
C        LREC (I)   INTEGER: RECORD LENGTH FOR DIRECT ACCESS FILE (NO. OF
C                   CHARACTERS FOR A FORMATTED FILE OR WORDS FOR
C                   AN UNFORMATTED FILE). NOT RELEVANT FOR A SEQUENTIAL
C                   FILE
C
C       IFAIL (I/O) INTEGER: ON INPUT     =0, STOP ON OPEN FAILURE
C                                         =1, CONTINUE AFTER OPEN FAILURE
C                                             (only on file not found)
C                                         =2, CONTINUE SILENTLY AFTER OPEN FAILURE
C                                         =-1, As 0, but silent on success
C                                             (equivalent to negative IUN)
C                            ON OUTPUT    UNCHANGED IF FILE OPEN OK
C                                         =-1, ERROR IN OPENING FILE
C_END_CCPDPN
C
C     .. Scalar Arguments ..
      INTEGER IFAIL,IUN,IUN1,LREC
      CHARACTER LOGNAM* (*),STATUS* (*),TYPE* (*)
C     ..
C     .. Local Scalars ..
      INTEGER ISTAT,ITYPE
      CHARACTER ERRSTR*80
C     ..
C     .. Local Arrays ..
      CHARACTER TYPES(4)*2,STATS(6)*8, STAT*8, TYP*2
C     ..
C     .. External Functions ..
      INTEGER CCPNUN,LENSTR
      EXTERNAL CCPNUN,LENSTR
C     ..
C     .. External Subroutines ..
ccccc      EXTERNAL CCPOPN
C     ..
C     .. Data statements ..
      DATA STATS/'UNKNOWN','SCRATCH','OLD','NEW','READONLY','PRINTER'/
      DATA TYPES/'F','U','DF','DU'/
C     ..
C
      IF (IUN .EQ. 0) IUN = CCPNUN()
      STAT = STATUS
      TYP = TYPE
      CALL CCPUPC(STAT)
      CALL CCPUPC(TYP)
      DO 10 ISTAT = 1,6
        IF (STAT.EQ.STATS(ISTAT)) GO TO 20
   10 CONTINUE
      ERRSTR = ' CCPDPN: illegal status : '
      ERRSTR(LENSTR(ERRSTR)+2:) = STATUS
      CALL CCPERR(1,ERRSTR)
C
   20 DO 30 ITYPE = 1,4
        IF (TYP.EQ.TYPES(ITYPE)) GO TO 40
   30 CONTINUE
      ERRSTR = ' CCPDPN: illegal type: '
      ERRSTR(LENSTR(ERRSTR)+2:) = TYPE
      CALL CCPERR(1,ERRSTR)
C
 40   CONTINUE
      IUN1 = IUN
C  If IFAIL lt 0 No open message from CCPOPN
      IF(IFAIL.LT.0 .AND. IUN.GT.0) THEN
        IUN1 = -IUN
        IFAIL = 0
      ENDIF
      CALL CCPOPN_test(IUN1,LOGNAM,ISTAT,ITYPE,LREC,IFAIL)
C
      END
C

C ========
C UNIX.FOR
C ========
C
C Subroutines:
C
C CCPOPN - open a file
C UBYTES - Returns number of bytes per word and 'words'/'bytes'
C          to indicate if byte handling is available
C UGERR  - Get error explanation
C UGTENV - Get value of env. variable
C UGTIUD - Get user id - it's name
C UISATT - Is file a terminal?
C CCPSPW - Spawns a new process to run shell command
C CEXIT  - Trivial interface to system dependent EXIT routine 
C TTSEND - Write string to terminal with various carriage control
C     options
C UGTARG - Get command-line argument
C hciftime - Time in cif format
C ccp4_fflush_stdout - Flush buffers to stdout
C
C Functions:
C
C VAXVMS - Logical function returns TRUE if VAX/VMS
C WINMVS - Logical function returns TRUE if WINMVS
C RTNBKS - Returns backslash for Windows.
C
      SUBROUTINE CCPOPN_test(IIUN,LOGNAM,KSTAT,ITYPE,LREC,IFAIL)
C     ====================================================
C
C---- This subroutine is used to open a file
C
C     The requirement to specify that leading carriage control
C     characters in the output records should be obeyed (or not) can't
C     be implemented portably; likewise specifying readonly opening.
C     Some compilers accept VAXtran `carriagecontrol=' and `readonly'
C     specifiers; if so we use them.  Others have IOINIT, which can be
C     used to specify the carriage control.  The HPUX compiler is said
C     not to have any means of doing this and AIX seems to be likewise,
C     sigh; they both seem to obey the normal Unix convention of
C     printing the format as-is rather than obeying the first character
C     as carriage control.  Concentrix does obey the first column a la
C     VMS and `traditional' Fortran; the MIPS compilers have a compile
C     (link?) option to do so.  [Unfortunately, carriagecontrol
C     specification isn't even defined in Fortan90, although
C     `ACTION="READ"' can be used.]
C
C PARAMETERS
C ==========
C
C        IIUN (I)   UNIT NUMBER
C      LOGNAM (I)   LOGICAL FILE NAME (UP TO 8 CHARACTERS)
C       KSTAT (I)   FILE STATUS FLAG =1, 'UNKNOWN'
C                                    =2, 'SCRATCH'
C                                    =3, 'OLD'
C                                    =4, 'NEW'
C                                    =5, 'READONLY'
C                                    =6, 'PRINTER'
C       ITYPE (I)   FILE TYPE FLAG =1, 'SEQUENTIAL' 'FORMATTED'
C                                  =2, 'SEQUENTIAL' 'UNFORMATTED'
C                                  =3, 'DIRECT'     'FORMATTED'
C                                  =4, 'DIRECT'     'UNFORMATTED'
C        LREC (I)   RECORD LENGTH FOR DIRECT ACCESS FILE (NO. OF
C                   CHARACTERS FOR A FORMATTED FILE OR WORDS FOR
C                   AN UNFORMATTED FILE). NOT RELEVANT FOR A SEQUENTIAL
C                   FILE
C       IFAIL (I/O) ON INPUT:     =0, STOP ON OPEN FAILURE
C                                 =1, CONTINUE AFTER OPEN FAILURE
C                                 =2, CONTINUE SILENTLY AFTER OPEN FAILURE
C                   ON OUTPUT:    UNCHANGED IF FILE OPEN OK
C                                 =-1, ERROR IN OPENING FILE
C
C     .. Scalar Arguments ..
      INTEGER IFAIL,KSTAT,ITYPE,IIUN,LREC
      CHARACTER LOGNAM* (*)
C     ..
C     .. Local Scalars ..
      INTEGER LLREC,IUN,IBYTES,ISTAT,L,IOS
      CHARACTER CCNTRL*7,ST*7,FRM*12,ERRSTR*500,
     +     NAMFIL*255,HANDLE*5,OPNVAR*20, access*10
      INTEGER UNKNWN, SCRTCH, OLD, NEW, RDONLY, PRINTR
      PARAMETER (UNKNWN=1, SCRTCH=2, OLD=3, NEW=4, RDONLY=5, PRINTR=6)
ifdef(_ioinit,[      LOGICAL JUNK])dnl
C     ..
C     .. Local Arrays ..
      CHARACTER STAT(6)*7, DISP*6
C     ..
C     .. External Functions ..
      INTEGER LENSTR
ifdef(_ioinit,[
      LOGICAL IOINIT])dnl
      EXTERNAL LENSTR
C     ..
C     .. External Subroutines ..
      EXTERNAL UGERR,UGTENV
C     ..
C     .. Data statements ..
C     NB mustn't have SCRATCH in here, because result is system
C     -dependent
      DATA STAT/'UNKNOWN','UNKNOWN','OLD','NEW','OLD','UNKNOWN'/
C     ..
C     
      ISTAT = KSTAT
C     Negative unit number means don't give messages for successful open
      IUN = IIUN
      IF (IIUN.LT.0) IUN = -IIUN
C     Check args:
      IF (ISTAT.LT.1 .OR. ISTAT.GT.6 .OR. ITYPE.LT.1 .OR. ITYPE.GT.4)
     +     THEN 
        IF (IFAIL.EQ.0) THEN
          CALL CCPERR(1,
     +         '**CCPOPN ERROR** Invalid parameters in call')
        ELSE
          WRITE (6,
     +         '('' **CCPOPN ERROR** Invalid parameters in call'',/)')
          IFAIL = -1
        END IF
        RETURN
      ENDIF 
C
C     Do nothing for pre-connected units (what's the significance of
C     `TERM...'?) 
      IF (LOGNAM.EQ.'DATA' .OR. LOGNAM.EQ.'PRINTER' .OR.
     $     LOGNAM(:4).EQ.'TERM') RETURN
C
C     if environment variable CCP4_OPEN has value `UNKNOWN', open files
C     with status UNKNOWN rather than new if they exist
      IF (ISTAT.EQ.NEW) THEN
        OPNVAR = ' '
        CALL UGTENV('CCP4_OPEN',OPNVAR)
        IF (OPNVAR.EQ.'UNKNOWN') ISTAT = 1
      END IF
C
C     check for `logical name' referencing real file
      NAMFIL = ' '
      CALL UGTENV(LOGNAM,NAMFIL)
      IF (NAMFIL.EQ.' ') NAMFIL = LOGNAM

C     check for blank filename
      IF (NAMFIL.EQ.' ') THEN
        WRITE (ERRSTR,FMT=6001) IUN
 6001   FORMAT (' Open failed on unit ',I4,
     +          ': CCPOPN has received a blank filename.')
        CALL CCPERR(1, ERRSTR)
      ENDIF

C     VMS null device (VMS code canonicalises /dev/null)
      IF (NAMFIL.EQ.'NL:' .OR. NAMFIL.EQ.'nl:') NAMFIL='/dev/null'
C     Special case:  /dev/null should be opened UNKNOWN
      IF ( NAMFIL.EQ.'/dev/null') ISTAT = 1
C
C     type of open
      ST = STAT(ISTAT)
      IF (ITYPE.EQ.2 .OR. ITYPE.EQ.4) THEN
        FRM = 'UNFORMATTED'
      ELSE
        FRM = 'FORMATTED'
      ENDIF 
      IF (ITYPE .EQ. 1 .OR. ITYPE.EQ.2) THEN
        ACCESS='SEQUENTIAL'
      ELSE
        ACCESS='DIRECT'
      ENDIF
C
      IF (ISTAT.EQ.SCRTCH) THEN
        DISP = 'DELETE'
      ELSE
        DISP = 'KEEP'
      ENDIF
C     
      IF (access.eq.'DIRECT') THEN
C       Need to check is record length in words or bytes and set LLREC
C       accordingly. 
        CALL UBYTES (IBYTES,HANDLE)
        LLREC = LREC*IBYTES
        IF (HANDLE.EQ.'WORDS'.AND.ITYPE.EQ.4) LLREC=LLREC/IBYTES
        IF (ISTAT.EQ.RDONLY) THEN
C         _readonly may be defined as null or as `READONLY,'
          OPEN(UNIT=IUN,STATUS='UNKNOWN',ACCESS='DIRECT',FORM=FRM
     +         _readonly
     +         ,FILE=NAMFIL,RECL=LLREC,IOSTAT=IOS,ERR=5)
        ELSE
          OPEN(UNIT=IUN,STATUS='UNKNOWN',ACCESS='DIRECT',FORM=FRM
     +         _dispose
     +         ,FILE=NAMFIL,RECL=LLREC,IOSTAT=IOS,ERR=5)
        ENDIF 
      ELSE
C       if available, carriagecontrol='fortran' for print file, else = 
C       'list'.  we can use ioinit instead where it's available (see e.g.
C       Sun manual). 
        IF (ISTAT.EQ.PRINTR) THEN
C         want to obey format characters in column 1
          CCNTRL = 'FORTRAN'
          FRM = 'FORMATTED'
ifdef(_ioinit,
[      JUNK = IOINIT(.TRUE., .FALSE., .FALSE., ' ' , .FALSE.)
])dnl
        ELSE
C         no special significance to column 1
          CCNTRL = 'LIST'
ifdef(_ioinit,
[      JUNK = IOINIT(.FALSE., .FALSE., .FALSE., ' ' , .FALSE.)
])dnl
        END IF
        IF (FRM .EQ. 'UNFORMATTED') THEN
C         (carriage control not relevant) THINK THIS IS THE FAILURE POINT MRH
          IF (ISTAT.EQ.RDONLY) THEN
            OPEN(UNIT=IUN, FILE=NAMFIL, STATUS=ST, ACCESS='SEQUENTIAL'
     +           _readonly
     +           ,FORM=FRM, ERR=5, IOSTAT=IOS)
          ELSE
            OPEN(UNIT=IUN, FILE=NAMFIL, STATUS=ST, ACCESS='SEQUENTIAL'
     +           _dispose
     +           ,FORM=FRM, ERR=5, IOSTAT=IOS)
          ENDIF
        ELSE
          IF (ISTAT.EQ.RDONLY) THEN
            OPEN(UNIT=IUN, FILE=NAMFIL, STATUS=ST, ACCESS='SEQUENTIAL'
     +           _readonly
     +           _carriagecontrol
     +           ,FORM=FRM, ERR=5, IOSTAT=IOS)
          ELSE
            OPEN(UNIT=IUN, FILE=NAMFIL, STATUS=ST, ACCESS='SEQUENTIAL'
     +           _carriagecontrol
     +           _dispose
     +           ,FORM=FRM, ERR=5, IOSTAT=IOS)
          ENDIF
        ENDIF
      ENDIF
C
C     Scratch files are immediately unlinked from the directory; they
C     become inaccessible only when closed, but don't appear in the
C     directory and the name can be re-used.
C     NB this may break with REWIND if that is implemented as close +
C     reopen, sigh.  See also _dispose above
ifelse(_cant_unlink,1,,[
      IF (ISTAT.EQ.SCRTCH) CALL CUNLINK (NAMFIL)]
)dnl
C
C     Error check
 5    CONTINUE
C     don't report UNKNOWN if actually SCRATCH
      IF (ISTAT.EQ.SCRTCH) ST = 'SCRATCH'
      IF (IOS.NE.0) THEN
        CALL UGERR(IOS,ERRSTR)
        IF (IFAIL.EQ.0) THEN
C         hard failure
          WRITE (6,FMT=6002) IUN, NAMFIL(1:LENSTR(NAMFIL)),
     +         LOGNAM(1:LENSTR(LOGNAM))
 6002     FORMAT (' Open failed: Unit:',I4,', File: ',A, ' (logical: ',
     +         A, ')')
          ERRSTR = ' Open failed: File: ' // NAMFIL
          CALL CCPERR(-1, ERRSTR)
        else
C         soft failure
          IF (IFAIL.EQ.1) WRITE (6,FMT=6004) FRM, ST, IUN, 
     +         LOGNAM(1:LENSTR(LOGNAM)), NAMFIL(1:LENSTR(NAMFIL)),
     +         ERRSTR(1:LENSTR(ERRSTR))
 6004     FORMAT (' **CCPOPN ERROR**  ',A,3X,A,
     +         ' file open failure on unit ',I3,/' Logical name: ',
     +         A,', ','File name: ',A/1X,A/)
          IFAIL = -1
          RETURN            
        ENDIF
      ELSE
        IF (IIUN.LE.0) RETURN 
        WRITE (ERRSTR,FMT=6000) FRM,ST,IUN
        CALL QPRINT (1, ' ')
        CALL QPRINT (1, ERRSTR)
        call ccp4h_summary_beg()
        ERRSTR = 'Logical name: '
        ERRSTR (15:) = LOGNAM
        L = MIN(LENSTR (ERRSTR) + 1, LEN (ERRSTR))
        ERRSTR (L:) = ', Filename: ' // NAMFIL
        CALL QPRINT (1, ERRSTR)
        call ccp4h_summary_end()
        CALL QPRINT (1, ' ')
 6000 FORMAT (A,3X,A,' file opened on unit ',I3)
      ENDIF 
      END
C
C 
