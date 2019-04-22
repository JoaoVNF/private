* *******************************************************************
* COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
* Licence, see http://hsl.rl.ac.uk/acuk/cou.html
*
* Please note that for a UK ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between AEA
*    Technology plc and the Licensee on suitable terms and conditions,
*    which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor AEA Technology plc shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
*######DATE 2 Apr. 1993
C****************************************************************
C Exploit zeros in front by setting ISAVE(20) > 1
C (HSL 12 version equivalent to ISAVE(20) = 1)
C This version does row and col. swaps to exploit zeros
C in both rows and columns.
C To exploit zeros in front, set ISAVE(20)>1 before first call to
C MA42B/BD.
C Also, PIVBLK = ISAVE(19) in MA42F/FD allows us to wait
C until sufficiently many variables can be eliminated at once
C (as in cache paper). PIVBLK ALSO USED IN MA42J/JD.
C In HSL 12, PIVBLK = 1. If block size other than 1 is wanted,
C user must set ISAVE(19) to min. pivot block before
C calls to MA42J/JD (no action needed for min. block size of 1)
C but in experiments found  16 or 32 better (16 used in MA62)
C This version of code should be used with MP42 (which allows user to
C choose min. pivot block size)
C Note that for equation problems, 16 may be too large.
C USER MUST NOT CHANGE ISAVE(19) BETWEEN CALLS
C TO MA42J/JD AND CALLS TO MA42B/BD (IE USER CAN ONLY RESET
C ISAVE(19) AND ISAVE(20) ONCE) AND ISAVE(20)  MUST NOT BE CHANGED
C BETWEEN CALLS TO MA42B/BD
C****************************************************************
C Minor changes made to code 18 Nov. 1993.
C A bug was found in the case of the error flag INFO(1)=-14 being
C returned (KPVLNK was not properly set in this case and this caused
C unassigned variable error messages to be returned when LAST was
C restored).
C Also, ISAVE(38) set on each return from MA42A/AD, MA42J/JD, MA42P/PD,
C MA42B/BD to hold a copy of INFO(1). On each entry, a check is made
C that an error was not previously encountered.
C
C 4 March 1994 bug found in MA42P/PD. If ISTRM(2)=0 then we need to
C assign NUMBLK(2)=0 (we had NUMBLK(2) unassigned in this case).
C 15 April 1994 :error return -18 not possible from MA42J/JD so
C test for this removed from MA42B/BD. Also, before calling MA42J/JD
C from MA42B/BD for the first time, it is not necessary to move
C data to the start of IW so this has been removed).
C 18 April 1994: argument list for MA42D/DD changed to pass lengths
C of the array Y as LY1 and LY2. Y is of length 1,1 when called from
C MA42B/BD.
C 27 April 1994: change to MA42B/BD so that ISAVE(1)-ISAVE(15) are set
C           on first call (unless MA42P/PD has been called).
C           Want to be able to recall MA42B/BD without recalling
C           MA42P/PD.
C 13 May 1994: -15 cannot be returned by MA42O/OD or by MA42H/HD.
C 20 May 1994: INFO(1)=3 replaced by INFO(1)=2.
C 2 June 1994: NUMPIV made a local variable in MA42O/OD
C         KX made a local variable in MA42N/ND
C         In MA42J/JD, IFSIZE(1) = MAX(INFO(8),KFRNT)
C                      IFSIZE(2) = MAX(INFO(9),LFRNT)
C        (to stop code return a lower bound for the frontsize which is
C         smaller than that provided by user).
C         Changed printed error message when -12 returned.
C 13 June 1994: Immediate return added to MA42B/BD if it is entered with
C          INFO(1)<0.
C 24 June 1994: Bug found in MA42N/ND. We do not have the best pivot on
C          hold if IFORCE>0 and IELIM=1
C 1 Feb 1995: Operation count (RINFO(2)) changed since BLAS
C             cannot take advantage of zeros in frontal matrix.
C 17 March 1995   CLOSE (IFILE(I)) replaced by
C          CLOSE (IFILE(I),STATUS='DELETE')
C 6 August 1995  Ensured INFO(2) = 0 and RINFO(1) = 0.0 if matrix
C          found to be singular (this was not always happening
C          if the computation continued)
C 10 August 1995 Bug in MA42B/BD. IF (ISAVE(31).EQ.31) should read
C          IF (ISAVE(31).EQ.1)
C 14 August 1995 Bug corrected in computation of sign of
C          the determinant.
C 4 October 1995. Filled in buffers on the last time they are
C          written out... June 1996 changed fill in to avoid
C          doing it unnecessarily. We now only do it
C          if the buffer is the first and last one to
C          be written to d.a. file.
C 17 April 1996. Moved statement OFDIAG = 0 from MA42N/ND into
C          MA42F/FD (if the problem has only one element,
C          all variables are static condensation variables
C          and OFDIAG would not be set).
C 4  June 1996.  In MA42B/BD, avoid zeroing the solution vector X
C          if the number of variables in the problem is
C          equal to the largest integer used to index a variable
C          (i.e if NDF = INFO(3)). N.B. Cannot avoid the zeroing
C          in MA42C/CD since we do not pass the number of
C          variables to MA42C/CD in an ISAVE entry (clearly, we should
C          have done!). July 2000 : still have to zero X if
C          problem found to be singular, otherwise components
C          of X may be undefined.
C 11 Nov 1996 After each call to MA42D/DD and MA42E/ED, check
C          the error flag and write error message if appropriate.
C 20 Jan 1997 In MA42G/GD and MA42H/HD, JFLAG should only
C          be checked for an error IF MA42L/LD has been called
C          (o.w. could be undefined). This has been corrected.
C 14 March 1997 In MA42F/FD added a test so that MA42N/ND is only
C          called if some variables have become fully summed since
C          the assembly of the most recent elt/equ.
C 4 Sept. 1997  In MA42J/JD, if called from within MA42B/BD
C          (because frontsize too small), use INFO values to
C          initialise IFSIZE and allow updated lower bounds on the
C          filesizes to be returned to the user.
C 27 October 1998. In MA42J/JD and MA42B/BD changed test on NMAXE
C          so that NMAXE does not have to have the same
C          value on each entry if elements are being used
C          (still check that user has not changed for element
C          to equation entry, or visa versa, by comparing NMAXE with
C          ISAVE(21)
C          Changed  IF (NMAXE.NE.ISAVE(21)) GO TO 170
C          to   IF (ISAVE(21).EQ.1 .AND. NMAXE.GT.1) GO TO 170
C               IF (ISAVE(21).GT.1 .AND. NMAXE.EQ.1) GO TO 170
C          Similar change in MA42J/JD.
C 12 December 1998. INFO(3) is now updated in MA42F/FD (not in MA42H/HD).
C          Change made so that when MA42 used with MA52, INFO(3)
C          contains correct info. on exit (otherwise, as we do not
C          pick all variable as pivots within a subdomain, INFO(3)
C          would not be a count of all variables in subdomain but
C          would hold number of variables eliminated within subdomain)
C 18 December 1998. In MA42F/FD changed
C         INFO(10) = MKEY(1)
C         INFO(11) = MKEY(2)
C         INFO(12) = MKEY(3)
C to
C         INFO(10) = MAX(1,MKEY(1))
C         INFO(11) = MAX(1,MKEY(2)) ... this only if L factored stored
C         INFO(12) = MAX(1,MKEY(3))
C   (since otherwise can get INFO(10)-INFO(12) equal to zero when error
C    -17 return ... INFO values then no use for resetting LENFLE)
C 14 April 1999 In MA42C/CD, bug found in workspace LW if direct
C          access files used.
C 28 April 1999 Error in flop count for static condensation variables
C          corrected. In MA42O/OD
C          OPS = OPS + DBLE(NUMPIV* (NVAR-1)* (MVAR-1)*2)
C          changed to
C          DO 456 J = 1,NUMPIV
C              OPS = OPS + DBLE((NVAR-J)* (MVAR-J)*2)
C    456   CONTINUE
C          NOTE: has no effect for equation entry (since NUMPIV=1)
C 20 May 1999  Error found in MA42D/DD.
C              Replace
C               DO 10 I = 1,IRECD
C                  IBUFR(DIMIBF-IRECD+I) = IBUFR(IR2+I)
C   10          CONTINUE
C              by
C               DO 10 I = IRECD,1,-1
C                  IBUFR(DIMIBF-IRECD+I) = IBUFR(IR2+I)
C   10          CONTINUE
C              and
C               DO 30 I = 1,IPNT
C                  BUFR(DIMBUF-IPNT+I) = BUFR(JR2+I)
C   30          CONTINUE
C              by
C               DO 30 I = IPNT,1,-1
C                  BUFR(DIMBUF-IPNT+I) = BUFR(JR2+I)
C   30          CONTINUE
C 6 August 1999 Error in MA42D/DD. If direct access files not used
C               then JFLAG not set so at end of subroutine jump
C               to return (in MA42E/ED jump already there.)
C 29 November 1999. Error in call to MA42H/HD from MA42N/ND when
C KPRE>0 or LFRE>0.
C Also, changed second dimension of FA in MA42G/GD from NFRONT
C to * (NFRONT causes error when MA42G/GD called from MA42N/ND
C when KPRE>0 or LFRE>0). Similarly, second dimension of FRHS
C changed to *
C
C Changes made to incorporate the option of forcing diagonal pivots
C until the last element has been assembled.
C ICNTL(7)  used to control whether diagonal pivoting is forced.
C ICNTL(7) =  0 (default) Static condensation + off-diag. pivoting
C ICNTL(7) =  1001 Static condensation + diag. pivoting
C ICNTL(7) = -1001 No static condensation + diag. pivoting
C ALL other values of ICNTL(7):
C             No static condensation + off-diag. pivoting.
C Diagonal pivoting only allowed for ELEMENT entry.
C If diagonal pivoting and no warnings issued, on exit
C INFO(1) set to 6 + (Number of off-diagonal pivots)
C ISAVE(38) set to hold abs(ICNTL(7)) UNLESS an error is
C encountered, in which case it holds INFO(1).
C
C Subroutine MA42N/ND altered to cope with choosing diagonal pivots.
C Subroutine MA42O/OD altered to cope with choosing diagonal pivots.
C
C Also made a change to pivot search (on or off diagonal).
C In MA42N/ND  added a loop to look only at largest entry in
C column to to see if it can be used as a pivot
C (trying to reduce search time by not looking for largest entry in
C fully summed part unless we have to). This may mean INFO(13),
C INFO(14), and INFO(15) having different values than for previous
C version of code. In particular, some changes to test deck output.
C
C
C *** Reference MA42 suite ***
C *** Any problems contact Iain S. Duff  or Jennifer A. Scott
C     at Atlas Centre, Rutherford Appleton Laboratory ***
C *** Although every effort has been made to ensure robustness and
C     reliability of the subroutines in this MA42 suite,  we
C     disclaim any liability arising through the use or misuse of
C     any of the subroutines in the MA42 suite ***

      SUBROUTINE MA42ID(ICNTL,CNTL,ISAVE)
      DOUBLE PRECISION CNTL(2)
      INTEGER ICNTL(8),ISAVE(45)
      INTEGER I
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 8
      ICNTL(4) = 4
      ICNTL(5) = 0
      ICNTL(6) = 0
      ICNTL(7) = 0
      ICNTL(8) = 0
      CNTL(1) = 0.0D0
      CNTL(2) = 0.1D0
      DO 10 I = 1,45
         ISAVE(I) = 0
   10 CONTINUE
      DO 20 I = 1,3
         ISAVE(3+I) = 1
   20 CONTINUE
      ISAVE(16) = 2
      ISAVE(31) = -1
      ISAVE(19) = 1
      ISAVE(20) = 1
      RETURN
      END
C*******************************************************************
      SUBROUTINE MA42AD(NVAR,IVAR,NDF,LAST,LENLST,ICNTL,ISAVE,INFO)
      INTEGER LENLST,NDF,NVAR
      INTEGER ICNTL(8),INFO(23),ISAVE(45),IVAR(NVAR),LAST(LENLST)
      INTEGER I,JVAR,K,LP
      INTRINSIC ABS,MAX
      LP = ICNTL(1)
      IF (ISAVE(16).EQ.2) THEN
         NDF = 0
         INFO(1) = 0
         ISAVE(38) = ABS(ICNTL(7))
         ISAVE(16) = -2
         ISAVE(18) = 0
         IF (LENLST.LE.0) GO TO 30
         DO 10 I = 1,LENLST
            LAST(I) = 0
   10    CONTINUE
      END IF
      IF (ISAVE(38).LT.0) GO TO 100
      ISAVE(18) = ISAVE(18) + 1
      IF (NVAR.LE.0) GO TO 40
      DO 20 I = 1,NVAR
         JVAR = IVAR(I)
         IF (JVAR.LT.1 .OR. JVAR.GT.LENLST) GO TO 50
         NDF = MAX(JVAR,NDF)
         IF (LAST(JVAR).EQ.ISAVE(18)) GO TO 60
         LAST(JVAR) = ISAVE(18)
   20 CONTINUE
      ISAVE(30) = NDF
      GO TO 90
   30 INFO(1) = -1
      ISAVE(18) = 1
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),ISAVE(18)
         WRITE (LP,FMT=9010) LENLST
      END IF
      GO TO 90
   40 INFO(1) = -2
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),ISAVE(18)
         WRITE (LP,FMT=9020) NVAR
      END IF
      GO TO 90
   50 INFO(1) = -3
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),ISAVE(18)
         WRITE (LP,FMT=9030) I,JVAR
      END IF
      GO TO 90
   60 INFO(1) = -4
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),ISAVE(18)
         WRITE (LP,FMT=9040) JVAR
         DO 70 K = 1,I
            IF (IVAR(K).EQ.JVAR) GO TO 80
   70    CONTINUE
   80    WRITE (LP,FMT=9050) K,I
      END IF
   90 IF (INFO(1).LT.0) ISAVE(38) = INFO(1)
  100 RETURN
 9000 FORMAT (' ***** Error return from MA42A/AD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9010 FORMAT (7X,'Length (LENLST) of array LAST is ',I8)
 9020 FORMAT (7X,'Number of variables in elt/eqn is ',I8)
 9030 FORMAT (7X,'Variable ',I8,' in elt/eqn has value ',I8)
 9040 FORMAT (7X,'More than one occurrence of variable ',I8)
 9050 FORMAT (7X,'Dual occurrences in positions ',I8,' and ',I8)
      END
C**********************************************************************
      SUBROUTINE MA42JD(NVAR,IVAR,NDF,LAST,NMAXE,IFSIZE,ICNTL,ISAVE,
     +                  INFO)
      INTEGER NDF,NMAXE,NVAR
      INTEGER ICNTL(8),IFSIZE(5),INFO(23),ISAVE(45),IVAR(NVAR),LAST(NDF)
      INTEGER IELL,ISTATC,KFRNT,L,LFRNT,LK,LP,MFR,NELL
      INTEGER PIVBLK,KR,J,J1
      INTRINSIC ABS,MAX
      IF (ISAVE(38).LT.0) GO TO 120
      LP = ICNTL(1)
      IELL = ISAVE(17)
      NELL = ISAVE(18)
      IELL = IELL + 1
      KFRNT = ISAVE(39)
      LFRNT = ISAVE(40)
      PIVBLK = ABS(ISAVE(19))
      IF (ISAVE(16).EQ.-2) THEN
         INFO(1) = 0
         DO 10 L = 1,4
            IFSIZE(L) = 0
   10    CONTINUE
         IFSIZE(5) = 1
         ISAVE(16) = -1
         ISAVE(21) = NMAXE
         IF (NDF.NE.ISAVE(30)) GO TO 70
         IF (NMAXE.LE.0) GO TO 50
         ISAVE(41) = 0
      ELSE IF (ISAVE(16).EQ.4) THEN
         IFSIZE(1) = MAX(INFO(8),KFRNT)
         IFSIZE(2) = MAX(INFO(9),LFRNT)
         IFSIZE(3) = INFO(4)
         IFSIZE(4) = INFO(6)
         IFSIZE(5) = INFO(7)
         ISAVE(16) = 1
         ISAVE(41) = 0
      END IF
      IF (NVAR.LE.0) GO TO 30
      IF (ISAVE(21).EQ.1 .AND. NMAXE.GT.1) GO TO 50
      IF (ISAVE(21).GT.1 .AND. NMAXE.EQ.1) GO TO 50
      IF (NMAXE.EQ.1) KFRNT = KFRNT + 1
      KR = ISAVE(41)
      ISTATC = 0
      DO 20 LK = 1,NVAR
         MFR = IVAR(LK)
         IF (MFR.LT.1 .OR. MFR.GT.NDF) GO TO 40
         IF (LAST(MFR).GE.0) THEN
            IF (LAST(MFR).LT.IELL) GO TO 60
            IF (LAST(MFR).EQ.IELL) ISTATC = ISTATC + 1
            LFRNT = LFRNT + 1
            IF (LAST(MFR).GT.IELL) LAST(MFR) = -LAST(MFR)
         ELSE
            IF (-LAST(MFR).EQ.IELL) THEN
               LAST(MFR) = -LAST(MFR)
               KR = KR + 1
            END IF
         END IF
   20 CONTINUE
      IF (NMAXE.GT.1) THEN
         IFSIZE(2) = MAX(IFSIZE(2),LFRNT)
         IFSIZE(1) = IFSIZE(2)
         KFRNT = LFRNT
      ELSE IF (NMAXE.EQ.1) THEN
         IFSIZE(1) = MAX(IFSIZE(1),KFRNT)
         IFSIZE(2) = MAX(IFSIZE(2),LFRNT)
      END IF
      IF (ISTATC.GT.1 .AND. NMAXE.EQ.1) THEN
         KR = KR + ISTATC
         ISTATC = 0
      END IF
      IF (ISTATC.GT.0) THEN
         IF (ISTATC.EQ.1 .AND. NMAXE.EQ.1) THEN
            IFSIZE(5) = IFSIZE(5) + 5 + 1 + NVAR
            IFSIZE(3) = IFSIZE(3) + NVAR
            IFSIZE(4) = IFSIZE(4) + 1
            LFRNT = LFRNT - 1
            KFRNT = KFRNT - 1
         ELSE
            IFSIZE(5) = IFSIZE(5) + 5 + LFRNT + KFRNT
            LFRNT = LFRNT - ISTATC
            KFRNT = KFRNT - ISTATC
            IFSIZE(3) = IFSIZE(3) + ISTATC*LFRNT +
     +                  (ISTATC* (ISTATC+1))/2
            IFSIZE(4) = IFSIZE(4) + ISTATC*KFRNT +
     +                  (ISTATC* (ISTATC+1))/2
         END IF
      END IF
      IF (KFRNT.EQ.0 .OR. LFRNT.EQ.0) GO TO 24
      IF (KR.LT.PIVBLK .AND. IELL.LT.NELL) GO TO 25
      IF (KR.GT.0) THEN
         IFSIZE(5) = IFSIZE(5) + 5 + LFRNT + KFRNT
         J1 = KR
         DO 23 J = 1,J1
            KFRNT = KFRNT - 1
            LFRNT = LFRNT - 1
            KR = KR - 1
            IF (KFRNT.EQ.0 .OR. LFRNT.EQ.0) THEN
              IFSIZE(3) = IFSIZE(3) + J*LFRNT + (J* (J+1))/2
              IFSIZE(4) = IFSIZE(4) + J*KFRNT + (J* (J+1))/2
              GO TO 24
            END IF
   23    CONTINUE
         IFSIZE(3) = IFSIZE(3) + J1*LFRNT + (J1* (J1+1))/2
         IFSIZE(4) = IFSIZE(4) + J1*KFRNT + (J1* (J1+1))/2
      END IF
  24  CONTINUE
      IF (IELL.EQ.NELL) GO TO 100
  25  ISAVE(17) = IELL
      ISAVE(39) = KFRNT
      ISAVE(40) = LFRNT
      ISAVE(41) = KR
      GO TO 110
   30 INFO(1) = -2
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),IELL
         WRITE (LP,FMT=9010) NVAR
      END IF
      GO TO 80
   40 INFO(1) = -3
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),IELL
         WRITE (LP,FMT=9020) LK,MFR
      END IF
      GO TO 80
   50 INFO(1) = -8
      IF (LP.GT.0) THEN
         IF (NMAXE.LE.0) THEN
            WRITE (LP,FMT=9070) INFO(1)
            WRITE (LP,FMT=9060) NMAXE
         ELSE
            WRITE (LP,FMT=9000) INFO(1),IELL
            WRITE (LP,FMT=9040) ISAVE(21),NMAXE
         END IF
      END IF
      GO TO 80
   60 INFO(1) = -13
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),IELL
         WRITE (LP,FMT=9050) MFR
      END IF
      GO TO 80
   70 INFO(1) = -15
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9070) INFO(1)
         WRITE (LP,FMT=9030) ISAVE(30),NDF
      END IF
   80 DO 90 L = 1,NDF
         IF (LAST(L).LT.0) LAST(L) = -LAST(L)
   90 CONTINUE
  100 ISAVE(17) = 0
      ISAVE(39) = 0
      ISAVE(40) = 0
      ISAVE(16) = 2
  110 IF (INFO(1).LT.0) ISAVE(38) = INFO(1)
  120 RETURN
 9000 FORMAT (' ***** Error return from MA42J/JD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9010 FORMAT (7X,'Number of variables in elt/eqn is',I8)
 9020 FORMAT (7X,'Variable',I8,' in elt/eqn has value',I8)
 9030 FORMAT (' NDF has been changed from ',I8,' to ',I8)
 9040 FORMAT (7X,'NMAXE has been changed from',I8,' to ',I8)
 9050 FORMAT (7X,'Variable',I8,' is already fully summed')
 9060 FORMAT (7X,'NMAXE is equal to ',I8)
 9070 FORMAT (' ***** Error return from MA42J/JD *****  INFO(1) = ',I3)
      END
C*******************************************************************
      SUBROUTINE MA42PD(ISTRM,LENBUF,LENFLE,ICNTL,ISAVE,INFO)
      INTEGER ICNTL(8),INFO(23),ISAVE(45),ISTRM(3),LENBUF(3),LENFLE(3)
      INTEGER IFILE(3),ISIZE(3),NUMBLK(3)
      INTEGER I,ICON,IOS,LP,MP
      IF (ISAVE(38).LT.0) GO TO 110
      INFO(1) = 0
      ISAVE(31) = 1
      LP = ICNTL(1)
      MP = ICNTL(2)
      DO 10 I = 1,3
         IF (ISTRM(I).EQ.LP .OR. ISTRM(I).EQ.MP .OR. ISTRM(I).LT.0 .OR.
     +       ISTRM(I).GT.99 .OR. ISTRM(I).EQ.6) GO TO 90
         IF (I.NE.2 .AND. ISTRM(I).EQ.0) GO TO 90
         IFILE(I) = ISTRM(I)
         IF (IFILE(I).GT.0) THEN
            IF (LENBUF(I).LE.0) GO TO 40
            NUMBLK(I) = LENFLE(I)/LENBUF(I)
            IF (NUMBLK(I).EQ.0) GO TO 50
         ELSE
            LENBUF(I) = 0
            NUMBLK(I) = 0
         END IF
         ISIZE(I) = LENBUF(I)
   10 CONTINUE
      IF (IFILE(1).EQ.IFILE(2)) GO TO 60
      IF (IFILE(1).EQ.IFILE(3)) GO TO 60
      IF (IFILE(2).EQ.IFILE(3)) GO TO 60
      IF (ICNTL(3).LT.0 .OR. ICNTL(4).LT.0) GO TO 80
      DO 20 I = 1,3
         ICON = ICNTL(3)
         IF (I.EQ.3) ICON = ICNTL(4)
         IF (IFILE(I).GT.0) THEN
            CLOSE (IFILE(I),STATUS='DELETE')
            OPEN (IFILE(I),ERR=70,ACCESS='DIRECT',
     +           RECL=LENBUF(I)*ICON,IOSTAT=IOS)
         END IF
   20 CONTINUE
      DO 30 I = 1,3
         ISAVE(I) = IFILE(I)
         ISAVE(3+I) = 1
         ISAVE(6+I) = ISIZE(I)
         ISAVE(9+I) = 0
         ISAVE(12+I) = NUMBLK(I)
   30 CONTINUE
      ISAVE(16) = 2
      ISAVE(17) = 0
      ISAVE(39) = 0
      ISAVE(40) = 0
      GO TO 100
   40 INFO(1) = -19
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9010) I,LENBUF(I)
      END IF
      GO TO 100
   50 INFO(1) = -20
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9020) I,LENFLE(I),I,LENBUF(I)
      END IF
      GO TO 100
   60 INFO(1) = -21
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         IF (IFILE(1).EQ.IFILE(2)) WRITE (LP,FMT=9030) ISTRM(1)
         IF (IFILE(1).EQ.IFILE(3)) WRITE (LP,FMT=9040) ISTRM(1)
         IF (IFILE(2).EQ.IFILE(3)) WRITE (LP,FMT=9050) ISTRM(2)
      END IF
      GO TO 100
   70 INFO(1) = -22
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9060) IFILE(I),IOS
      END IF
      GO TO 100
   80 INFO(1) = -23
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9080) ICNTL(3),ICNTL(4)
      END IF
      GO TO 100
   90 INFO(1) = -24
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9070) I,ISTRM(I)
      END IF
  100 IF (INFO(1).LT.0) ISAVE(38) = INFO(1)
  110 RETURN
 9000 FORMAT (' ***** Error return from MA42P/PD ***** INFO(1) = ',I3)
 9010 FORMAT (7X,'Non-positive buffer size (LENBUF(',I1,')) of',I8)
 9020 FORMAT (7X,'File length (LENFLE(',I1,')) ',I8,/7X,'is less than ',
     +       'buffer size (LENBUF(',I1,')) of ',I8)
 9030 FORMAT (7X,'Stream numbers for L and U direct access data sets ',
     +       /7X,'both equal to ',I3)
 9040 FORMAT (7X,'Stream numbers for U and integer direct access data',
     +       ' sets',/7X,'both equal to ',I3)
 9050 FORMAT (7X,'Stream numbers for L and integer direct access data',
     +       ' sets',/7X,'both equal to ',I3)
 9060 FORMAT (7X,'Failure in direct access open',/7X,'file on stream',
     +       I3,5X,'IOSTAT = ',I10)
 9070 FORMAT (7X,'Illegal stream number. ISTRM(',I1,') is set to ',I3)
 9080 FORMAT (7X,'ICNTL(3) and ICNTL(4) are set to ',I3,',',I3)
      END
C**********************************************************************
      SUBROUTINE MA42BD(NVAR,IVAR,NDF,LAST,NMAXE,AVAR,NRHS,RHS,LRHS,LX,
     +                  X,NFRONT,LENBUF,LW,W,LIW,IW,ICNTL,CNTL,ISAVE,
     +                  INFO,RINFO)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER LIW,LRHS,LW,LX,NDF,NMAXE,NRHS,NVAR
      DOUBLE PRECISION AVAR(NMAXE,NVAR),CNTL(2),RHS(NMAXE,LRHS),
     +                 RINFO(2),W(LW),X(LX,LRHS)
      INTEGER ICNTL(8),INFO(23),ISAVE(45),IVAR(NVAR),IW(LIW),LAST(NDF),
     +        LENBUF(3),NFRONT(2)
      INTEGER DIMBUF,DIMIBF,FA,FRHS,I,IELL,JBUFL,JBUFR,JBUFU,K,KDEST,
     +        KFRNT,KHED,KPIV,KPVLNK,L,LDEST,LFL,LFRNT,LHED,LLIW,LLW,LP,
     +        MFR,NELL,NFR,PIVBLK,BZERO
      DOUBLE PRECISION Y(1,1)
      INTEGER IFILE(3),IREC(3),ISIZE(3),MKEY(3),NUMBLK(3)
      EXTERNAL MA42DD,MA42FD,MA42JD
      INTRINSIC ABS,MAX
      IELL = ISAVE(17)
      NELL = ISAVE(18)
      IELL = IELL + 1
      IF (ISAVE(38).LT.0) GO TO 280
      DO 10 I = 1,3
         IF (ISAVE(31).EQ.1) THEN
            IFILE(I) = ISAVE(I)
            IF (ICNTL(1).EQ.IFILE(I) .AND. IFILE(I).NE.0) ICNTL(1) = 6
            IF (ICNTL(2).EQ.IFILE(I) .AND. IFILE(I).NE.0) ICNTL(2) = 6
         END IF
   10 CONTINUE
      LP = ICNTL(1)
      IF (IELL.EQ.1) THEN
         IF (NMAXE.EQ.1 .AND. ABS(ICNTL(7)).EQ.1001) THEN
            IF (ICNTL(7).EQ.1001) ICNTL(7) = 0
            IF (ICNTL(7).EQ.-1001) ICNTL(7) = 1
         END IF
         ISAVE(38) = ABS(ICNTL(7))
         INFO(1) = 0
         INFO(2) = 1
         DO 20 I = 3,23
            INFO(I) = 0
   20    CONTINUE
         DO 30 I = 1,2
            RINFO(I) = ZERO
   30    CONTINUE
         DO 40 I = 1,3
            ISAVE(3+I) = 1
            ISAVE(9+I) = 0
   40    CONTINUE
         IF (ISAVE(31).EQ.-1) THEN
            DO 50 I = 1,3
               ISIZE(I) = LENBUF(I)
               ISAVE(6+I) = ISIZE(I)
               ISAVE(I) = 0
               ISAVE(12+I) = 0
   50       CONTINUE
            IF (LENBUF(1).LE.0 .OR. LENBUF(3).LE.0) GO TO 250
            IF (LENBUF(2).LT.0) GO TO 250
         ELSE
            DO 60 I = 1,3
               IF (LENBUF(I).NE.ISAVE(6+I)) GO TO 250
               ISIZE(I) = ISAVE(6+I)
   60       CONTINUE
         END IF
         IF (NRHS.GE.1 .AND. LX.LT.NDF) GO TO 140
         IF (LRHS.LT.NRHS .OR. LRHS.LT.1) GO TO 190
         IF (NDF.NE.ISAVE(30)) GO TO 230
         IF (NRHS.LT.0) GO TO 240
         IF (NMAXE.LE.0) GO TO 170
         ISAVE(19) = ABS(ISAVE(19))
         ISAVE(20) = ABS(ISAVE(20))
         IF (NMAXE.GT.1) NFRONT(2) = NFRONT(1)
         ISAVE(21) = NMAXE
         ISAVE(22) = NRHS
         ISAVE(23) = NFRONT(1)
         ISAVE(24) = NFRONT(2)
         JBUFU = ISIZE(1)
         JBUFL = MAX(1,ISIZE(2))
         JBUFR = 1 + ISIZE(2)
         FA = JBUFR + JBUFU
         FRHS = FA + NFRONT(1)*NFRONT(2)
         LLW = FRHS + MAX(NFRONT(1)*LRHS,NFRONT(2)*NRHS)
         IF (LW.LT.LLW) GO TO 150
         ISAVE(25) = JBUFU
         ISAVE(26) = JBUFL
         ISAVE(27) = JBUFR
         ISAVE(28) = FA
         ISAVE(29) = FRHS
         LHED = 1 + ISIZE(3)
         KHED = LHED + NFRONT(2)
         KPIV = KHED + NFRONT(1)
         KPVLNK = KPIV + NFRONT(2)
         LDEST = KPVLNK + NFRONT(2)
         KDEST = LDEST + NFRONT(2)
         LLIW = KDEST + NFRONT(1) - 1
         IF (LIW.LT.LLIW) GO TO 160
         ISAVE(32) = LHED
         ISAVE(33) = KHED
         ISAVE(34) = KPIV
         ISAVE(35) = KPVLNK
         ISAVE(36) = LDEST
         ISAVE(37) = KDEST
      END IF
      DO 70 I = 1,3
         IFILE(I) = ISAVE(I)
         IREC(I) = ISAVE(3+I)
         ISIZE(I) = ISAVE(6+I)
         MKEY(I) = ISAVE(9+I)
         NUMBLK(I) = ISAVE(12+I)
   70 CONTINUE
      IF (INFO(1).EQ.4) GO TO 120
      IF (NVAR.LE.0) GO TO 130
      IF (ISAVE(21).EQ.1 .AND. NMAXE.GT.1) GO TO 170
      IF (ISAVE(21).GT.1 .AND. NMAXE.EQ.1) GO TO 170
      IF (NMAXE.GT.1 .AND. NVAR.GT.NMAXE) GO TO 130
      IF (NRHS.NE.ISAVE(22)) GO TO 180
      IF (NFRONT(1).NE.ISAVE(23) .OR. NFRONT(1).LE.0) GO TO 200
      IF (NFRONT(2).NE.ISAVE(24) .OR. NFRONT(2).LE.0) GO TO 200
      IF (IELL.GT.1) THEN
         JBUFU = ISAVE(25)
         JBUFL = ISAVE(26)
         JBUFR = ISAVE(27)
         FA = ISAVE(28)
         FRHS = ISAVE(29)
         LHED = ISAVE(32)
         KHED = ISAVE(33)
         KPIV = ISAVE(34)
         KPVLNK = ISAVE(35)
         LDEST = ISAVE(36)
         KDEST = ISAVE(37)
      END IF
      KFRNT = ISAVE(39)
      LFRNT = ISAVE(40)
      PIVBLK = ABS(ISAVE(19))
      BZERO = ABS(ISAVE(20))
      CALL MA42FD(AVAR,RHS,LRHS,NRHS,IVAR,LAST,NDF,NMAXE,NVAR,NFRONT(1),
     +            NFRONT(2),W,JBUFL,W(JBUFR),JBUFU,IW,ISIZE(3),W(FA),
     +            W(FRHS),IW(LHED),IW(KHED),IW(KPIV),IW(KPVLNK),
     +            IW(LDEST),IW(KDEST),ICNTL,CNTL,IFILE,IREC,ISIZE,MKEY,
     +            NUMBLK,IELL,NELL,KFRNT,LFRNT,ISAVE(41),PIVBLK,
     +            BZERO,INFO,RINFO)
      IF (INFO(1).EQ.1 .OR. INFO(1).EQ.-14) THEN
         RINFO(1) = ZERO
         INFO(2) = 0
      END IF
      DO 80 I = 1,3
         ISAVE(3+I) = IREC(I)
         ISAVE(9+I) = MKEY(I)
   80 CONTINUE
      IF (INFO(1).LT.0) GO TO 260
      IF (INFO(1).EQ.4) THEN
         ISAVE(16) = 4
         LHED = ISAVE(32)
         NFR = 0
         DO 90 L = 1,LFRNT
            MFR = IW(LHED+L-1)
            IF (LAST(MFR).GE.IELL) THEN
               NFR = NFR + 1
               LAST(MFR) = -LAST(MFR)
            END IF
   90    CONTINUE
         KFRNT = KFRNT - (LFRNT-NFR)
         LFRNT = NFR
         ISAVE(39) = KFRNT
         ISAVE(40) = LFRNT
         GO TO 120
      END IF
      ISAVE(17) = IELL
      ISAVE(39) = KFRNT
      ISAVE(40) = LFRNT
      IF (INFO(1).EQ.5 .OR. INFO(1).EQ.6) GO TO 220
      IF (IELL.EQ.NELL) THEN
         IF (NRHS.GT.0) THEN
            IF (NDF.NE.INFO(3) .OR. INFO(2).EQ.0) THEN
               DO 110 L = 1,NRHS
                  DO 100 K = 1,NDF
                     X(K,L) = ZERO
  100             CONTINUE
  110          CONTINUE
            END IF
            IF (IFILE(1).EQ.0) THEN
               DIMBUF = JBUFU
            ELSE
               DIMBUF = JBUFU + NFRONT(1)*NFRONT(2)
            END IF
            IF (IFILE(3).EQ.0) THEN
               DIMIBF = ISIZE(3)
            ELSE
               DIMIBF = LIW
            END IF
            CALL MA42DD(1,NRHS,LX,X,.FALSE.,1,1,Y,W(JBUFR),DIMBUF,IW,
     +                  DIMIBF,W(FRHS),NFRONT(2),LP,NRHS,IFILE,IREC,
     +                  ISIZE,MKEY,INFO(1))
            IF (INFO(1).LT.0) THEN
               IF (LP.GE.0) WRITE (LP,FMT=9110) INFO(1)
               GO TO 260
            END IF
         END IF
         GO TO 260
      END IF
      GO TO 270
  120 CALL MA42JD(NVAR,IVAR,NDF,LAST,NMAXE,IW,ICNTL,ISAVE,INFO)
      IF (IELL.EQ.NELL) GO TO 210
      GO TO 270
  130 INFO(1) = -2
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9100) INFO(1),IELL
         WRITE (LP,FMT=9020) NVAR
         IF (NVAR.GT.NMAXE) WRITE (LP,FMT=9220) NMAXE
      END IF
      GO TO 260
  140 INFO(1) = -5
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9110) INFO(1)
         WRITE (LP,FMT=9040) LX,ISAVE(30)
      END IF
      GO TO 260
  150 INFO(1) = -6
      IELL = 1
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9100) INFO(1),IELL
         WRITE (LP,FMT=9030) LW,LLW
      END IF
      GO TO 260
  160 INFO(1) = -7
      IELL = 1
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9100) INFO(1),IELL
         WRITE (LP,FMT=9050) LIW,LLIW
      END IF
      GO TO 260
  170 INFO(1) = -8
      IF (LP.GT.0) THEN
         IF (NMAXE.LE.0) THEN
            WRITE (LP,FMT=9110) INFO(1)
            WRITE (LP,FMT=9060) NMAXE
         ELSE
            WRITE (LP,FMT=9100) INFO(1),IELL
            WRITE (LP,FMT=9000) ISAVE(21),NMAXE
         END IF
      END IF
      GO TO 260
  180 INFO(1) = -9
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9100) INFO(1),IELL
         WRITE (LP,FMT=9070) ISAVE(22),NRHS
      END IF
      GO TO 260
  190 INFO(1) = -10
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9110) INFO(1)
         WRITE (LP,FMT=9080) NRHS,LRHS
      END IF
      GO TO 260
  200 INFO(1) = -11
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9100) INFO(1),IELL
         IF (ISAVE(23).GT.0 .AND. ISAVE(24).GT.0) WRITE (LP,
     +       FMT=9090) ISAVE(23),ISAVE(24),NFRONT
         IF (ISAVE(23).LE.0) WRITE (LP,FMT=9170) ISAVE(23)
         IF (ISAVE(24).LE.0) WRITE (LP,FMT=9180) ISAVE(24)
      END IF
      GO TO 260
  210 INFO(1) = -12
      IF (LP.GT.0) THEN
         IF (NFRONT(1).GE.IW(1) .AND. NFRONT(2).GE.IW(2)) WRITE (LP,
     +       FMT=9270) NFRONT
      END IF
      NFRONT(1) = IW(1)
      NFRONT(2) = IW(2)
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9110) INFO(1)
         WRITE (LP,FMT=9130) NFRONT
      END IF
       INFO(4) = MAX(INFO(4),IW(3))
       INFO(5) = MAX(INFO(5),IW(3) + NRHS*NDF)
       INFO(6) = MAX(INFO(6),IW(4))
       INFO(7) = MAX(INFO(7),IW(5))
      GO TO 270
  220 IF (INFO(1).EQ.5 .AND. ISAVE(19).GT.0) ISAVE(19) = -ISAVE(19)
      IF (INFO(1).EQ.6 .AND. ISAVE(20).GT.0) ISAVE(20) = -ISAVE(20)
      IF (IELL.EQ.NELL) THEN
         IF (INFO(1).EQ.5) INFO(1) = -16
         IF (INFO(1).EQ.6) INFO(1) = -17
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9110) INFO(1)
            LFL = INFO(10)*ISIZE(1)
            IF (INFO(10).GT.NUMBLK(1) .AND. ISAVE(20).LT.0) WRITE (LP,
     +          FMT=9140) INFO(5),ISIZE(1),LFL
            IF (ISAVE(19).LT.0 .AND. INFO(5).GT.ISIZE(1)) WRITE (LP,
     +          FMT=9190) INFO(5),ISIZE(1),LFL
            LFL = INFO(11)*ISIZE(2)
            IF (INFO(11).GT.NUMBLK(2) .AND. ISAVE(20).LT.0) WRITE (LP,
     +          FMT=9150) INFO(6),ISIZE(2),LFL
            IF (ISAVE(19).LT.0 .AND. INFO(6).GT.ISIZE(2)) WRITE (LP,
     +          FMT=9200) INFO(6),ISIZE(2),LFL
            LFL = INFO(12)*ISIZE(3)
            IF (INFO(12).GT.NUMBLK(3) .AND. ISAVE(20).LT.0) WRITE (LP,
     +          FMT=9160) INFO(7),ISIZE(3),LFL
            IF (ISAVE(19).LT.0 .AND. INFO(7).GT.ISIZE(3)) WRITE (LP,
     +          FMT=9210) INFO(7),ISIZE(3),LFL
         END IF
         GO TO 260
      ELSE
         ISAVE(17) = IELL
         GO TO 270
      END IF
  230 INFO(1) = -15
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9110) INFO(1)
         WRITE (LP,FMT=9120) ISAVE(30),NDF
      END IF
      GO TO 260
  240 INFO(1) = -18
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9110) INFO(1)
         WRITE (LP,FMT=9010) NRHS
      END IF
      GO TO 260
  250 INFO(1) = -19
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9110) INFO(1)
         IF (ISAVE(31).EQ.1) THEN
            WRITE (LP,FMT=9250) I
         ELSE
            IF (LENBUF(1).LE.0) WRITE (LP,FMT=9230) LENBUF(1)
            IF (LENBUF(3).LE.0) WRITE (LP,FMT=9240) LENBUF(3)
            IF (LENBUF(2).LT.0) WRITE (LP,FMT=9260) LENBUF(2)
         END IF
      END IF
      GO TO 260
  260 ISAVE(17) = 0
      ISAVE(39) = 0
      ISAVE(40) = 0
      ISAVE(16) = 2
      ISAVE(19) = ABS(ISAVE(19))
      ISAVE(20) = ABS(ISAVE(20))
  270 IF (INFO(1).LT.0 .OR. INFO(1).GT.3) ISAVE(38) = INFO(1)
  280 RETURN
 9000 FORMAT (7X,'NMAXE has been changed from ',I8,' to ',I8)
 9010 FORMAT (7X,'Number of right hand sides (NRHS) is',I8)
 9020 FORMAT (7X,'Number of variables in elt/eqn is ',I8)
 9030 FORMAT (7X,'Length of real workspace too small.',/7X,'Increase ',
     +       'LW from ',I8,' to ',I8)
 9040 FORMAT (7X,'First dimension of array X too small.',/7X,'Increase',
     +       ' LX from ',I8,' to ',I8)
 9050 FORMAT (7X,'Length of integer workspace too small.',/7X,'Increas',
     +       'e LIW from ',I8,' to ',I8)
 9060 FORMAT (7X,'NMAXE is equal to ',I8)
 9070 FORMAT (7X,'Change in number of right-hand sides from ',I8,' to ',
     +       I8)
 9080 FORMAT (7X,'You wish to solve for ',I8,' right hand sides',/7X,
     +       'however arrays allow for ',I8,' right hand sides.',/7X,
     +       'Check values of NRHS and LRHS')
 9090 FORMAT (7X,'NFRONT(1) and NFRONT(2) have been changed',/7X,
     +       'from ',I8,' and ',I8,' to ',I8,' and ',I8,
     +       ' respectively.')
 9100 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9110 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3)
 9120 FORMAT (' NDF has been changed from ',I8,' to ',I8)
 9130 FORMAT (7X,'Lower bound on size of NFRONT required is ',2I8)
 9140 FORMAT (7X,'Total storage needed for U-factors is ',I8,/7X,'At ',
     +       'present LENBUF(1) size of ',I8,/7X,'LENFLE(1) must be at',
     +       ' least ',I8)
 9150 FORMAT (7X,'Total storage needed for L-factors is ',I8,/7X,'At ',
     +       'present LENBUF(2) size of ',I8,/7X,'LENFLE(2) must be at',
     +       ' least ',I8)
 9160 FORMAT (7X,'Total storage needed for integers is ',I8,/7X,'At ',
     +       'present LENBUF(3) size of ',I8,/7X,'LENFLE(3) must be at',
     +       ' least ',I8)
 9170 FORMAT (7X,'Front size (NFRONT(1)) non-positive equal to ',I8)
 9180 FORMAT (7X,'Front size (NFRONT(2)) non-positive equal to ',I8)
 9190 FORMAT (7X,'Total storage needed for U-factors is ',I8,/7X,'Set ',
     +       'LENBUF(1) to this or use direct access data sets.',/7X,
     +       'At present LENBUF(1) size of ',I8,/7X,'LENFLE(1) must be',
     +       ' at least ',I8)
 9200 FORMAT (7X,'Total storage needed for L-factors is ',I8,/7X,'Set ',
     +       'LENBUF(2) to this or use direct access data sets.',/7X,
     +       'At present LENBUF(2) size of ',I8,/7X,'LENFLE(2) must ',
     +       'be at least',I8)
 9210 FORMAT (7X,'Total storage needed for integers is ',I8,/7X,'Set ',
     +       'LENBUF(3) to this or use direct access data sets.',/7X,
     +       'At present LENBUF(3) size of ',I8,/7X,'LENFLE(3) must ',
     +       'be at least',I8)
 9220 FORMAT (7X,'But first dimension of element arrays is only ',I8)
 9230 FORMAT (7X,'Non-positive buffer size (LENBUF(1)) of ',I8)
 9240 FORMAT (7X,'Non-positive buffer size (LENBUF(3)) of ',I8)
 9250 FORMAT (7X,'LENBUF(',I1,') has been changed since the call to',
     +       ' MA42P/PD')
 9260 FORMAT (7X,'Negative buffer size (LENBUF(2)) of ',I8)
 9270 FORMAT (7X,'Frontsize of ',2I8,' is too small',/7X,'but unable',
     +       ' to provide more useful information on a suitable',/7X,
     +       'value for NFRONT.')
      END
C*********************************************************************
      SUBROUTINE MA42CD(TRANS,NRHS,LX,B,X,LW,W,LIW,IW,ICNTL,ISAVE,INFO)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER LIW,LW,LX,NRHS
      LOGICAL TRANS
      DOUBLE PRECISION B(LX,NRHS),W(LW),X(LX,NRHS)
      INTEGER ICNTL(8),INFO(23),ISAVE(45),IW(LIW)
      INTEGER IFILE(3),IREC(3),ISIZE(3),MKEY(3)
      INTEGER DIMIBF,I,L,L1,L2,L3,LBUFR,LLW,LP,NFMAX,NRHSB
      EXTERNAL MA42DD,MA42ED
      INTRINSIC MAX
      IF (ISAVE(38).LT.0) GO TO 100
      INFO(1) = 0
      DO 10 I = 1,3
         IFILE(I) = ISAVE(I)
         IF (ICNTL(1).EQ.IFILE(I) .AND. IFILE(I).NE.0) ICNTL(1) = 6
         IREC(I) = ISAVE(3+I)
         ISIZE(I) = ISAVE(6+I)
         MKEY(I) = ISAVE(9+I)
   10 CONTINUE
      LP = ICNTL(1)
      IF (NRHS.LT.1) GO TO 70
      NRHSB = ISAVE(22)
      IF (LX.LT.ISAVE(30)) GO TO 40
      DO 30 L = 1,NRHS
         DO 20 I = 1,ISAVE(30)
            X(I,L) = ZERO
   20    CONTINUE
   30 CONTINUE
      IF (ISIZE(2).EQ.0) GO TO 80
      L1 = ISIZE(1)
      L2 = ISIZE(2)
      L3 = ISIZE(3)
      NFMAX = MAX(INFO(8),INFO(9))
      IF (ISAVE(31).EQ.-1) THEN
         DIMIBF = L3
         IF (LIW.LT.DIMIBF) GO TO 60
         LLW = NRHS*NFMAX + L1 + L2
         IF (LW.LT.LLW) GO TO 50
         IF (TRANS) THEN
            CALL MA42ED(1,NRHS,LX,B,W(L2+1),L1,IW,DIMIBF,W(L2+L1+1),
     +                 NFMAX,LP,NRHSB,IFILE,ISIZE,MKEY,INFO(1))
            IF (INFO(1).LT.0) GO TO 90
            CALL MA42DD(2,NRHS,LX,X,.TRUE.,LX,NRHS,B,W,L2,IW,DIMIBF,
     +                 W(L2+L1+1),NFMAX,LP,NRHSB,IFILE,IREC,ISIZE,MKEY,
     +                 INFO(1))
            IF (INFO(1).LT.0) GO TO 90
         ELSE
            CALL MA42ED(2,NRHS,LX,B,W,L2,IW,DIMIBF,W(L2+L1+1),NFMAX,LP,
     +                 NRHSB,IFILE,ISIZE,MKEY,INFO(1))
            IF (INFO(1).LT.0) GO TO 90
            CALL MA42DD(1,NRHS,LX,X,.TRUE.,LX,NRHS,B,W(L2+1),L1,IW,
     +                 DIMIBF,W(L2+L1+1),NFMAX,LP,NRHSB,IFILE,IREC,
     +                 ISIZE,MKEY,INFO(1))
            IF (INFO(1).LT.0) GO TO 90
         END IF
      ELSE
         DIMIBF = L3 + 5 + INFO(8) + INFO(9)
         IF (LIW.LT.DIMIBF) GO TO 60
         LBUFR = MAX(L1,L2) + INFO(22)*MAX(INFO(8),INFO(9)+ISAVE(22))
         LLW = NRHS*NFMAX + LBUFR
         IF (LW.LT.LLW) GO TO 50
         IF (TRANS) THEN
            CALL MA42ED(1,NRHS,LX,B,W,LBUFR,IW,DIMIBF,W(LBUFR+1),NFMAX,
     +                 LP,NRHSB,IFILE,ISIZE,MKEY,INFO(1))
            IF (INFO(1).LT.0) GO TO 90
            CALL MA42DD(2,NRHS,LX,X,.TRUE.,LX,NRHS,B,W,LBUFR,IW,DIMIBF,
     +                 W(LBUFR+1),NFMAX,LP,NRHSB,IFILE,IREC,ISIZE,MKEY,
     +                 INFO(1))
            IF (INFO(1).LT.0) GO TO 90
         ELSE
            CALL MA42ED(2,NRHS,LX,B,W,LBUFR,IW,DIMIBF,W(LBUFR+1),NFMAX,
     +                 LP,NRHSB,IFILE,ISIZE,MKEY,INFO(1))
            IF (INFO(1).LT.0) GO TO 90
            CALL MA42DD(1,NRHS,LX,X,.TRUE.,LX,NRHS,B,W,LBUFR,IW,DIMIBF,
     +                 W(LBUFR+1),NFMAX,LP,NRHSB,IFILE,IREC,ISIZE,MKEY,
     +                 INFO(1))
            IF (INFO(1).LT.0) GO TO 90
         END IF
      END IF
      GO TO 100
   40 INFO(1) = -5
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9060) LX,ISAVE(30)
      END IF
      GO TO 100
   50 INFO(1) = -6
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9030) LW,LLW
      END IF
      GO TO 100
   60 INFO(1) = -7
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9020) LIW,DIMIBF
      END IF
      GO TO 100
   70 INFO(1) = -18
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9010) NRHS
      END IF
      GO TO 100
   80 INFO(1) = -27
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         IF (ISAVE(31).EQ.1) WRITE (LP,FMT=9040) IFILE(2)
         IF (ISAVE(31).EQ.-1) WRITE (LP,FMT=9050) ISIZE(2)
      END IF
      GO TO 100
   90 IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
  100 RETURN
 9000 FORMAT (' ***** Error return from MA42C/CD ***** INFO(1) = ',I3)
 9010 FORMAT (7X,'Number of right hand sides (NRHS) is ',I8)
 9020 FORMAT (7X,'Length of integer workspace too small.',/7X,'Increas',
     +       'e LIW from ',I8,' to ',I8)
 9030 FORMAT (7X,'Length of real workspace too small.',/7X,'Increase ',
     +       'LW from ',I8,' to ',I8)
 9040 FORMAT (7X,'No record of L-factor being stored.',/7X,'MA42P/PD ',
     +       'was called with ISTRM(2) set to ',I3)
 9050 FORMAT (7X,'No record of L-factor being stored.',/7X,'MA42B/BD ',
     +       'was called with LENBUF(2) set to ',I8)
 9060 FORMAT (7X,'First dimension of arrays B and X too small.',/7X,
     +       'Increase LX from ',I8,' to ',I8)
      END
C**********************************************************************
      SUBROUTINE MA42DD(IND,NRHS,NDF,X,LYINBF,LY1,LY2,Y,BUFR,DIMBUF,
     +                  IBUFR,DIMIBF,W,LDW,LP,NRHSB,IFILE,IREC,ISIZE,
     +                  MKEY,INFO)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      INTEGER DIMBUF,DIMIBF,IND,INFO,LDW,LP,LY1,LY2,NDF,NRHS,NRHSB
      LOGICAL LYINBF
      DOUBLE PRECISION BUFR(DIMBUF),W(LDW,NRHS),X(NDF,NRHS),Y(LY1,LY2)
      INTEGER IBUFR(DIMIBF),IFILE(3),IREC(3),ISIZE(3),MKEY(3)
      EXTERNAL DGEMM,DGEMV,DTPSV,MA42LD
      INTRINSIC ABS,MAX
      INTEGER FIL,I,IFIL,IKEY,ILNGTH,IPNT,IR1,IR2,IREAD,IRECD,J1,J2,
     +        JFLAG,JFRNT1,JR1,JR2,K,K1,K1PK,K2,K2PK,K3,K3PK,KEY,KFRNT,
     +        KFRNT1,KRO,L,LBUFR,LCO,LENGTH,LFRNT,LFRNT1,LIBUFR,NLOOP,
     +        NPIV,NPIV1,NREAD,NREC
      LBUFR = ISIZE(IND)
      LIBUFR = ISIZE(3)
      KEY = MKEY(IND)
      IKEY = MKEY(3)
      FIL = IFILE(IND)
      IFIL = IFILE(3)
      IF (FIL.NE.0 .AND. KEY.NE.0) THEN
         CALL MA42LD(2,FIL,KEY,BUFR(DIMBUF-LBUFR+1),LBUFR,IBUFR,LIBUFR,
     +               LP,JFLAG)
         IF (JFLAG.LT.0) GO TO 150
      END IF
      IF (IFIL.NE.0) THEN
         CALL MA42LD(-2,IFIL,IKEY,BUFR,LBUFR,IBUFR(DIMIBF-LIBUFR+1),
     +               LIBUFR,LP,JFLAG)
         IF (JFLAG.LT.0) GO TO 150
      ELSE
         IKEY = 0
      END IF
      JR2 = DIMBUF - LBUFR
      IR2 = DIMIBF - LIBUFR
      IPNT = IREC(IND) - 1
      IF (IREC(IND).EQ.1) IPNT = ISIZE(IND)
      IRECD = IREC(3) - 2
      NLOOP = MAX(1,MKEY(3))*LIBUFR
      DO 140 NREC = 1,NLOOP
         IF (IRECD.EQ.0 .AND. IKEY.EQ.0) GO TO 160
         IF (IRECD.EQ.0) THEN
            IF (IR2-LIBUFR.LT.0) IR2 = DIMIBF
            CALL MA42LD(-2,IFIL,IKEY,BUFR,LBUFR,IBUFR(IR2-LIBUFR+1),
     +                  LIBUFR,LP,JFLAG)
            IF (JFLAG.LT.0) GO TO 150
            IRECD = LIBUFR
            IR2 = IR2 - LIBUFR
         END IF
         ILNGTH = IBUFR(IR2+IRECD)
         NREAD = 0
         J1 = ABS(ILNGTH) - IRECD
         J2 = J1/LIBUFR
         IF (J1.GT.0) THEN
            NREAD = 1 + J2
            IF (J2*LIBUFR.EQ.J1) NREAD = NREAD - 1
            IF (IR2-NREAD*LIBUFR.LT.0) THEN
               DO 10 I = IRECD,1,-1
                  IBUFR(DIMIBF-IRECD+I) = IBUFR(IR2+I)
   10          CONTINUE
               IR2 = DIMIBF - IRECD
            END IF
            DO 20 IREAD = 1,NREAD
               CALL MA42LD(-2,IFIL,IKEY,BUFR,LBUFR,
     +                     IBUFR(IR2-IREAD*LIBUFR+1),LIBUFR,LP,JFLAG)
               IF (JFLAG.LT.0) GO TO 150
   20       CONTINUE
            IRECD = LIBUFR - (J1-J2*LIBUFR)
            IF (J2*LIBUFR.EQ.J1) IRECD = 0
         ELSE
            IRECD = IRECD - ABS(ILNGTH)
         END IF
         IR2 = IR2 - NREAD*LIBUFR
         IF (ILNGTH.LT.0) GO TO 140
         IR1 = IR2 + IRECD + 1
         NPIV = IBUFR(IR1+1)
         KFRNT = IBUFR(IR1+2)
         LFRNT = IBUFR(IR1+3)
         LFRNT1 = LFRNT - NPIV
         KFRNT1 = KFRNT - NPIV
         NPIV1 = (NPIV* (NPIV+1))/2
         IF (IND.EQ.1) THEN
            JFRNT1 = LFRNT1
            K1 = IR1 + 3 + KFRNT1
            K2 = IR1 + 3 + KFRNT
            K3 = IR1 + 3 + KFRNT + LFRNT1
            LENGTH = NPIV* (LFRNT1+NRHSB) + NPIV1
         ELSE
            JFRNT1 = KFRNT1
            K1 = IR1 + 3 + KFRNT + LFRNT1
            K2 = IR1 + 3
            K3 = IR1 + 3 + KFRNT1
            LENGTH = NPIV*KFRNT1 + NPIV1
         END IF
         NREAD = 0
         J1 = LENGTH - IPNT
         J2 = J1/LBUFR
         IF (J1.GT.0) THEN
            NREAD = 1 + J2
            IF (J2*LBUFR.EQ.J1) NREAD = NREAD - 1
            IF (JR2-NREAD*LBUFR.LT.0) THEN
               DO 30 I = IPNT,1,-1
                  BUFR(DIMBUF-IPNT+I) = BUFR(JR2+I)
   30          CONTINUE
               JR2 = DIMBUF - IPNT
            END IF
            DO 40 IREAD = 1,NREAD
               CALL MA42LD(2,FIL,KEY,BUFR(JR2-IREAD*LBUFR+1),LBUFR,
     +                     IBUFR,LIBUFR,LP,JFLAG)
               IF (JFLAG.LT.0) GO TO 150
   40       CONTINUE
            IPNT = LBUFR - (J1-J2*LBUFR)
            IF (J2*LBUFR.EQ.J1) IPNT = 0
         ELSE
            IPNT = IPNT - LENGTH
         END IF
         JR2 = JR2 - NREAD*LBUFR
         IF (LYINBF) THEN
            DO 60 L = 1,NRHS
               DO 50 K = 1,NPIV
                  K1PK = K1 + K
                  KRO = IBUFR(K1PK)
                  W(K,L) = Y(KRO,L)
   50          CONTINUE
   60       CONTINUE
         ELSE IF (IND.EQ.1) THEN
            JR1 = JR2 + IPNT + 1 + JFRNT1*NPIV + NPIV1
            DO 80 L = 1,NRHSB
               DO 70 K = 1,NPIV
                  W(K,L) = BUFR(JR1)
                  JR1 = JR1 + 1
   70          CONTINUE
   80       CONTINUE
         END IF
         IF (JFRNT1.NE.0) THEN
            DO 100 L = 1,NRHS
               DO 90 K = 1,JFRNT1
                  K2PK = K2 + K
                  LCO = IBUFR(K2PK)
                  W(K+NPIV,L) = X(LCO,L)
   90          CONTINUE
  100       CONTINUE
            JR1 = JR2 + IPNT + 1
            IF (NRHS.NE.1) CALL DGEMM('T','N',NPIV,NRHS,JFRNT1,-ONE,
     +                                BUFR(JR1),JFRNT1,W(1+NPIV,1),LDW,
     +                                ONE,W,LDW)
            IF (NRHS.EQ.1) CALL DGEMV('T',JFRNT1,NPIV,-ONE,BUFR(JR1),
     +                                JFRNT1,W(1+NPIV,1),1,ONE,W,1)
         END IF
         JR1 = JR2 + IPNT + 1 + JFRNT1*NPIV
         DO 110 L = 1,NRHS
            IF (IND.EQ.1) CALL DTPSV('U','T','U',NPIV,BUFR(JR1),W(1,L),
     +                               1)
            IF (IND.EQ.2) CALL DTPSV('U','T','N',NPIV,BUFR(JR1),W(1,L),
     +                               1)
  110    CONTINUE
         DO 130 L = 1,NRHS
            DO 120 K = 1,NPIV
               K3PK = K3 + K
               LCO = IBUFR(K3PK)
               X(LCO,L) = W(K,L)
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
      GO TO 160
  150 INFO = JFLAG
  160 RETURN
      END
C**********************************************************************
      SUBROUTINE MA42ED(IND,NRHS,NDF,R1,BUFR,DIMBUF,IBUFR,DIMIBF,W,LDW,
     +                  LP,NRHSB,IFILE,ISIZE,MKEY,INFO)
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
      INTEGER DIMBUF,DIMIBF,IND,INFO,LDW,LP,NDF,NRHS,NRHSB
      DOUBLE PRECISION BUFR(DIMBUF),R1(NDF,NRHS),W(LDW,NRHS)
      INTEGER IBUFR(DIMIBF),IFILE(3),ISIZE(3),MKEY(3)
      INTEGER FIL,I,IFIL,IKEY,ILNGTH,IP1,IR2,IREAD,IRECD,ISPACE,J1,J2,
     +        JFLAG,JFRNT,JFRNT1,JR2,JRECD,JSPACE,K,K1,K1PK,KEY,KFRNT,
     +        KFRNT1,KRO,L,LBUFR,LENGTH,LFRNT,LFRNT1,LIBUFR,NLOOP,NPIV,
     +        NPIV1,NREAD,NREC
      EXTERNAL DGEMM,DGEMV,DTPSV,MA42LD
      INTRINSIC ABS,MAX
      LBUFR = ISIZE(IND)
      LIBUFR = ISIZE(3)
      KEY = 1
      IKEY = 1
      FIL = IFILE(IND)
      IFIL = IFILE(3)
      IF (FIL.NE.0 .AND. MKEY(IND).NE.0) THEN
         CALL MA42LD(1,FIL,KEY,BUFR,LBUFR,IBUFR,LIBUFR,LP,JFLAG)
         IF (JFLAG.LT.0) GO TO 150
      END IF
      IF (IFIL.NE.0) THEN
         CALL MA42LD(-1,IFIL,IKEY,BUFR,LBUFR,IBUFR,LIBUFR,LP,JFLAG)
         IF (JFLAG.LT.0) GO TO 150
      END IF
      JRECD = 1
      IRECD = 1
      JR2 = LBUFR
      JSPACE = LBUFR
      IR2 = LIBUFR
      ISPACE = LIBUFR
      NLOOP = MAX(1,MKEY(3))*LIBUFR
      DO 110 NREC = 1,NLOOP
         NREAD = 0
         IF (ISPACE.EQ.0) THEN
            IF (IR2+LIBUFR.GT.DIMIBF) THEN
               IR2 = 0
               IRECD = 1
            END IF
            CALL MA42LD(-1,IFIL,IKEY,BUFR,LBUFR,IBUFR(IR2+1),LIBUFR,LP,
     +                  JFLAG)
            IF (JFLAG.LT.0) GO TO 150
            IR2 = IR2 + LIBUFR
            ISPACE = LIBUFR
         END IF
         ILNGTH = IBUFR(IRECD)
         IF (ILNGTH.EQ.0) GO TO 160
         J1 = ABS(ILNGTH) - ISPACE
         J2 = J1/LIBUFR
         IF (J1.GT.0) THEN
            NREAD = 1 + J2
            IF (J2*LIBUFR.EQ.J1) NREAD = NREAD - 1
            IF (IR2+NREAD*LIBUFR.GT.DIMIBF) THEN
               DO 10 I = 1,ISPACE
                  IBUFR(I) = IBUFR(IR2-ISPACE+I)
   10          CONTINUE
               IRECD = 1
               IR2 = ISPACE
            END IF
            DO 20 IREAD = 1,NREAD
               CALL MA42LD(-1,IFIL,IKEY,BUFR,LBUFR,
     +                     IBUFR(IR2+ (IREAD-1)*LIBUFR+1),LIBUFR,LP,
     +                     JFLAG)
               IF (JFLAG.LT.0) GO TO 150
   20       CONTINUE
         END IF
         IF (J1.LE.0) THEN
            ISPACE = ISPACE - ABS(ILNGTH)
         ELSE IF (J2*LIBUFR.EQ.J1) THEN
            ISPACE = 0
         ELSE
            ISPACE = LIBUFR - J1 + J2*LIBUFR
         END IF
         IR2 = IR2 + NREAD*LIBUFR
         NPIV = IBUFR(IRECD+1)
         KFRNT = IBUFR(IRECD+2)
         LFRNT = IBUFR(IRECD+3)
         KFRNT1 = KFRNT - NPIV
         LFRNT1 = LFRNT - NPIV
         NPIV1 = (NPIV* (NPIV+1))/2
         IF (IND.EQ.1) THEN
            JFRNT = LFRNT
            JFRNT1 = LFRNT1
            K1 = IRECD + 3 + KFRNT
            LENGTH = NPIV* (LFRNT1+NRHSB) + NPIV1
         ELSE
            JFRNT = KFRNT
            JFRNT1 = KFRNT1
            K1 = IRECD + 3
            LENGTH = NPIV*KFRNT1 + NPIV1
         END IF
         IF (ILNGTH.LT.0) GO TO 120
         NREAD = 0
         J1 = LENGTH - JSPACE
         J2 = J1/LBUFR
         IF (J1.GT.0) THEN
            NREAD = 1 + J2
            IF (J2*LBUFR.EQ.J1) NREAD = NREAD - 1
            IF (JR2+NREAD*LBUFR.GT.DIMBUF) THEN
               DO 30 I = 1,JSPACE
                  BUFR(I) = BUFR(JR2-JSPACE+I)
   30          CONTINUE
               JRECD = 1
               JR2 = JSPACE
            END IF
            DO 40 IREAD = 1,NREAD
               CALL MA42LD(1,FIL,KEY,BUFR(JR2+ (IREAD-1)*LBUFR+1),LBUFR,
     +                     IBUFR,LIBUFR,LP,JFLAG)
               IF (JFLAG.LT.0) GO TO 150
   40       CONTINUE
         END IF
         IF (J1.LE.0) THEN
            JSPACE = JSPACE - LENGTH
         ELSE IF (J2*LBUFR.EQ.J1) THEN
            JSPACE = 0
         ELSE
            JSPACE = LBUFR - J1 + J2*LBUFR
         END IF
         JR2 = JR2 + NREAD*LBUFR
         DO 60 L = 1,NRHS
            DO 50 K = 1,JFRNT
               K1PK = K1 + K
               KRO = IBUFR(K1PK)
               W(K,L) = R1(KRO,L)
   50       CONTINUE
   60    CONTINUE
         IP1 = JRECD + JFRNT1*NPIV
         DO 70 L = 1,NRHS
            IF (IND.EQ.1) CALL DTPSV('U','N','U',NPIV,BUFR(IP1),
     +                               W(JFRNT1+1,L),1)
            IF (IND.EQ.2) CALL DTPSV('U','N','N',NPIV,BUFR(IP1),
     +                               W(JFRNT1+1,L),1)
   70    CONTINUE
         IF (JFRNT1.NE.0) THEN
            IF (NRHS.NE.1) CALL DGEMM('N','N',JFRNT1,NRHS,NPIV,-ONE,
     +                                BUFR(JRECD),JFRNT1,W(JFRNT1+1,1),
     +                                LDW,ONE,W,LDW)
            IF (NRHS.EQ.1) CALL DGEMV('N',JFRNT1,NPIV,-ONE,BUFR(JRECD),
     +                                JFRNT1,W(JFRNT1+1,1),1,ONE,W,1)
         END IF
         DO 100 L = 1,NRHS
            DO 80 K = 1,JFRNT1
               K1PK = K1 + K
               KRO = IBUFR(K1PK)
               R1(KRO,L) = W(K,L)
   80       CONTINUE
            DO 90 K = JFRNT1 + 1,JFRNT
               K1PK = K1 + K
               KRO = IBUFR(K1PK)
               R1(KRO,L) = -W(K,L)
   90       CONTINUE
  100    CONTINUE
         JRECD = JRECD + LENGTH
         IRECD = IRECD + ABS(ILNGTH)
  110 CONTINUE
  120 DO 140 L = 1,NRHS
         DO 130 K = 1,JFRNT
            K1PK = K1 + K
            KRO = IBUFR(K1PK)
            R1(KRO,L) = ZERO
  130    CONTINUE
  140 CONTINUE
      GO TO 160
  150 INFO = JFLAG
  160 RETURN
      END
C**********************************************************************
      SUBROUTINE MA42FD(AVAR,RHS,LRHS,NRHS,IVAR,LAST,NDF,NMAXE,NVAR,
     +                  MFRONT,NFRONT,BUFRL,LLB,BUFRU,LUB,IBUFR,LIBUFR,
     +                  FA,FRHS,LHED,KHED,KPIV,KPVLNK,LDEST,KDEST,ICNTL,
     +                  CNTL,IFILE,IREC,ISIZE,MKEY,NUMBLK,IELL,NELL,
     +                  KFRNT,LFRNT,ICOM,PIVBLK,BZERO,INFO,RINFO)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER IELL,KFRNT,LFRNT,LIBUFR,LLB,LRHS,LUB,MFRONT,NDF,NELL,
     +        NFRONT,NMAXE,NRHS,NVAR,PIVBLK,BZERO
      DOUBLE PRECISION AVAR(NMAXE,NVAR),BUFRL(LLB),BUFRU(LUB),CNTL(2),
     +                 FA(MFRONT,NFRONT),FRHS(MFRONT,LRHS),
     +                 RHS(NMAXE,LRHS),RINFO(2)
      INTEGER IBUFR(LIBUFR),ICNTL(8),ICOM(5),IFILE(3),INFO(23),IREC(3),
     +        ISIZE(3),IVAR(NVAR),KDEST(MFRONT),KHED(MFRONT),
     +        KPIV(NFRONT),KPVLNK(NFRONT),LAST(NDF),LDEST(NFRONT),
     +        LHED(NFRONT),MKEY(3),NUMBLK(3)
      DOUBLE PRECISION DET,OPS,PIVOT
      INTEGER I,IFORCE,ISTATC,KPIVRX,KR,KS,L,LFREE,LFRNTT,LK,LP,LPIVCX,
     +        MFR,MP,MVAR,NELIM,OFDIAG
      INTEGER NSIZE,NEW
      LOGICAL LSTAT
      EXTERNAL MA42GD,MA42HD,MA42MD,MA42ND,MA42OD
      INTRINSIC ABS,MAX
      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (IFILE(1).EQ.MP .OR. IFILE(2).EQ.MP .OR. IFILE(3).EQ.MP) MP = 0
      IF (IELL.EQ.1) THEN
         ICOM(1) = 0
         ICOM(2) = 0
         ICOM(3) = 1
         ICOM(4) = 0
         ICOM(5) = 0
         DO 10 I = 1,NFRONT
            LDEST(I) = I + 1
   10    CONTINUE
      END IF
         LSTAT = .FALSE.
         IF (BZERO.GT.1) LSTAT = .TRUE.
      KR = ICOM(1)
      KS = ICOM(2)
      LFREE = ICOM(3)
      KPIVRX = ICOM(4)
      LPIVCX = ICOM(5)
      DET = RINFO(1)
      OPS = RINFO(2)
      OFDIAG = 0
      MVAR = NVAR
      IF (NMAXE.EQ.1) MVAR = 1
      IFORCE = 0
      NEW = 0
      ISTATC = 0
      DO 20 LK = 1,NVAR
         MFR = IVAR(LK)
         IF (MFR.GT.NDF .OR. MFR.LE.0) GO TO 50
         IF (LAST(MFR).GE.0) THEN
            IF (LAST(MFR).LT.IELL) GO TO 70
            NEW = NEW + 1
            IF (LAST(MFR).EQ.IELL) ISTATC = ISTATC + 1
         END IF
   20 CONTINUE
      IF (NMAXE.EQ.1 .AND. ISTATC.GT.1) THEN
         INFO(2) = 0
         DET = ZERO
         IF (ICNTL(8).EQ.0) GO TO 80
         IF (INFO(1).EQ.0) THEN
            INFO(1) = 1
            IF (MP.GT.0) THEN
               WRITE (MP,FMT=9060) INFO(1),IELL
               WRITE (MP,FMT=9040)
            END IF
         END IF
         ISTATC = 0
      END IF
      LFRNTT = LFRNT + NEW
      INFO(3) = INFO(3) + NEW
      IFORCE = MAX(0,LFRNTT-NFRONT)
      IF (NMAXE.EQ.1 .AND. KFRNT.EQ.MFRONT) IFORCE = MAX(IFORCE,1)
      IF (IFORCE.GT.KR) GO TO 60
      IF (IFORCE.GT.0 .AND. PIVBLK.GT.1) THEN
         NELIM = KR
         NSIZE = NELIM
         CALL MA42ND(IELL,NDF,LAST,MFRONT,NFRONT,FA,LHED,KHED,LRHS,NRHS,
     +               FRHS,KDEST,LDEST,KPIV,KPVLNK,KPIVRX,LPIVCX,PIVOT,
     +               KFRNT,LFRNT,KR,KS,NMAXE,NELIM,0,LIBUFR,IBUFR,
     +               BUFRL,LLB,BUFRU,LUB,IFILE,IREC,ISIZE,MKEY,NUMBLK,
     +               DET,OPS,NELL,CNTL,ICNTL,INFO,OFDIAG,NSIZE,LSTAT)
         IF (INFO(1).LT.0) GO TO 90
         LFRNTT = LFRNT + NEW
         IFORCE = MAX(0,LFRNTT-NFRONT)
         IF (NMAXE.EQ.1 .AND. KFRNT.EQ.MFRONT) IFORCE = MAX(IFORCE,1)
      END IF
      IF (IFORCE.GT.0) THEN
         NELIM = IFORCE
         NSIZE = NELIM
         CALL MA42ND(IELL,NDF,LAST,MFRONT,NFRONT,FA,LHED,KHED,LRHS,NRHS,
     +               FRHS,KDEST,LDEST,KPIV,KPVLNK,KPIVRX,LPIVCX,PIVOT,
     +               KFRNT,LFRNT,KR,KS,NMAXE,NELIM,IFORCE,LIBUFR,IBUFR,
     +               BUFRL,LLB,BUFRU,LUB,IFILE,IREC,ISIZE,MKEY,NUMBLK,
     +               DET,OPS,NELL,CNTL,ICNTL,INFO,OFDIAG,NSIZE,LSTAT)
         IF (INFO(1).LT.0) GO TO 90
         IF (INFO(1).EQ.4) GO TO 60
      END IF
      IF (ICNTL(7).NE.0 .AND. ICNTL(7).NE.1001) ISTATC = 0
      CALL MA42OD(IELL,AVAR,LRHS,NMAXE,RHS,NRHS,IVAR,NDF,LAST,NVAR,
     +            MFRONT,NFRONT,LHED,KHED,KR,KPIV,KPVLNK,LDEST,KDEST,
     +            ISTATC,LFREE,FRHS,BUFRL,LLB,BUFRU,LUB,IBUFR,LIBUFR,
     +            MVAR,KFRNT,LFRNT,FA,IFILE,IREC,ISIZE,MKEY,NUMBLK,DET,
     +            OPS,NELL,CNTL,ICNTL,INFO)
      IF (INFO(1).LT.0) GO TO 90
      CALL MA42MD(IELL,AVAR,LRHS,NMAXE,RHS,NRHS,IVAR,LAST,NDF,NVAR,
     +            MFRONT,NFRONT,FA,FRHS,LHED,KHED,KPIV,KPVLNK,LDEST,
     +            KDEST,ISTATC,LFREE,KR,KFRNT,LFRNT,INFO)
      IF (KR.EQ.0) GO TO 40
      IF (KR.LT.PIVBLK .AND. IELL.NE.NELL) GO TO 40
      IF (IFORCE.EQ.0 .AND. KR.EQ.ICOM(1) .AND. IELL.NE.NELL) GO TO 40
      NELIM = KR
      IFORCE = 0
      NSIZE = KR
      DO 30 I = 1,KR - NSIZE + 1
         NELIM = NSIZE
      CALL MA42ND(IELL,NDF,LAST,MFRONT,NFRONT,FA,LHED,KHED,LRHS,NRHS,
     +            FRHS,KDEST,LDEST,KPIV,KPVLNK,KPIVRX,LPIVCX,PIVOT,
     +            KFRNT,LFRNT,KR,KS,NMAXE,NELIM,IFORCE,LIBUFR,IBUFR,
     +            BUFRL,LLB,BUFRU,LUB,IFILE,IREC,ISIZE,MKEY,NUMBLK,DET,
     +            OPS,NELL,CNTL,ICNTL,INFO,OFDIAG,NSIZE,LSTAT)
      IF (INFO(1).LT.0) GO TO 90
   30 CONTINUE
   40 CONTINUE
      IF (IELL.EQ.NELL) THEN
         IF (KR.GT.0) THEN
            CALL MA42HD(-1,KFRNT,LFRNT,KR,IBUFR,LIBUFR,KHED,LHED,MFRONT,
     +                  NFRONT,LP,IFILE(3),IREC(3),ISIZE(3),MKEY(3),
     +                  NUMBLK(3),IELL,NELL,INFO)
            INFO(23) = KR
         END IF
         CALL MA42HD(2,1,1,0,IBUFR,LIBUFR,KHED,LHED,MFRONT,NFRONT,LP,
     +               IFILE(3),IREC(3),ISIZE(3),MKEY(3),NUMBLK(3),IELL,
     +               NELL,INFO)
         IF (INFO(1).LT.0) GO TO 90
         CALL MA42GD(1,2,FA,FRHS,NRHS,1,1,1,BUFRU,LUB,MFRONT,
     +               LP,IFILE(1),IREC(1),ISIZE(1),MKEY(1),
     +               NUMBLK(1),IELL,NELL,INFO)
         IF (INFO(1).LT.0) GO TO 90
         CALL MA42GD(2,2,FA,FRHS,NRHS,1,1,1,BUFRL,LLB,MFRONT,
     +               LP,IFILE(2),IREC(2),ISIZE(2),MKEY(2),
     +               NUMBLK(2),IELL,NELL,INFO)
         IF (INFO(1).LT.0) GO TO 90
         IF (INFO(16).GT.0) THEN
            IF (INFO(1).EQ.0 .OR. INFO(1).EQ.1) INFO(1) = 2
            IF (MP.GT.0) WRITE (MP,FMT=9000) INFO(16)
         END IF
         INFO(10) = MAX(1,MKEY(1))
         INFO(11) = MKEY(2)
         IF (IFILE(2).NE.0) INFO(11) = MAX(1,MKEY(2))
         INFO(12) = MAX(1,MKEY(3))
         IF (INFO(1).EQ.0 .AND. ABS(ICNTL(7)).EQ.1001) THEN
            IF (OFDIAG.GT.0) THEN
               INFO(1) = 6 + OFDIAG
               IF (MP.GT.0) WRITE (MP,FMT=9070) OFDIAG
            END IF
         END IF
      END IF
      GO TO 110
   50 INFO(1) = -3
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9050) INFO(1),IELL
         WRITE (LP,FMT=9010) LK,MFR
      END IF
      GO TO 90
   60 INFO(1) = 4
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9050) INFO(1),IELL
         WRITE (LP,FMT=9020) MFRONT,NFRONT
      END IF
      GO TO 90
   70 INFO(1) = -13
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9050) INFO(1),IELL
         WRITE (LP,FMT=9030) MFR
      END IF
      GO TO 90
   80 INFO(1) = -14
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9050) INFO(1),IELL
         WRITE (LP,FMT=9040)
      END IF
   90 DO 100 L = 1,LFRNT
         IF (KPVLNK(L).GT.0) GO TO 100
         MFR = LHED(L)
         LAST(MFR) = -KPVLNK(L)
  100 CONTINUE
  110 ICOM(1) = KR
      ICOM(2) = KS
      ICOM(3) = LFREE
      ICOM(4) = KPIVRX
      ICOM(5) = LPIVCX
      RINFO(1) = DET
      RINFO(2) = OPS
      RETURN
 9000 FORMAT (' ***** Warning from MA42B/BD *****  INFO(1) = 2',/7X,
     +       'Numerical criterion not satisfied by ',I8,' pivots')
 9010 FORMAT (7X,'Variable',I8,' in elt/eqn has value ',I8)
 9020 FORMAT (7X,'NFRONT not large enough, currently equals ',2I8)
 9030 FORMAT (7X,'Variable',I8,' is already fully summed')
 9040 FORMAT (7X,'Matrix found to be singular')
 9050 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9060 FORMAT (' ***** Warning from MA42B/BD *****  INFO(1) = ',I3,/7X,
     +       'after input of elt/equ ',I8)
 9070 FORMAT (' ***** Warning from MA42B/BD *****  ',/7X,I8,
     +       '  off-diagonal pivots were selected.')
      END
C**********************************************************************
      SUBROUTINE MA42GD(IND,JFL,FA,FRHS,NRHS,KFRNT,LFRNT,NPIV,BUFR,
     +                 LBUFR,MFRONT,LP,IFILE,IREC,ISIZE,MKEY,
     +                 NUMBLK,IELL,NELL,INFO)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER IELL,IFILE,IND,IREC,ISIZE,JFL,KFRNT,LBUFR,LFRNT,LP,
     +        MFRONT,MKEY,NELL,NPIV,NRHS,NUMBLK
      DOUBLE PRECISION BUFR(LBUFR),FA(MFRONT,*),FRHS(MFRONT,*)
      INTEGER INFO(23)
      INTEGER IBUFR(1)
      INTEGER I,I1,IBLOCK,IFRNT1,ISPACE,IUP,J,J1,J2,JBUFR,JFLAG,JFRNT1,
     +        K,K1,KFRNT1,LENGTH,LFRNT1,LIBUFR,NBUFR,NPIV1
      EXTERNAL MA42KD,MA42LD
      INTRINSIC MIN
      INTEGER NZEROL,NZEROU
      SAVE NZEROL,NZEROU
      IF (INFO(6).EQ.0) NZEROL = 0
      IF (INFO(4).EQ.0) NZEROU = 0
      IF (JFL.EQ.1) THEN
         KFRNT1 = KFRNT - NPIV
         LFRNT1 = LFRNT - NPIV
         NPIV1 = (NPIV* (NPIV+1))/2
         IF (IND.EQ.1) THEN
            IFRNT1 = KFRNT1
            JFRNT1 = LFRNT1
            LENGTH = NPIV* (LFRNT1+NRHS) + NPIV1
            INFO(4) = INFO(4) + NPIV*LFRNT1 + NPIV1
            INFO(5) = INFO(5) + LENGTH
         ELSE
            IFRNT1 = LFRNT1
            JFRNT1 = KFRNT1
            LENGTH = NPIV*KFRNT1 + NPIV1
            INFO(6) = INFO(6) + LENGTH
         END IF
      ELSE
         IF (IREC.EQ.1) THEN
            GO TO 110
         END IF
      END IF
      IF (ISIZE.EQ.0) THEN
         DO 20 I = 1,NPIV
            I1 = IFRNT1 + I
            DO 10 J = 1,JFRNT1 + I
               FA(J,I1) = ZERO
   10       CONTINUE
   20    CONTINUE
         GO TO 110
      END IF
      CALL MA42KD(JFL,IND,LENGTH,IFILE,ISIZE,MKEY,NUMBLK,IREC,IELL,NELL,
     +           LP,INFO)
      IF (INFO(1).EQ.5 .OR. INFO(1).EQ.6) GO TO 110
      LIBUFR = 1
      IF (JFL.EQ.1) THEN
         DO 80 IBLOCK = 1,3
            IF (IBLOCK.EQ.3 .AND. IND.EQ.2) GO TO 80
            IUP = NPIV
            IF (IBLOCK.EQ.3) IUP = NRHS
            DO 70 I = 1,IUP
               I1 = IFRNT1 + I
               ISPACE = ISIZE - IREC + 1
               IF (IBLOCK.EQ.1) THEN
                  J = 0
                  J2 = JFRNT1
               ELSE IF (IBLOCK.EQ.2) THEN
                  J = JFRNT1
                  J2 = I
               ELSE
                  J = IFRNT1
                  J2 = NPIV
               END IF
               J1 = J2 - ISPACE
               NBUFR = 1
               IF (J1.GT.0) NBUFR = 2 + (J1/ISIZE)
               DO 60 JBUFR = 1,NBUFR
                  K1 = ISIZE
                  IF (JBUFR.EQ.NBUFR) K1 = J1 - (J1/ISIZE)*ISIZE
                  IF (JBUFR.EQ.1) K1 = MIN(J2,ISPACE)
                  IF (IND.EQ.1) THEN
                     IF (IBLOCK.EQ.1 .OR. IBLOCK.EQ.2) THEN
                        DO 30 K = 1,K1
                           IF (FA(I1,J+K).EQ.ZERO) NZEROU = NZEROU + 1
                           BUFR(IREC+K-1) = FA(I1,J+K)
                           IF (J+K.NE.LFRNT1+I) FA(I1,J+K) = ZERO
   30                   CONTINUE
                     ELSE
                        DO 40 K = 1,K1
                           BUFR(IREC+K-1) = FRHS(J+K,I)
                           FRHS(J+K,I) = ZERO
   40                   CONTINUE
                     END IF
                  ELSE
                     DO 50 K = 1,K1
                        IF (FA(J+K,I1).EQ.ZERO) NZEROL = NZEROL + 1
                        BUFR(IREC+K-1) = FA(J+K,I1)
                        FA(J+K,I1) = ZERO
   50                CONTINUE
                  END IF
                  J = J + K1
                  IF (JBUFR.EQ.NBUFR .AND. IREC+K1.LE.ISIZE) THEN
                     IREC = IREC + K1
                  ELSE
                     MKEY = MKEY + 1
                     CALL MA42LD(3,IFILE,MKEY,BUFR,LBUFR,IBUFR,LIBUFR,
     +                          LP,JFLAG)
                     IF (JFLAG.LT.0) GO TO 100
                     IREC = 1
                  END IF
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
      ELSE
         IF (IREC.NE.1) THEN
            MKEY = MKEY + 1
            IF (IFILE.NE.0 .AND. MKEY.EQ.1) THEN
               DO 90 I = IREC,LBUFR
                  BUFR(I) = ZERO
   90          CONTINUE
            END IF
            IF (IFILE.NE.0) THEN
               CALL MA42LD(3,IFILE,MKEY,BUFR,LBUFR,IBUFR,LIBUFR,
     +                     LP,JFLAG)
               IF (JFLAG.LT.0) GO TO 100
            END IF
         END IF
      END IF
      GO TO 110
  100 INFO(1) = JFLAG
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
  110 RETURN
 9000 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3)
      END
C*********************************************************************
      SUBROUTINE MA42HD(JFL,KFRNT,LFRNT,NPIV,IBUFR,LIBUFR,KHED,LHED,
     +                 MFRONT,NFRONT,LP,IFILE,IREC,ISIZE,MKEY,NUMBLK,
     +                 IELL,NELL,INFO)
      INTEGER IELL,IFILE,IREC,ISIZE,JFL,KFRNT,LFRNT,LIBUFR,LP,MFRONT,
     +        MKEY,NELL,NFRONT,NPIV,NUMBLK
      INTEGER IBUFR(LIBUFR),INFO(23),KHED(MFRONT),LHED(NFRONT)
      INTEGER I,ISPACE,J,J1,J2,JBUFR,JFLAG,K,K1,LBUFR,LENGTH,NBUFR
      DOUBLE PRECISION BUFR(1)
      INTEGER ITEMP(4)
      EXTERNAL MA42KD,MA42LD
      INTRINSIC ABS,MAX,MIN
      IF (JFL.EQ.1) THEN
         LENGTH = 5 + LFRNT + KFRNT
      ELSE IF (JFL.EQ.-1) THEN
         LENGTH = - (5+LFRNT+KFRNT)
      ELSE
         LENGTH = 1
      END IF
      INFO(7) = INFO(7) + ABS(LENGTH)
      INFO(22) = MAX(INFO(22),NPIV)
      CALL MA42KD(JFL,3,ABS(LENGTH),IFILE,ISIZE,MKEY,NUMBLK,IREC,IELL,
     +           NELL,LP,INFO)
      IF (INFO(1).EQ.5 .OR. INFO(1).EQ.6) GO TO 90
      LBUFR = 1
      IF (ABS(JFL).EQ.1) THEN
         ITEMP(1) = LENGTH
         ITEMP(2) = NPIV
         ITEMP(3) = KFRNT
         ITEMP(4) = LFRNT
         DO 60 I = 1,4
            ISPACE = ISIZE - IREC + 1
            IF (I.EQ.1) J2 = 4
            IF (I.EQ.2) J2 = KFRNT
            IF (I.EQ.3) J2 = LFRNT
            IF (I.EQ.4) J2 = 1
            J1 = J2 - ISPACE
            NBUFR = 1
            IF (J1.GT.0) NBUFR = 2 + (J1/ISIZE)
            J = 0
            DO 50 JBUFR = 1,NBUFR
               K1 = ISIZE
               IF (JBUFR.EQ.NBUFR) K1 = J1 - (J1/ISIZE)*ISIZE
               IF (JBUFR.EQ.1) K1 = MIN(J2,ISPACE)
               IF (I.EQ.1) THEN
                  DO 10 K = 1,K1
                     IBUFR(IREC+K-1) = ITEMP(J+K)
   10             CONTINUE
               ELSE IF (I.EQ.2) THEN
                  DO 20 K = 1,K1
                     IBUFR(IREC+K-1) = KHED(J+K)
   20             CONTINUE
               ELSE IF (I.EQ.3) THEN
                  DO 30 K = 1,K1
                     IBUFR(IREC+K-1) = LHED(J+K)
   30             CONTINUE
               ELSE
                  DO 40 K = 1,K1
                     IBUFR(IREC) = LENGTH
   40             CONTINUE
               END IF
               J = J + K1
               IF (JBUFR.EQ.NBUFR .AND. IREC+K1.LE.ISIZE) THEN
                  IREC = IREC + K1
               ELSE IF (IFILE.EQ.0) THEN
                  INFO(1) = 5
                  IF (LP.GT.0) WRITE (LP,FMT=9010) INFO(1),IELL
                  IF (LP.GT.0) WRITE (LP,FMT=9020)
                  MKEY = MKEY + 1
                  IREC = 1
                  INFO(21) = 1
               ELSE
                  MKEY = MKEY + 1
                  CALL MA42LD(-3,IFILE,MKEY,BUFR,LBUFR,IBUFR,LIBUFR,LP,
     +                       JFLAG)
                  IF (JFLAG.LT.0) GO TO 80
                  IREC = 1
               END IF
   50       CONTINUE
   60    CONTINUE
      ELSE
         IBUFR(IREC) = 0
         MKEY = MKEY + 1
         IF (IFILE.NE.0 .AND. MKEY.EQ.1) THEN
            DO 70 I = IREC + 1,LIBUFR
               IBUFR(I) = -1
   70       CONTINUE
         END IF
         IREC = IREC + 1
         IF (IFILE.NE.0) THEN
            CALL MA42LD(-3,IFILE,MKEY,BUFR,LBUFR,IBUFR,LIBUFR,LP,JFLAG)
            IF (JFLAG.LT.0) GO TO 80
         END IF
      END IF
      GO TO 90
   80 INFO(1) = JFLAG
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
   90 RETURN
 9000 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3)
 9010 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9020 FORMAT (7X,'LENBUF(3) too small.',/7X,'Direct access data',
     +       ' sets not requested.',/7X,'Continue computation to find',
     +       ' space required.')
      END
C*******************************************************************
      SUBROUTINE MA42KD(JFL,IND,LENGTH,IFILE,ISIZE,MKEY,NUMBLK,IREC,
     +                  IELL,NELL,LP,INFO)
      INTEGER IELL,IFILE,IND,IREC,ISIZE,JFL,LENGTH,LP,MKEY,NELL,NUMBLK
      INTEGER INFO(23)
      INTEGER INFOIN,ISPACE,JFIT,JFIT1,NBUFR,NWRITE
      INTRINSIC ABS,MAX
      INFOIN = INFO(1)
      IF (ABS(JFL).EQ.1) THEN
         ISPACE = ISIZE - IREC + 1
         JFIT = LENGTH - ISPACE
         NBUFR = 1
         JFIT1 = JFIT/ISIZE
         IF (JFIT.GT.0) NBUFR = 2 + JFIT1
         IF (JFIT.EQ.0 .AND. IELL.LT.NELL) NBUFR = 2
         NWRITE = NBUFR
         IF (JFIT.GT.0 .AND. ISIZE*JFIT1.EQ.JFIT) NWRITE = NWRITE - 1
         IF (JFIT.EQ.0 .AND. IELL.LT.NELL) NWRITE = 1
         IF (IND.EQ.1) INFO(19) = MAX(INFO(19),NWRITE)
         IF (IND.EQ.2) INFO(20) = MAX(INFO(20),NWRITE)
         IF (IND.EQ.3) INFO(21) = MAX(INFO(21),NWRITE)
         IF (NBUFR.GT.1) THEN
            IF (IFILE.EQ.0) THEN
               INFO(1) = 5
            ELSE IF (MKEY+NBUFR-1.GT.NUMBLK) THEN
               INFO(1) = 6
            END IF
         END IF
         IF (INFO(1).EQ.5 .OR. INFO(1).EQ.6) THEN
            IF (NBUFR.EQ.1) THEN
               IREC = IREC + LENGTH
            ELSE
               MKEY = MKEY + NBUFR - 1
               IREC = 1 + JFIT - JFIT1*ISIZE
            END IF
         END IF
      ELSE
         IF (INFO(1).EQ.5) THEN
            MKEY = MKEY + 1
         ELSE IF (MKEY+1.GT.NUMBLK .AND. IFILE.NE.0) THEN
            INFO(1) = 6
            MKEY = MKEY + 1
         END IF
      END IF
      IF (INFO(1).NE.INFOIN .AND. LP.GT.0) THEN
         WRITE (LP,FMT=9010) INFO(1),IELL
         IF (INFO(1).EQ.5) WRITE (LP,FMT=9020) IND
         IF (INFO(1).EQ.6) WRITE (LP,FMT=9000) IND
      END IF
      RETURN
 9000 FORMAT (7X,'LENFLE(',I1,') too small.',/7X,'Continue computation',
     +       ' to find space required.')
 9010 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9020 FORMAT (7X,'LENBUF(',I1,') too small.',/7X,'Direct access data',
     +       ' sets not requested.',/7X,'Continue computation to find',
     +       ' space required.')
      END
C**********************************************************************
      SUBROUTINE MA42LD(IOPT,IFILE,MKEY,BUFR,LBUFR,IBUFR,LIBUFR,LP,
     +                  JFLAG)
      INTEGER IFILE,IOPT,JFLAG,LBUFR,LIBUFR,LP,MKEY
      DOUBLE PRECISION BUFR(LBUFR)
      INTEGER IBUFR(LIBUFR)
      INTEGER IOS
      INTRINSIC ABS
      JFLAG = 0
      IF (IFILE.EQ.0) GO TO 30
      IF (IOPT.EQ.1 .OR. IOPT.EQ.2) THEN
         READ (IFILE,REC=MKEY,ERR=10,IOSTAT=IOS) BUFR
      ELSE IF (IOPT.EQ.3) THEN
         WRITE (IFILE,REC=MKEY,ERR=20,IOSTAT=IOS) BUFR
      ELSE IF (IOPT.EQ.-1 .OR. IOPT.EQ.-2) THEN
         READ (IFILE,REC=MKEY,ERR=10,IOSTAT=IOS) IBUFR
      ELSE IF (IOPT.EQ.-3) THEN
         WRITE (IFILE,REC=MKEY,ERR=20,IOSTAT=IOS) IBUFR
      END IF
      IF (ABS(IOPT).EQ.1) MKEY = MKEY + 1
      IF (ABS(IOPT).EQ.2) MKEY = MKEY - 1
      GO TO 30
   10 JFLAG = -25
      IF (LP.GT.0) WRITE (LP,FMT=9000) MKEY,IFILE,IOS
      GO TO 30
   20 JFLAG = -26
      IF (LP.GT.0) WRITE (LP,FMT=9010) MKEY,IFILE,IOS
   30 RETURN
 9000 FORMAT (7X,'Error in direct access read.',/7X,'Record ',I8,' in ',
     +       'file ',I2,' IOSTAT ',I10)
 9010 FORMAT (7X,'Error in direct access write.',/7X,'Record ',I8,' in',
     +       ' file ',I2,' IOSTAT ',I10)
      END
C*********************************************************************
      SUBROUTINE MA42MD(IELL,AVAR,LRHS,NMAXE,RHS,NRHS,IVAR,LAST,NDF,
     +                  NVAR,MFRONT,NFRONT,FA,FRHS,LHED,KHED,KPIV,
     +                  KPVLNK,LDEST,KDEST,ISTATC,LFREE,KR,KFRNT,LFRNT,
     +                  INFO)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER IELL,ISTATC,KFRNT,KR,LFREE,LFRNT,LRHS,MFRONT,NDF,NFRONT,
     +        NMAXE,NRHS,NVAR
      DOUBLE PRECISION AVAR(NMAXE,NVAR),FA(MFRONT,NFRONT),
     +                 FRHS(MFRONT,LRHS),RHS(NMAXE,LRHS)
      INTEGER INFO(23),IVAR(NVAR),KDEST(MFRONT),KHED(MFRONT),
     +        KPIV(NFRONT),KPVLNK(NFRONT),LAST(NDF),LDEST(NFRONT),
     +        LHED(NFRONT)
      INTEGER IPOS,J,JVAR,K,KC1,KDST,KPOS,KR1,L,LDST,LELL,LFR,LPOS,MFR,
     +        NFREE,NUVAR
      INTRINSIC MAX,MIN
      IF (ISTATC.EQ.0) THEN
         IF (NMAXE.EQ.1) THEN
            KFRNT = KFRNT + 1
            KHED(KFRNT) = IELL
         END IF
         DO 10 L = 1,NVAR
            MFR = IVAR(L)
            IF (LAST(MFR).GE.0) THEN
               LFRNT = LFRNT + 1
               NFREE = LDEST(LFREE)
               LDEST(LFREE) = LFRNT
               LHED(LFRNT) = MFR
               IF (NMAXE.GT.1) THEN
                  KFRNT = KFRNT + 1
                  KHED(KFRNT) = MFR
                  KDEST(LFREE) = KFRNT
               END IF
               KPVLNK(LFRNT) = -LAST(MFR)
               LAST(MFR) = -LFREE
               LFREE = NFREE
            END IF
   10    CONTINUE
         KR1 = INFO(8) + 1
         KC1 = INFO(9) + 1
         IF (KR1.LE.KFRNT) THEN
            INFO(8) = KFRNT
            LFR = MAX(INFO(9),LFRNT)
            DO 40 K = KR1,KFRNT
               DO 20 L = 1,LFR
                  FA(K,L) = ZERO
   20          CONTINUE
               DO 30 J = 1,NRHS
                  FRHS(K,J) = ZERO
   30          CONTINUE
   40       CONTINUE
         END IF
         IF (KC1.LE.LFRNT) INFO(9) = LFRNT
         DO 60 L = KC1,LFRNT
            DO 50 K = 1,MIN(MFRONT,KR1)
               FA(K,L) = ZERO
   50       CONTINUE
   60    CONTINUE
      END IF
      NUVAR = NVAR - ISTATC
      DO 70 L = 1,NUVAR
         MFR = IVAR(L)
         IPOS = -LAST(MFR)
         LDST = LDEST(IPOS)
         IVAR(L) = LDST
   70 CONTINUE
      IF (NMAXE.EQ.1 .AND. ISTATC.EQ.1) GO TO 130
      IF (NMAXE.EQ.1) THEN
         DO 80 K = 1,NVAR
            LDST = IVAR(K)
            FA(KFRNT,LDST) = FA(KFRNT,LDST) + AVAR(1,K)
   80    CONTINUE
         DO 90 J = 1,NRHS
            FRHS(KFRNT,J) = FRHS(KFRNT,J) + RHS(1,J)
   90    CONTINUE
      ELSE
         DO 120 L = 1,NUVAR
            LDST = IVAR(L)
            JVAR = LHED(LDST)
            KPOS = -LAST(JVAR)
            KDST = KDEST(KPOS)
            DO 100 K = 1,NUVAR
               LDST = IVAR(K)
               FA(KDST,LDST) = FA(KDST,LDST) + AVAR(L,K)
  100       CONTINUE
            DO 110 J = 1,NRHS
               FRHS(KDST,J) = FRHS(KDST,J) + RHS(L,J)
  110       CONTINUE
  120    CONTINUE
      END IF
  130 DO 140 L = 1,NUVAR
         LDST = IVAR(L)
         JVAR = LHED(LDST)
         LELL = -KPVLNK(LDST)
         IF (LELL.LE.IELL) THEN
            KR = KR + 1
            KPIV(KR) = LDST
            KPVLNK(LDST) = KR
            LPOS = -LAST(JVAR)
            LAST(JVAR) = IELL
            LDEST(LPOS) = LFREE
            LFREE = LPOS
         END IF
         IVAR(L) = JVAR
  140 CONTINUE
      END
C**********************************************************************
      SUBROUTINE MA42ND(IELL,NDF,LAST,MFRONT,NFRONT,FA,LHED,KHED,LRHS,
     +                  NRHS,FRHS,KDEST,LDEST,KPIV,KPVLNK,KPIVRX,LPIVCX,
     +                  PIVOT,KFRNT,LFRNT,KR,KS,NMAXE,NELIM,IFORCE,
     +                  LIBUFR,IBUFR,BUFRL,LLB,BUFRU,LUB,IFILE,IREC,
     +                  ISIZE,MKEY,NUMBLK,DET,OPS,NELL,CNTL,ICNTL,INFO,
     +                  OFDIAG,NSIZE,LSTAT)
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
      DOUBLE PRECISION DET,OPS,PIVOT
      INTEGER IELL,IFORCE,KFRNT,KPIVRX,KR,KS,LFRNT,LIBUFR,LLB,LPIVCX,
     +        LRHS,LUB,MFRONT,NDF,NELIM,NELL,NFRONT,NMAXE,NRHS,NSIZE,
     +        OFDIAG
      LOGICAL LSTAT
      DOUBLE PRECISION BUFRL(LLB),BUFRU(LUB),CNTL(2),FA(MFRONT,NFRONT),
     +                 FRHS(MFRONT,LRHS)
      INTEGER IBUFR(LIBUFR),ICNTL(8),IFILE(3),INFO(23),IREC(3),ISIZE(3),
     +        KDEST(MFRONT),KHED(MFRONT),KPIV(NFRONT),KPVLNK(NFRONT),
     +        LAST(NDF),LDEST(NFRONT),LHED(NFRONT),MKEY(3),NUMBLK(3)
      DOUBLE PRECISION GMON,GMONX,PIVINV,PIVX,RMAX,SWAP
      INTEGER CHFRNT,I,IAUTO,IELIM,ISRCH,J,JJ,JSRCH,K,KFR,KFRO,KH,KK,
     +        KMAX,KPIVKR,KPIVRO,KPOS,KPRE,KRTEMP,KX,L,LCOL,LFR,LFRO,LK,
     +        LL,LP,LPIVCO,LPOS,LPRE,LTEMP,MFR,MP,L1
      INTEGER IDAMAX
      EXTERNAL IDAMAX
      EXTERNAL DAXPY,DGEMM,DGER,DTRSM,MA42GD,MA42HD
      INTRINSIC ABS,DBLE,LOG,MIN
      IF (NELIM.EQ.0) GO TO 420
      LP = ICNTL(1)
      MP = ICNTL(2)
      IAUTO = ICNTL(5)
      ISRCH = ICNTL(6)
      KFR = KFRNT
      LFR = LFRNT
      KRTEMP = KR
      DO 140 IELIM = 1,NELIM
         IF (KFR.EQ.0 .OR. LFR.EQ.0) THEN
            INFO(2) = 0
            DET = ZERO
            IF (ICNTL(8).EQ.0) GO TO 410
            IF (INFO(1).EQ.0) THEN
               INFO(1) = 1
               IF (MP.GT.0) THEN
                  WRITE (MP,FMT=9020) INFO(1),IELL
                  WRITE (MP,FMT=9000)
               END IF
            END IF
            GO TO 150
         END IF
         JSRCH = KR
         IF (IAUTO.GT.0 .AND. KR.GE.IAUTO .AND.
     +       ISRCH.GT.0) JSRCH = MIN(KR,ISRCH)
         IF (ABS(ICNTL(7)).NE.1001) GO TO 20
         GMONX = ZERO
         KX = 0
         DO 10 L = 1,JSRCH
            KS = KS + 1
            IF (KS.GT.KR) KS = 1
            LPIVCO = KPIV(KS)
            KPIVRO = LPIVCO
            PIVOT = FA(KPIVRO,LPIVCO)
            INFO(13) = INFO(13) + 1
            INFO(15) = INFO(15) + KFR
            KMAX = IDAMAX(KFR,FA(1,LPIVCO),1)
            RMAX = ABS(FA(KMAX,LPIVCO))
            IF (RMAX.LE.CNTL(1)) THEN
               INFO(2) = 0
               DET = ZERO
               IF (ICNTL(8).EQ.0) GO TO 410
               IF (INFO(1).EQ.0) THEN
                  INFO(1) = 1
                  IF (MP.GT.0) THEN
                     WRITE (MP,FMT=9020) INFO(1),IELL
                     WRITE (MP,FMT=9000)
                  END IF
               END IF
               GO TO 10
            END IF
            INFO(14) = INFO(14) + 1
            KH = KX
            KX = KS
            GMON = ABS(PIVOT)/RMAX
            IF (GMON.GT.CNTL(2)) GO TO 60
            KX = KH
            IF (GMON.GT.GMONX) THEN
               PIVX = PIVOT
               KX = KS
               GMONX = GMON
               LPIVCX = LPIVCO
               KPIVRX = KPIVRO
            END IF
   10    CONTINUE
         IF (IELL.EQ.NELL .AND. IFORCE.EQ.0) GO TO 20
         IF (ABS(PIVOT).LE.CNTL(1)) THEN
            IF (IFORCE.EQ.0) GO TO 150
            INFO(1) = 4
            GO TO 420
         ELSE IF (IFORCE.EQ.0 .AND. (IAUTO.EQ.0.OR.KR.LT.IAUTO)) THEN
            GO TO 150
         ELSE
            INFO(16) = INFO(16) + 1
            PIVOT = PIVX
            LPIVCO = LPIVCX
            KPIVRO = KPIVRX
            GO TO 60
         END IF
   20    CONTINUE
         GMONX = ZERO
         KX = 0
         PIVOT = ZERO
         DO 30 L = 1,JSRCH
            KS = KS + 1
            IF (KS.GT.KR) KS = 1
            LPIVCO = KPIV(KS)
            INFO(13) = INFO(13) + 1
            INFO(15) = INFO(15) + KFR
            KMAX = IDAMAX(KFR,FA(1,LPIVCO),1)
            RMAX = ABS(FA(KMAX,LPIVCO))
            IF (RMAX.LE.CNTL(1)) THEN
               INFO(2) = 0
               DET = ZERO
               IF (ICNTL(8).EQ.0) GO TO 410
               IF (INFO(1).EQ.0) THEN
                  INFO(1) = 1
                  IF (MP.GT.0) THEN
                     WRITE (MP,FMT=9020) INFO(1),IELL
                     WRITE (MP,FMT=9000)
                  END IF
               END IF
               GO TO 30
            END IF
            KH = KX
            KX = KS
            MFR = KHED(KMAX)
            IF (NMAXE.EQ.1 .OR. LAST(MFR).GE.0) THEN
               KPIVRO = KMAX
               PIVOT = FA(KPIVRO,LPIVCO)
               GO TO 60
            END IF
   30    CONTINUE
         DO 50 L = 1,JSRCH
            KS = KS + 1
            IF (KS.GT.KR) KS = 1
            LPIVCO = KPIV(KS)
            INFO(13) = INFO(13) + 1
            INFO(15) = INFO(15) + KFR
            KMAX = IDAMAX(KFR,FA(1,LPIVCO),1)
            RMAX = ABS(FA(KMAX,LPIVCO))
            IF (RMAX.LE.CNTL(1)) THEN
               INFO(2) = 0
               DET = ZERO
               IF (ICNTL(8).EQ.0) GO TO 410
               IF (INFO(1).EQ.0) THEN
                  INFO(1) = 1
                  IF (MP.GT.0) THEN
                     WRITE (MP,FMT=9020) INFO(1),IELL
                     WRITE (MP,FMT=9000)
                  END IF
               END IF
               GO TO 50
            END IF
            KH = KX
            KX = KS
            MFR = KHED(KMAX)
            IF (NMAXE.EQ.1 .OR. LAST(MFR).GE.0) THEN
               KPIVRO = KMAX
               PIVOT = FA(KPIVRO,LPIVCO)
               GO TO 60
            END IF
            PIVOT = ZERO
            DO 40 K = 1,KFR
               MFR = KHED(K)
               IF (LAST(MFR).GE.0) THEN
                  INFO(14) = INFO(14) + 1
                  IF (ABS(FA(K,LPIVCO)).GE.ABS(PIVOT)) THEN
                     KPIVRO = K
                     PIVOT = FA(KPIVRO,LPIVCO)
                  END IF
               END IF
   40       CONTINUE
            GMON = ABS(PIVOT)/RMAX
            IF (GMON.GT.CNTL(2)) GO TO 60
            KX = KH
            IF (GMON.GT.GMONX) THEN
               PIVX = PIVOT
               KX = KS
               GMONX = GMON
               LPIVCX = LPIVCO
               KPIVRX = KPIVRO
            END IF
   50    CONTINUE
         IF (ABS(PIVOT).LE.CNTL(1)) THEN
            IF (IFORCE.EQ.0) GO TO 150
            INFO(1) = 4
            GO TO 420
         ELSE IF (IFORCE.EQ.0 .AND. (IAUTO.EQ.0.OR.KR.LT.IAUTO)) THEN
            GO TO 150
         ELSE
            INFO(16) = INFO(16) + 1
            PIVOT = PIVX
            LPIVCO = LPIVCX
            KPIVRO = KPIVRX
         END IF
   60    CONTINUE
         IF (PIVOT.LT.ZERO) INFO(2) = -INFO(2)
         IF (KPIVRO.NE.LPIVCO) OFDIAG = OFDIAG + 1
         IF (KPVLNK(LFR).GT.0) KX = KPVLNK(LFR)
         KPIVKR = KPIV(KR)
         KPIV(KX) = KPIVKR
         KPVLNK(KPIVKR) = KX
         IF (KPVLNK(LFR).LT.0) KPVLNK(LPIVCO) = KPVLNK(LFR)
         LFRO = LHED(LFR)
         LPOS = -LAST(LFRO)
         IF (LPOS.GT.0) LDEST(LPOS) = LPIVCO
         IF (NMAXE.GT.1) THEN
            KFRO = KHED(KFR)
            KPOS = -LAST(KFRO)
            IF (KPOS.GT.0) KDEST(KPOS) = KPIVRO
         END IF
         IF (KFR.NE.KPIVRO) THEN
            LTEMP = KHED(KFR)
            KHED(KFR) = KHED(KPIVRO)
            KHED(KPIVRO) = LTEMP
            DO 90 L = 1,LFRNT
               SWAP = FA(KFR,L)
               FA(KFR,L) = FA(KPIVRO,L)
               FA(KPIVRO,L) = SWAP
   90       CONTINUE
            INFO(2) = -INFO(2)
         END IF
         IF (LFR.NE.LPIVCO) THEN
            LTEMP = LHED(LFR)
            LHED(LFR) = LHED(LPIVCO)
            LHED(LPIVCO) = LTEMP
            DO 100 K = 1,KFRNT
               SWAP = FA(K,LFR)
               FA(K,LFR) = FA(K,LPIVCO)
               FA(K,LPIVCO) = SWAP
  100       CONTINUE
            INFO(2) = -INFO(2)
         END IF
         DO 110 K = 1,KFR - 1
            FA(K,LFR) = -FA(K,LFR)
  110    CONTINUE
         FA(KFR,LFR) = -PIVOT
         KR = KR - 1
         IF (NSIZE.EQ.1) THEN
            KRTEMP = KR
            KR = 0
         END IF
         PIVINV = ONE/FA(KFR,LFR)
         DO 120 L = 1,KR
            LCOL = KPIV(L)
            FA(KFR,LCOL) = -FA(KFR,LCOL)*PIVINV
            CALL DAXPY(KFR-1,FA(KFR,LCOL),FA(1,LFR),1,FA(1,LCOL),1)
  120    CONTINUE
         IF (KFR.NE.KPIVRO) THEN
            DO 130 L = 1,NRHS
               SWAP = FRHS(KFR,L)
               FRHS(KFR,L) = FRHS(KPIVRO,L)
               FRHS(KPIVRO,L) = SWAP
  130       CONTINUE
         END IF
         DET = DET + LOG(ABS(PIVOT))
         KFR = KFR - 1
         LFR = LFR - 1
  140 CONTINUE
  150 CHFRNT = LFRNT - LFR
      DO 200 L = 1,KR
         LCOL = KPIV(L)
         IF (LCOL.GE.LFR-KR+1) GO TO 200
         DO 190 LK = LFR,LFR - KR + 1,-1
            IF (KPVLNK(LK).LT.0) THEN
               KPVLNK(LCOL) = KPVLNK(LK)
               KPIV(L) = LK
               KPVLNK(LK) = L
               IF (LPIVCX.EQ.LCOL) LPIVCX = LK
               LFRO = LHED(LK)
               LPOS = -LAST(LFRO)
               LDEST(LPOS) = LCOL
               LTEMP = LHED(LK)
               LHED(LK) = LHED(LCOL)
               LHED(LCOL) = LTEMP
               DO 160 K = 1,KFRNT
                  SWAP = FA(K,LK)
                  FA(K,LK) = FA(K,LCOL)
                  FA(K,LCOL) = SWAP
  160          CONTINUE
               INFO(2) = -INFO(2)
               IF (ABS(ICNTL(7)).EQ.1001 .AND. KPIVRO.EQ.LPIVCO) THEN
                  IF (NMAXE.GT.1) THEN
                     KFRO = KHED(LK)
                     KPOS = -LAST(KFRO)
                     KDEST(KPOS) = LCOL
                  END IF
                  LTEMP = KHED(LK)
                  KHED(LK) = KHED(LCOL)
                  KHED(LCOL) = LTEMP
                  DO 170 K = 1,LFRNT
                     SWAP = FA(LK,K)
                     FA(LK,K) = FA(LCOL,K)
                     FA(LCOL,K) = SWAP
  170             CONTINUE
                  INFO(2) = -INFO(2)
                  DO 180 K = 1,NRHS
                     SWAP = FRHS(LK,K)
                     FRHS(LK,K) = FRHS(LCOL,K)
                     FRHS(LCOL,K) = SWAP
  180             CONTINUE
               END IF
               GO TO 200
            END IF
  190    CONTINUE
  200 CONTINUE
      IF ((CHFRNT/2)*2.NE.CHFRNT) THEN
         IF ((KFR/2)*2.NE.KFR) INFO(2) = -INFO(2)
         IF ((LFR/2)*2.NE.LFR) INFO(2) = -INFO(2)
      END IF
      KPRE = 0
      LPRE = 0
      IF (.NOT.LSTAT) GO TO 370
      DO 240 L = LFR - KR,1,-1
          IF (L.EQ.LPRE) GO TO 250
          DO 210 LK = 1,CHFRNT
             KPIVRO = KFR + LK
             IF (FA(KPIVRO,L).NE.ZERO) GO TO 240
  210     CONTINUE
          LPRE = LPRE + 1
            L1 = LPRE
            DO 220 LL = L1,L - 1
               DO 215 LK = 1,CHFRNT
                  KPIVRO = KFR + LK
                  IF (FA(KPIVRO,LL).NE.ZERO) GO TO 225
  215          CONTINUE
               LPRE = LPRE + 1
 220       CONTINUE
            GO TO 250
  225       CONTINUE
             J = LHED(L)
             JJ = LHED(LL)
            IF (KPVLNK(L).LT.0) THEN
               IF (KPVLNK(LL).LT.0) THEN
                  I = KPVLNK(L)
                  KPVLNK(L) = KPVLNK(LL)
                  KPVLNK(LL) = I
                  LDEST(-LAST(J)) = LL
                  LDEST(-LAST(JJ)) = L
               ELSE
                  I = KPVLNK(LL)
                  KPVLNK(LL) = KPVLNK(L)
                  KPVLNK(L) = I
                  KPIV(I) = L
                  LDEST(-LAST(J)) = LL
               END IF
            ELSE
               IF (KPVLNK(LL).LT.0) THEN
                  I = KPVLNK(L)
                  KPVLNK(L) = KPVLNK(LL)
                  KPVLNK(LL) = I
                  KPIV(I) = LL
                  LDEST(-LAST(JJ)) = L
               ELSE
                  I = KPVLNK(L)
                  LCOL = KPVLNK(LL)
                  KPVLNK(L) = LCOL
                  KPIV(LCOL) = L
                  KPVLNK(LL) = I
                  KPIV(I) = LL
               END IF
            END IF
            LHED(L) = JJ
            LHED(LL) = J
            INFO(2) = -INFO(2)
            DO 230 K = 1,KFRNT
               SWAP = FA(K,L)
               FA(K,L) = FA(K,LL)
               FA(K,LL) = SWAP
  230       CONTINUE
  240 CONTINUE
  250 CONTINUE
      DO 310 K = KFR - KR,1,-1
         IF (K.EQ.KPRE) GO TO 320
         DO 260 LK = 1,CHFRNT
             LPIVCO = LFR + LK
             IF (FA(K,LPIVCO).NE.ZERO) GO TO 310
  260    CONTINUE
         KPRE = KPRE + 1
            L1 = KPRE
            DO 280 KK = L1,K - 1
               DO 270 LK = 1,CHFRNT
                  LPIVCO = LFR + LK
                  IF (FA(KK,LPIVCO).NE.ZERO) GO TO 290
  270          CONTINUE
               KPRE = KPRE + 1
  280       CONTINUE
            GO TO 320
  290       CONTINUE
            J = KHED(K)
            JJ = KHED(KK)
            IF (NMAXE.GT.1) THEN
               IF (LAST(J).LT.0) KDEST(-LAST(J)) = KK
               IF (LAST(JJ).LT.0) KDEST(-LAST(JJ)) = K
            END IF
            KHED(K) = JJ
            KHED(KK) = J
            INFO(2) = -INFO(2)
            DO 300 L = 1,LFRNT
               SWAP = FA(K,L)
               FA(K,L) = FA(KK,L)
               FA(KK,L) = SWAP
  300       CONTINUE
            DO 305 LK = 1,NRHS
               SWAP = FRHS(K,LK)
               FRHS(K,LK) = FRHS(KK,LK)
               FRHS(KK,LK) = SWAP
  305       CONTINUE
  310 CONTINUE
  320 CONTINUE
  370 CONTINUE
      IF (CHFRNT.GT.0) THEN
         IF (CHFRNT.EQ.1) THEN
            PIVINV = -ONE/FA(KFRNT,LFRNT)
            DO 380 L = 1,NRHS
               FRHS(KFRNT,L) = FRHS(KFRNT,L)*PIVINV
  380       CONTINUE
            IF (KFR.NE.KPRE .AND. NRHS.NE.0) CALL DGER(KFR-KPRE,NRHS,
     +          ONE,FA(KPRE+1,LFRNT),1,FRHS(KFRNT,1),MFRONT,
     +          FRHS(KPRE+1,1),MFRONT)
            IF (LFR.NE.KR+LPRE) THEN
               DO 390 K = 1,LFR - KR - LPRE
                  FA(KFRNT,LPRE+K) = FA(KFRNT,LPRE+K)*PIVINV
  390          CONTINUE
               IF (KFR.NE.KPRE) CALL DGER(KFR-KPRE,LFR-KR-LPRE,ONE,
     +                                    FA(KPRE+1,LFR+1),1,
     +                                    FA(KFR+1,LPRE+1),MFRONT,
     +                                    FA(KPRE+1,LPRE+1),MFRONT)
            END IF
         ELSE
            CALL DTRSM('L','U','N','N',CHFRNT,NRHS,-ONE,FA(KFR+1,LFR+1),
     +                 MFRONT,FRHS(KFR+1,1),MFRONT)
            IF (KFR.NE.KPRE .AND. NRHS.NE.0) CALL DGEMM('N','N',
     +          KFR-KPRE,NRHS,CHFRNT,ONE,FA(KPRE+1,LFR+1),MFRONT,
     +          FRHS(KFR+1,1),MFRONT,ONE,FRHS(1+KPRE,1),MFRONT)
            IF (LFR.NE.KR+LPRE) THEN
               CALL DTRSM('L','U','N','N',CHFRNT,LFR-KR-LPRE,-ONE,
     +                    FA(KFR+1,LFR+1),MFRONT,FA(KFR+1,LPRE+1),
     +                    MFRONT)
               IF (KFR.NE.KPRE) CALL DGEMM('N','N',KFR-KPRE,LFR-KR-LPRE,
     +                                     CHFRNT,ONE,FA(KPRE+1,LFR+1),
     +                                     MFRONT,FA(KFR+1,LPRE+1),
     +                                     MFRONT,ONE,FA(KPRE+1,LPRE+1),
     +                                     MFRONT)
            END IF
         END IF
         CALL MA42HD(1,KFRNT-KPRE,LFRNT-LPRE,CHFRNT,IBUFR,LIBUFR,
     +               KHED(1+KPRE),LHED(1+LPRE),MFRONT-KPRE,NFRONT-LPRE,
     +               LP,IFILE(3),IREC(3),ISIZE(3),MKEY(3),NUMBLK(3),
     +               IELL,NELL,INFO)
         IF (INFO(1).LT.0) GO TO 420
         CALL MA42GD(1,1,FA(KPRE+1,LPRE+1),FRHS(1+KPRE,1),NRHS,
     +               KFRNT-KPRE,LFRNT-LPRE,CHFRNT,BUFRU,LUB,MFRONT,
     +               LP,IFILE(1),IREC(1),ISIZE(1),MKEY(1),
     +               NUMBLK(1),IELL,NELL,INFO)
         IF (INFO(1).LT.0) GO TO 420
         CALL MA42GD(2,1,FA(KPRE+1,LPRE+1),FRHS(1+KPRE,1),NRHS,
     +               KFRNT-KPRE,LFRNT-LPRE,CHFRNT,BUFRL,LLB,MFRONT,
     +               LP,IFILE(2),IREC(2),ISIZE(2),MKEY(2),
     +               NUMBLK(2),IELL,NELL,INFO)
         IF (INFO(1).LT.0) GO TO 420
         DO 400 J = 1,CHFRNT
            OPS = OPS + DBLE((KFRNT-J-KPRE)* (LFRNT-J-LPRE)*2)
  400    CONTINUE
         LFRNT = LFR
         KFRNT = KFR
      END IF
      IF (NSIZE.EQ.1) KR = KRTEMP
      GO TO 420
  410 INFO(1) = -14
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9010) INFO(1),IELL
         WRITE (LP,FMT=9000)
      END IF
  420 RETURN
 9000 FORMAT (7X,'Matrix found to be singular')
 9010 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9020 FORMAT (' ***** Warning from MA42B/BD *****  INFO(1) = ',I3,/7X,
     +       'after input of elt/equ ',I8)
      END
C**********************************************************************
      SUBROUTINE MA42OD(IELL,AVAR,LRHS,NMAXE,RHS,NRHS,IVAR,NDF,LAST,
     +                  NVAR,MFRONT,NFRONT,LHED,KHED,KR,KPIV,KPVLNK,
     +                  LDEST,KDEST,ISTATC,LFREE,FRHS,BUFRL,LLB,BUFRU,
     +                  LUB,IBUFR,LIBUFR,MVAR,KFRNT,LFRNT,FA,IFILE,IREC,
     +                  ISIZE,MKEY,NUMBLK,DET,OPS,NELL,CNTL,ICNTL,INFO)
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
      DOUBLE PRECISION DET,OPS
      INTEGER IELL,ISTATC,KFRNT,KR,LFREE,LFRNT,LIBUFR,LLB,LRHS,LUB,
     +        MFRONT,MVAR,NDF,NELL,NFRONT,NMAXE,NRHS,NVAR
      DOUBLE PRECISION AVAR(NMAXE,NVAR),BUFRL(LLB),BUFRU(LUB),CNTL(2),
     +                 FA(MFRONT,NFRONT),FRHS(MFRONT,LRHS),
     +                 RHS(NMAXE,LRHS)
      INTEGER IBUFR(LIBUFR),ICNTL(8),IFILE(3),INFO(23),IREC(3),ISIZE(3),
     +        IVAR(NVAR),KDEST(MFRONT),KHED(MFRONT),KPIV(NFRONT),
     +        KPVLNK(NFRONT),LAST(NDF),LDEST(NFRONT),LHED(NFRONT),
     +        MKEY(3),NUMBLK(3)
      DOUBLE PRECISION PIVINV,PIVOT,RMAX,SMAX,SWAP
      INTEGER I,IMAX,IPIV,IPOS,IRHS,ISWAP,J,J1,JCOL,JJ,JSTATC,JVAR,K,
     +        KDST,KFR,KFRNEW,KFROLD,KPOS,KR1,L,LC1,LDST,LFR,LFRNEW,
     +        LFROLD,LK,LP,MFR,MOVE,MP,NFREE,NUMPIV,NUVAR,NUVRP1
      LOGICAL LDET
      INTEGER IDAMAX
      EXTERNAL IDAMAX
      EXTERNAL DGER,MA42GD,MA42HD
      INTRINSIC ABS,LOG,MAX,MIN,DBLE
      IF (ISTATC.EQ.0) GO TO 490
      LP = ICNTL(1)
      MP = ICNTL(2)
      NUVAR = NVAR - ISTATC
      JCOL = NVAR
      J1 = 1
      DO 60 MOVE = 1,ISTATC
         MFR = IVAR(JCOL)
         IF (LAST(MFR).NE.IELL) THEN
            LDET = LAST(MFR) .GT. 0
            DO 10 J = J1,JCOL - 1
               MFR = IVAR(J)
               IF (LAST(MFR).EQ.IELL) GO TO 20
   10       CONTINUE
   20       ISWAP = IVAR(J)
            IVAR(J) = IVAR(JCOL)
            IVAR(JCOL) = ISWAP
            J1 = J + 1
            IF (LDET) INFO(2) = -INFO(2)
            DO 30 I = 1,MVAR
               SWAP = AVAR(I,J)
               AVAR(I,J) = AVAR(I,JCOL)
               AVAR(I,JCOL) = SWAP
   30       CONTINUE
            IF (NMAXE.GT.1) THEN
               IF (LDET) INFO(2) = -INFO(2)
               DO 40 I = 1,NVAR
                  SWAP = AVAR(J,I)
                  AVAR(J,I) = AVAR(JCOL,I)
                  AVAR(JCOL,I) = SWAP
   40          CONTINUE
               DO 50 IRHS = 1,NRHS
                  SWAP = RHS(J,IRHS)
                  RHS(J,IRHS) = RHS(JCOL,IRHS)
                  RHS(JCOL,IRHS) = SWAP
   50          CONTINUE
            END IF
         END IF
         JCOL = JCOL - 1
   60 CONTINUE
      DO 70 LK = 1,NUVAR
         MFR = IVAR(LK)
         IF (LAST(MFR).GE.0) THEN
            LFRNT = LFRNT + 1
            NFREE = LDEST(LFREE)
            LDEST(LFREE) = LFRNT
            LHED(LFRNT) = MFR
            IF (NMAXE.GT.1) THEN
               KFRNT = KFRNT + 1
               KHED(KFRNT) = MFR
               KDEST(LFREE) = KFRNT
            END IF
            KPVLNK(LFRNT) = -LAST(MFR)
            LAST(MFR) = -LFREE
            LFREE = NFREE
         END IF
   70 CONTINUE
      IF (NMAXE.EQ.1) THEN
         MFR = IVAR(NVAR)
         KFRNT = KFRNT + 1
         LFRNT = LFRNT + 1
         KHED(KFRNT) = IELL
         LHED(LFRNT) = MFR
         PIVOT = AVAR(1,NVAR)
         IF (ABS(PIVOT).LE.CNTL(1)) THEN
            NUMPIV = 0
            INFO(2) = 0
            DET = ZERO
            IF (ICNTL(8).EQ.0) GO TO 460
            IF (INFO(1).EQ.0) THEN
               INFO(1) = 1
               IF (MP.GT.0) THEN
                  WRITE (MP,FMT=9020) INFO(1),IELL
                  WRITE (MP,FMT=9010)
               END IF
            END IF
            GO TO 300
         END IF
         DET = DET + LOG(ABS(PIVOT))
         IF (PIVOT.LT.ZERO) INFO(2) = -INFO(2)
         AVAR(1,NVAR) = -PIVOT
         PIVINV = ONE/PIVOT
         DO 80 IRHS = 1,NRHS
            RHS(1,IRHS) = RHS(1,IRHS)*PIVINV
   80    CONTINUE
         DO 90 J = 1,NUVAR
            AVAR(1,J) = AVAR(1,J)*PIVINV
   90    CONTINUE
      ELSE
         NUVRP1 = NUVAR + 1
         DO 100 I = NUVRP1,NVAR
            MFR = IVAR(I)
            LFRNT = LFRNT + 1
            LHED(LFRNT) = MFR
            KFRNT = KFRNT + 1
            KHED(KFRNT) = MFR
  100    CONTINUE
         JCOL = NVAR
         IF (ABS(ICNTL(7)).EQ.1001) GO TO 200
         DO 190 JSTATC = 1,ISTATC
            DO 110 J = JCOL,NUVRP1,-1
               INFO(13) = INFO(13) + 1
               INFO(15) = INFO(15) + JCOL
               IMAX = IDAMAX(JCOL,AVAR(1,J),1)
               RMAX = ABS(AVAR(IMAX,J))
               IF (RMAX.LE.CNTL(1)) THEN
                  INFO(2) = 0
                  DET = ZERO
                  IF (ICNTL(8).EQ.0) THEN
                     NUMPIV = JSTATC - 1
                     GO TO 460
                  ELSE IF (INFO(1).EQ.0) THEN
                     INFO(1) = 1
                     IF (MP.GT.0) THEN
                        WRITE (MP,FMT=9020) INFO(1),IELL
                        WRITE (MP,FMT=9010)
                     END IF
                  END IF
                  IF (J.EQ.NUVRP1) THEN
                     NUMPIV = JSTATC - 1
                     GO TO 300
                  END IF
                  GO TO 110
               END IF
               INFO(14) = INFO(14) + 1
               IF (ABS(AVAR(JCOL,J)).GE.CNTL(2)*RMAX) GO TO 150
               IF (IMAX.GE.NUVRP1) GO TO 120
               INFO(14) = INFO(14) + JCOL - NUVAR
               IMAX = IDAMAX(JCOL-NUVAR,AVAR(NUVRP1,J),1) + NUVAR
               SMAX = ABS(AVAR(IMAX,J))
               IF (SMAX.GE.CNTL(2)*RMAX) GO TO 120
               IF (J.EQ.NUVRP1) THEN
                  NUMPIV = JSTATC - 1
                  GO TO 300
               END IF
  110       CONTINUE
  120       DO 130 JJ = 1,NVAR
               SWAP = AVAR(JCOL,JJ)
               AVAR(JCOL,JJ) = AVAR(IMAX,JJ)
               AVAR(IMAX,JJ) = SWAP
  130       CONTINUE
            INFO(2) = -INFO(2)
            KFRNEW = KFRNT + JCOL - NVAR
            KFROLD = KFRNT + IMAX - NVAR
            ISWAP = KHED(KFRNEW)
            KHED(KFRNEW) = KHED(KFROLD)
            KHED(KFROLD) = ISWAP
            DO 140 IRHS = 1,NRHS
               SWAP = RHS(JCOL,IRHS)
               RHS(JCOL,IRHS) = RHS(IMAX,IRHS)
               RHS(IMAX,IRHS) = SWAP
  140       CONTINUE
  150       IF (J.NE.JCOL) THEN
               DO 160 I = 1,NVAR
                  SWAP = AVAR(I,J)
                  AVAR(I,J) = AVAR(I,JCOL)
                  AVAR(I,JCOL) = SWAP
  160          CONTINUE
               INFO(2) = -INFO(2)
               LFRNEW = LFRNT + JCOL - NVAR
               LFROLD = LFRNT + J - NVAR
               ISWAP = LHED(LFRNEW)
               LHED(LFRNEW) = LHED(LFROLD)
               LHED(LFROLD) = ISWAP
            END IF
            PIVOT = AVAR(JCOL,JCOL)
            DET = DET + LOG(ABS(PIVOT))
            IF (PIVOT.LT.ZERO) INFO(2) = -INFO(2)
            AVAR(JCOL,JCOL) = -PIVOT
            PIVINV = ONE/PIVOT
            DO 170 I = 1,JCOL - 1
               AVAR(I,JCOL) = -AVAR(I,JCOL)
               AVAR(JCOL,I) = AVAR(JCOL,I)*PIVINV
  170       CONTINUE
            DO 180 IRHS = 1,NRHS
               RHS(JCOL,IRHS) = RHS(JCOL,IRHS)*PIVINV
  180       CONTINUE
            IF (JCOL.NE.1) CALL DGER(JCOL-1,JCOL-1,ONE,AVAR(1,JCOL),1,
     +                               AVAR(JCOL,1),NMAXE,AVAR(1,1),NMAXE)
            IF (JCOL.NE.1 .AND. NRHS.NE.0) CALL DGER(JCOL-1,NRHS,ONE,
     +          AVAR(1,JCOL),1,RHS(JCOL,1),NMAXE,RHS(1,1),NMAXE)
            JCOL = JCOL - 1
  190    CONTINUE
         GO TO 290
  200    CONTINUE
         DO 280 JSTATC = 1,ISTATC
            DO 210 J = JCOL,NUVRP1,-1
               INFO(13) = INFO(13) + 1
               INFO(15) = INFO(15) + JCOL
               IMAX = IDAMAX(JCOL,AVAR(1,J),1)
               RMAX = ABS(AVAR(IMAX,J))
               IF (RMAX.LE.CNTL(1)) THEN
                  INFO(2) = 0
                  DET = ZERO
                  IF (ICNTL(8).EQ.0) THEN
                     NUMPIV = JSTATC - 1
                     GO TO 460
                  ELSE IF (INFO(1).EQ.0) THEN
                     INFO(1) = 1
                     IF (MP.GT.0) THEN
                        WRITE (MP,FMT=9020) INFO(1),IELL
                        WRITE (MP,FMT=9010)
                     END IF
                  END IF
                  IF (J.EQ.NUVRP1) THEN
                     NUMPIV = JSTATC - 1
                     GO TO 300
                  END IF
                  GO TO 210
               END IF
               INFO(14) = INFO(14) + 1
               IF (ABS(AVAR(J,J)).GE.CNTL(2)*RMAX) GO TO 220
               IF (J.EQ.NUVRP1) THEN
                  NUMPIV = JSTATC - 1
                  GO TO 300
               END IF
  210       CONTINUE
  220       IF (J.NE.JCOL) THEN
               DO 230 JJ = 1,NVAR
                  SWAP = AVAR(JCOL,JJ)
                  AVAR(JCOL,JJ) = AVAR(J,JJ)
                  AVAR(J,JJ) = SWAP
  230          CONTINUE
               KFRNEW = KFRNT + JCOL - NVAR
               KFROLD = KFRNT + J - NVAR
               ISWAP = KHED(KFRNEW)
               KHED(KFRNEW) = KHED(KFROLD)
               KHED(KFROLD) = ISWAP
               DO 240 IRHS = 1,NRHS
                  SWAP = RHS(JCOL,IRHS)
                  RHS(JCOL,IRHS) = RHS(J,IRHS)
                  RHS(J,IRHS) = SWAP
  240          CONTINUE
               DO 250 I = 1,NVAR
                  SWAP = AVAR(I,J)
                  AVAR(I,J) = AVAR(I,JCOL)
                  AVAR(I,JCOL) = SWAP
  250          CONTINUE
               LFRNEW = LFRNT + JCOL - NVAR
               LFROLD = LFRNT + J - NVAR
               ISWAP = LHED(LFRNEW)
               LHED(LFRNEW) = LHED(LFROLD)
               LHED(LFROLD) = ISWAP
            END IF
            PIVOT = AVAR(JCOL,JCOL)
            DET = DET + LOG(ABS(PIVOT))
            IF (PIVOT.LT.ZERO) INFO(2) = -INFO(2)
            AVAR(JCOL,JCOL) = -PIVOT
            PIVINV = ONE/PIVOT
            DO 260 I = 1,JCOL - 1
               AVAR(I,JCOL) = -AVAR(I,JCOL)
               AVAR(JCOL,I) = AVAR(JCOL,I)*PIVINV
  260       CONTINUE
            DO 270 IRHS = 1,NRHS
               RHS(JCOL,IRHS) = RHS(JCOL,IRHS)*PIVINV
  270       CONTINUE
            IF (JCOL.NE.1) CALL DGER(JCOL-1,JCOL-1,ONE,AVAR(1,JCOL),1,
     +                               AVAR(JCOL,1),NMAXE,AVAR(1,1),NMAXE)
            IF (JCOL.NE.1 .AND. NRHS.NE.0) CALL DGER(JCOL-1,NRHS,ONE,
     +          AVAR(1,JCOL),1,RHS(JCOL,1),NMAXE,RHS(1,1),NMAXE)
            JCOL = JCOL - 1
  280    CONTINUE
  290    CONTINUE
      END IF
      NUMPIV = ISTATC
  300 INFO(17) = INFO(17) + NUMPIV
      INFO(18) = INFO(18) + ISTATC
      KR1 = INFO(8) + 1
      LC1 = INFO(9) + 1
      IF (KR1.LE.KFRNT) THEN
         INFO(8) = KFRNT
         LFR = MAX(INFO(9),LFRNT)
         DO 330 K = KR1,KFRNT
            DO 310 L = 1,LFR
               FA(K,L) = ZERO
  310       CONTINUE
            DO 320 J = 1,NRHS
               FRHS(K,J) = ZERO
  320       CONTINUE
  330    CONTINUE
      END IF
      IF (LC1.LE.LFRNT) INFO(9) = LFRNT
      DO 350 L = LC1,LFRNT
         DO 340 K = 1,MIN(MFRONT,KR1)
            FA(K,L) = ZERO
  340    CONTINUE
  350 CONTINUE
      KFR = KFRNT
      LFR = LFRNT
      IF (NMAXE.EQ.1 .AND. NUMPIV.EQ.1) GO TO 455
      DO 360 LK = 1,NUVAR
         MFR = IVAR(LK)
         IPOS = -LAST(MFR)
         LDST = LDEST(IPOS)
         IVAR(LK) = LDST
  360 CONTINUE
      IF (NMAXE.EQ.1) THEN
         IF (NUMPIV.EQ.0) THEN
          DO 370 J = 1,NVAR - 1
            LDST = IVAR(J)
            FA(KFR,LDST) = AVAR(1,J)
  370     CONTINUE
          FA(KFR,LFR) = AVAR(1,NVAR)
          DO 380 J = 1,NRHS
            FRHS(KFR,J) = RHS(1,J)
  380     CONTINUE
         END IF
      ELSE
         JCOL = NVAR
         DO 440 IPIV = 1,ISTATC
            DO 390 J = 1,NUVAR
               LDST = IVAR(J)
               FA(KFR,LDST) = AVAR(JCOL,J)
  390       CONTINUE
            DO 400 J = NUVAR + 1,JCOL
               FA(KFR,LFR+J-JCOL) = AVAR(JCOL,J)
  400       CONTINUE
            DO 410 J = 1,NUVAR
               LDST = IVAR(J)
               JVAR = LHED(LDST)
               KPOS = -LAST(JVAR)
               KDST = KDEST(KPOS)
               FA(KDST,LFR) = AVAR(J,JCOL)
  410       CONTINUE
            DO 420 J = NUVAR + 1,JCOL
               FA(KFR+J-JCOL,LFR) = AVAR(J,JCOL)
  420       CONTINUE
            DO 430 J = 1,NRHS
               FRHS(KFR,J) = RHS(JCOL,J)
  430       CONTINUE
            KFR = KFR - 1
            LFR = LFR - 1
            JCOL = JCOL - 1
  440    CONTINUE
      END IF
      DO 450 LK = 1,NUVAR
         LDST = ABS(IVAR(LK))
         JVAR = LHED(LDST)
         IVAR(LK) = JVAR
  450 CONTINUE
  455 CONTINUE
      IF ((NUMPIV/2)*2.NE.NUMPIV) THEN
         KFR = KFRNT - NUMPIV
         LFR = LFRNT - NUMPIV
         IF ((KFR/2)*2.NE.KFR) INFO(2) = -INFO(2)
         IF ((LFR/2)*2.NE.LFR) INFO(2) = -INFO(2)
      END IF
      IF (NUMPIV.GT.0) THEN
          IF (NMAXE.EQ.1) THEN
             CALL MA42HD(1,1,NVAR,NUMPIV,IBUFR,LIBUFR,
     +                   KHED(KFRNT),IVAR,
     +                   1,NVAR,LP,IFILE(3),IREC(3),ISIZE(3),MKEY(3),
     +                   NUMBLK(3),IELL,NELL,INFO)
             IF (INFO(1).EQ.-26) GO TO 490
             CALL MA42GD(1,1,AVAR,RHS,NRHS,MVAR,NVAR,NUMPIV,BUFRU,
     +               LUB,NMAXE,LP,IFILE(1),IREC(1),ISIZE(1),
     +               MKEY(1),NUMBLK(1),IELL,NELL,INFO)
             IF (INFO(1).EQ.-26) GO TO 490
             CALL MA42GD(2,1,AVAR,RHS,NRHS,MVAR,NVAR,NUMPIV,BUFRL,
     +               LLB,NMAXE,LP,IFILE(2),IREC(2),ISIZE(2),
     +               MKEY(2),NUMBLK(2),IELL,NELL,INFO)
             IF (INFO(1).EQ.-26) GO TO 490
          ELSE
             CALL MA42HD(1,KFRNT,LFRNT,NUMPIV,IBUFR,LIBUFR,KHED,LHED,
     +               MFRONT,NFRONT,LP,IFILE(3),IREC(3),ISIZE(3),MKEY(3),
     +               NUMBLK(3),IELL,NELL,INFO)
             IF (INFO(1).EQ.-26) GO TO 490
             CALL MA42GD(1,1,FA,FRHS,NRHS,KFRNT,LFRNT,NUMPIV,BUFRU,
     +               LUB,MFRONT,LP,IFILE(1),IREC(1),ISIZE(1),
     +               MKEY(1),NUMBLK(1),IELL,NELL,INFO)
             IF (INFO(1).EQ.-26) GO TO 490
             CALL MA42GD(2,1,FA,FRHS,NRHS,KFRNT,LFRNT,NUMPIV,BUFRL,
     +               LLB,MFRONT,LP,IFILE(2),IREC(2),ISIZE(2),
     +               MKEY(2),NUMBLK(2),IELL,NELL,INFO)
             IF (INFO(1).EQ.-26) GO TO 490
         END IF
         DO 456 J = 1,NUMPIV
            OPS = OPS + DBLE((NVAR-J)* (MVAR-J)*2)
  456    CONTINUE
         KFRNT = KFRNT - NUMPIV
         LFRNT = LFRNT - NUMPIV
      END IF
      GO TO 470
  460 INFO(1) = -14
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1),IELL
         WRITE (LP,FMT=9010)
      END IF
      NUMPIV = 0
  470 DO 480 J = 1,ISTATC - NUMPIV
         KR = KR + 1
         LDST = LFRNT - J + 1
         KPIV(KR) = LDST
         KPVLNK(LDST) = KR
  480 CONTINUE
  490 RETURN
 9000 FORMAT (' ***** Error return from MA42B/BD *****  INFO(1) = ',I3,
     +       /' ***** Error discovered after input of elt/eqn *****',I8)
 9010 FORMAT (7X,'Matrix found to be singular')
 9020 FORMAT (' ***** Warning from MA42B/BD *****  INFO(1) = ',I3,/7X,
     +       'after input of elt/equ ',I8)
      END



C *******************************************************************
C COPYRIGHT (c) 1998 Council for the Central Laboratory
*                    of the Research Councils
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
C Licence, see http://hsl.rl.ac.uk/acuk/cou.html
C
C Please note that for a UK ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between AEA
C    Technology plc and the Licensee on suitable terms and conditions,
C    which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither CCLRC nor AEA Technology plc shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C DATE 16 Jan 1998
      SUBROUTINE MC60AD(N,LIRN,IRN,ICPTR,ICNTL,IW,INFO)
      INTEGER N
      INTEGER LIRN
      INTEGER IRN(LIRN)
      INTEGER ICPTR(N+1)
      INTEGER ICNTL(2)
      INTEGER IW(N)
      INTEGER INFO(4)
      INTEGER CKP1,I,I1,I2,II,IOUT,IPOS,IREP,J,KZ,LP,NDIAG,NEWTAU
      LP = ICNTL(2)
      DO 5 J = 1,4
         INFO(J) = 0
    5 CONTINUE
      IF (N.LT.1) THEN
          INFO(1) = -1
          IF (LP.GT.0) WRITE (LP,'(/,A,I3/,A,I6)')
     *         ' MC60A/AD error: INFO(1) =', INFO(1), ' N =',N
          RETURN
      END IF
      IF (LIRN.LT.ICPTR(N+1)-1) THEN
          INFO(1) = -1
          IF (LP.GT.0) WRITE (LP,'(/,A,I3/,A)')
     *         ' MC60A/AD error:  INFO(1) =', INFO(1),
     *         ' LIRN is less than ICPTR(N+1)-1'
          RETURN
      END IF
      DO 10 I = 1,N
        IW(I) = 0
   10 CONTINUE
      IOUT = 0
      IREP = 0
      KZ = 0
      IF (ICNTL(1).EQ.1) THEN
         I1 = ICPTR(1)
         ICPTR(1) = 1
         DO 16 J = 1,N
            DO 15 II = I1,ICPTR(J+1) - 1
               I = IRN(II)
               IF (I.GT.N .OR. I.LT.J) THEN
                  IOUT = IOUT + 1
               ELSE IF (IW(I).EQ.J) THEN
                  IREP = IREP + 1
               ELSE
                  KZ = KZ + 1
                  IRN(KZ) = I
                  IW(I) = J
               END IF
   15       CONTINUE
            I1 = ICPTR(J+1)
            ICPTR(J+1) = KZ + 1
   16    CONTINUE
         IF (IOUT.GT.0) THEN
             INFO(1) = 1
             IF (LP.GT.0) WRITE (LP,'(/,A,I6,A)')
     *         ' MC60A/AD warning:',IOUT,' out-of-range entries ignored'
         END IF
         IF (IREP.GT.0) THEN
             INFO(1) = 1
             IF (LP.GT.0) WRITE (LP,'(/,A,I6,A)')
     *         ' MC60A/AD warning:',IREP,' duplicated entries ignored'
         END IF
         INFO(2) = IOUT
         INFO(3) = IREP
      ELSE
         I1 = ICPTR(1)
         DO 26 J = 1,N
            DO 25 II = I1,ICPTR(J+1) - 1
               I = IRN(II)
               IF (I.GT.N .OR. I.LT.J) THEN
                  IOUT = IOUT + 1
               ELSE IF (IW(I).EQ.J) THEN
                  IREP = IREP + 1
               ELSE
                  KZ = KZ + 1
                  IW(I) = J
               END IF
   25       CONTINUE
            I1 = ICPTR(J+1)
   26    CONTINUE
         IF (IOUT.GT.0 .OR. IREP.GT.0) THEN
            INFO(1) = -3
            IF (LP.GT.0) THEN
              WRITE (LP,'(/,A,I3)')
     *            ' MC60A/AD error:  INFO(1) =', INFO(1)
            IF (IOUT.GT.0) WRITE (LP,'(A,I6,A)')
     *         ' There are ',IOUT,' out-of-range entries'
            IF (IREP.GT.0) WRITE (LP,'(A,I6,A)')
     *         ' There are ',IREP,' duplicated entries'
            END IF
            INFO(2) = IOUT
            INFO(3) = IREP
            RETURN
         END IF
      END IF
      DO 30 J = 1,N
        IW(J) = 0
   30 CONTINUE
      NDIAG = 0
      DO 40 J = 1,N
        I1 = ICPTR(J)
        I2 = ICPTR(J+1) - 1
        IW(J) = IW(J) + I2 - I1 + 1
        DO 35 II = I1,I2
          I = IRN(II)
          IF (I.NE.J) THEN
            IW(I) = IW(I) + 1
          ELSE
            NDIAG = NDIAG + 1
          END IF
   35   CONTINUE
   40 CONTINUE
      NEWTAU = 2*KZ - NDIAG
      INFO(4) = NEWTAU
      IF (NEWTAU.GT.LIRN) THEN
          INFO(1) = -2
          IF (LP.GT.0) WRITE (LP,'(/,A)')
     *         ' MC60A/AD error: LIRN is too small'
          RETURN
      END IF
      I1 = KZ + 1
      CKP1 = NEWTAU + 1
      DO 60 J = N,1,-1
        I2 = I1 - 1
        I1 = ICPTR(J)
        IPOS = CKP1
          DO 50 II = I2,I1,-1
            IPOS = IPOS - 1
            IRN(IPOS) = IRN(II)
   50     CONTINUE
        ICPTR(J) = IPOS
        CKP1 = CKP1 - IW(J)
        IW(J) = I2 - I1 + 1
   60 CONTINUE
      DO 80 J = N,1,-1
        I1 = ICPTR(J)
        I2 = ICPTR(J) + IW(J) - 1
        IF(I1.LE.I2) THEN
          DO 70 II = I1,I2
            I = IRN(II)
            IF(I.EQ.J) GO TO 70
            ICPTR(I) = ICPTR(I) - 1
            IRN(ICPTR(I)) = J
   70     CONTINUE
        END IF
   80 CONTINUE
      ICPTR(N+1) = NEWTAU + 1
      END
      SUBROUTINE MC60BD(N,LIRN,IRN,ICPTR,NSUP,SVAR,VARS,IW)
      INTEGER N
      INTEGER LIRN
      INTEGER IRN(LIRN)
      INTEGER ICPTR(N+1)
      INTEGER NSUP
      INTEGER SVAR(N)
      INTEGER VARS(N)
      INTEGER IW(2*N)
      INTEGER FLAG
      EXTERNAL MC60OD,MC60PD
      FLAG = N+1
      CALL MC60OD(N,N,LIRN,IRN,ICPTR,SVAR,NSUP,IW,VARS,IW(FLAG))
      CALL MC60PD(N,NSUP,LIRN,IRN,ICPTR,SVAR,VARS,IW,IW(FLAG))
      END
      SUBROUTINE MC60CD(N,NSUP,LIRN,IRN,ICPTR,VARS,JCNTL,
     +                  PERMSV,WEIGHT,PAIR,INFO,IW,W)
      INTEGER N,NSUP,LIRN,IRN(LIRN),ICPTR(NSUP+1),VARS(NSUP),JCNTL(2)
      INTEGER PERMSV(NSUP),PAIR(2,NSUP/2),INFO(4),IW(3*NSUP+1)
      DOUBLE PRECISION WEIGHT(2),W(NSUP)
      INTEGER DEGREE,I,IL,HINFO(6),J,K,LIST,LSTNUM,LWIDTH,LWDTH1,
     *        LZNUM,MAXPSV,MINPSV,NLVL,NLVL1,NODES,NSTOP,NSTRT,NVARS,XLS
      DOUBLE PRECISION NWGHT(2)
      EXTERNAL MC60HD,MC60JD,MC60LD
      XLS = NSUP
      LIST = 2*NSUP + 1
      NWGHT(1) = WEIGHT(1)
      NWGHT(2) = WEIGHT(2)
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      INFO(4) = 0
      NVARS = 0
      LSTNUM = 0
      LZNUM = NSUP+1
      IF(JCNTL(2).NE.2) THEN
         DO 5 I = 1,NSUP
            PERMSV(I) = 1
    5    CONTINUE
      END IF
      DO 6 I = 1,NSUP
         K = ICPTR(I)
         DEGREE = ICPTR(I+1) - K
         IF (DEGREE.LE.1) THEN
           IF (DEGREE.EQ.0) THEN
             LZNUM = LZNUM - 1
             PERMSV(I) = -LZNUM
           ELSE IF (IRN(K).EQ.I) THEN
              LSTNUM = LSTNUM + 1
              PERMSV(I) = -LSTNUM
            END IF
         END IF
    6 CONTINUE
      DO 30 I = 1, NSUP
        IF (LSTNUM.GE.LZNUM-1) GO TO 35
        IF(JCNTL(2).EQ.0) THEN
          CALL MC60HD(N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +               IW(XLS+1),IW(LIST+1),HINFO)
          NSTRT = HINFO(1)
          PAIR(1,I) = HINFO(1)
          PAIR(2,I) = HINFO(2)
          NLVL = HINFO(3)
          LWIDTH = HINFO(4)
          NVARS = HINFO(5)
          NODES = HINFO(6)
        ELSE IF(JCNTL(2).EQ.1) THEN
          NSTOP = PAIR(1,I)
          CALL MC60LD(NSTOP,N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +            IW(XLS+1),NLVL,LWIDTH,NVARS)
          LWDTH1 = LWIDTH
          NLVL1 = NLVL
          NODES = IW(XLS+NLVL+1)-1
          NSTRT = NSTOP
          NSTOP = PAIR(2,I)
          CALL MC60LD(NSTOP,N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +            IW(XLS+1),NLVL,LWIDTH,NVARS)
          IF (NLVL1.GT.NLVL .OR.
     +             (NLVL1.EQ.NLVL .AND. LWDTH1.LT.LWIDTH)) THEN
             NSTRT = NSTOP
             NSTOP = PAIR(1,I)
             CALL MC60LD(NSTOP,N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +            IW(XLS+1),NLVL,LWIDTH,NVARS)
          END IF
        ELSE
          MAXPSV = 0
          DO 10 J = 1,NSUP
            IF (PERMSV(J).GT.0) THEN
               IF (PERMSV(J).GT.MAXPSV) THEN
                  IF (MAXPSV.EQ.0) THEN
                    MINPSV = PERMSV(J)
                    NSTRT = J
                  END IF
                  MAXPSV = PERMSV(J)
               END IF
               IF (PERMSV(J).LT.MINPSV) THEN
                  MINPSV = PERMSV(J)
                  NSTRT = J
               END IF
            END IF
   10     CONTINUE
          CALL MC60LD(NSTRT,N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,
     +            IW(XLS+1),NLVL,LWIDTH,NVARS)
          NODES = IW(XLS+NLVL+1)-1
          IF(MAXPSV.NE.MINPSV) THEN
            NWGHT(2) = (WEIGHT(2)*(NLVL-1))/(MAXPSV-MINPSV)
          ELSE
            NWGHT(2) = WEIGHT(2)
          END IF
        END IF
        INFO(1) = INFO(1) + 1
        IF (NVARS.GT.INFO(2)) THEN
           INFO(2) = NVARS
           INFO(3) = NLVL
           INFO(4) = LWIDTH
        END IF
        IF(JCNTL(1).EQ.1) THEN
           DO 11 J = NODES,1,-1
             LSTNUM = LSTNUM + 1
             PERMSV(IW(J)) = -LSTNUM
   11      CONTINUE
        ELSE
          IF(JCNTL(2).NE.2) THEN
            DO 15 IL = 1,NLVL
              DO 12 J = IW(XLS+IL), IW(XLS+IL+1) - 1
                PERMSV(IW(J)) = NLVL - IL
   12         CONTINUE
   15       CONTINUE
          END IF
          CALL MC60JD(NSUP,LIRN,NODES,NSTRT,LSTNUM,IRN,ICPTR,
     +         VARS,PERMSV,NWGHT,IW,IW(LIST+1),IW(XLS+1),W)
        END IF
   30 CONTINUE
   35 DO 40 I = 1,NSUP
          PERMSV(I) = -PERMSV(I)
   40 CONTINUE
      END
      SUBROUTINE MC60DD(N,NSUP,SVAR,VARS,PERMSV,PERM,POSSV)
      INTEGER N
      INTEGER NSUP
      INTEGER SVAR(N)
      INTEGER VARS(NSUP)
      INTEGER PERMSV(NSUP)
      INTEGER PERM(N)
      INTEGER POSSV(NSUP)
      INTEGER I,IS,L
      DO 10 IS = 1,NSUP
         PERM(PERMSV(IS)) = IS
   10 CONTINUE
      L = 1
      DO 20 I = 1,NSUP
         L = L + VARS(PERM(I))
         POSSV(PERM(I)) = L
   20 CONTINUE
      DO 30 I = 1,N
         IS = SVAR(I)
         L = POSSV(IS)-1
         POSSV(IS) = L
         PERM(I) = L
   30 CONTINUE
      END
      SUBROUTINE MC60ED(N,NSUP,LIRN,IRN,ICPTR,SVAR,VARS,PERMSV,PERM,IW)
      INTEGER N,NSUP,LIRN
      INTEGER IRN(LIRN),ICPTR(NSUP+1),PERM(N),PERMSV(NSUP),IW(NSUP),
     *        SVAR(N),VARS(NSUP)
      INTEGER I,IS,JS,K,L,M
      DO 10 IS = 1,NSUP
        JS = PERMSV(IS)
        IW(JS) = IS
        PERM(IS) = 0
   10 CONTINUE
      K = 0
      DO 30 L = 1,NSUP
        IS = IW(L)
        DO 20 M = ICPTR(IS), ICPTR(IS+1)-1
          JS = IRN(M)
          IF(PERM(JS).EQ.0) THEN
            K = K + 1
            PERM(JS) = K
            PERMSV(K) = JS
          END IF
   20   CONTINUE
        IF(PERM(IS).EQ.0) THEN
            K = K + 1
            PERM(IS) = K
            PERMSV(K) = IS
          END IF
   30 CONTINUE
      IF (N.EQ.NSUP) THEN
         DO 40 I = 1,N
            PERM(I) = PERMSV(I)
   40    CONTINUE
         RETURN
      END IF
      L = 1
      DO 45 IS = 1,NSUP
         JS = PERMSV(IS)
         L = L + VARS(JS)
         IW(JS) = L
   45 CONTINUE
      DO 50 I = 1,N
         IS = SVAR(I)
         L = IW(IS) - 1
         IW(IS) = L
         PERM(L) = I
   50 CONTINUE
      END
      SUBROUTINE MC60FD(N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,RINFO)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER LIRN,N,NSUP
      INTEGER IRN(LIRN),IW(2*NSUP+1),PERMSV(NSUP),ICPTR(NSUP+1),
     *        VARS(NSUP)
      DOUBLE PRECISION RINFO(4)
      INTEGER I,IMIN,J,JSTOP,JSTRT,K,NACTIV,NBR,NV
      INTRINSIC ABS,DBLE,MAX,MIN,SQRT
      DO 10 I = 1,4
         RINFO(I) = ZERO
  10  CONTINUE
      NACTIV = 0
      IW(NSUP+1) = 0
      DO 30 I = 1,NSUP
         J = PERMSV(I)
         IW(J) = I
   30 CONTINUE
      DO 80 I = 1,NSUP
         J = ABS(IW(I))
         NV = VARS(J)
         IW(NSUP+I+1) = IW(NSUP+I) + NV
         JSTRT = ICPTR(J)
         JSTOP = ICPTR(J+1) - 1
         IF(JSTRT.GE.JSTOP) THEN
           IF(JSTRT.GT.JSTOP) THEN
             GO TO 80
           ELSE
             IF(IRN(JSTRT).EQ.J)THEN
               RINFO(2) = MAX(RINFO(2),DBLE(NV))
               RINFO(3) = MAX(RINFO(3),DBLE(NV))
               DO 40 K = 1,NV
                 RINFO(1) = RINFO(1) + DBLE(K)
                 RINFO(4) = RINFO(4) + DBLE(K)**2
  40           CONTINUE
               GO TO 80
             END IF
           END IF
         END IF
         IMIN = I
         DO 50 K = JSTRT,JSTOP
            NBR = IRN(K)
            IMIN = MIN(IMIN,PERMSV(NBR))
            IF (IW(NBR).GT.0) THEN
               NACTIV = NACTIV +  VARS(NBR)
               IW(NBR) = -IW(NBR)
            END IF
   50    CONTINUE
         RINFO(3) = MAX(RINFO(3),DBLE(IW(NSUP+I+1)-IW(NSUP+IMIN)))
         IF (IW(J).GT.0) THEN
           IW(J) = -IW(J)
           RINFO(2) = MAX(RINFO(2),DBLE(NACTIV+1))
           RINFO(1) = RINFO(1) + NV*DBLE(NACTIV+1)
           RINFO(4) = RINFO(4) + NV*DBLE(NACTIV+1)**2
         ELSE
           RINFO(2) = MAX(RINFO(2),DBLE(NACTIV))
           DO 70 J = 1, NV
             RINFO(1) = RINFO(1) + DBLE(NACTIV)
             RINFO(4) = RINFO(4) + DBLE(NACTIV)**2
             NACTIV = NACTIV - 1
   70     CONTINUE
        END IF
   80 CONTINUE
      RINFO(4) = SQRT(RINFO(4)/DBLE(N))
      END
      SUBROUTINE MC60GD(N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,RINFO)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER LIRN,N,NSUP
      DOUBLE PRECISION RINFO(4)
      INTEGER ICPTR(NSUP+1),IRN(LIRN),IW(NSUP),PERMSV(NSUP),VARS(NSUP)
      INTEGER I,IEQ,J,JSTOP,JSTRT,K,KFRNT,LFRNT,MFR,NV
      INTRINSIC DBLE,MAX,SQRT
      DO 10 I = 1,4
         RINFO(I) = ZERO
   10 CONTINUE
      DO 20 I = 1,NSUP
         IW(I) = 0
   20 CONTINUE
      KFRNT = 0
      LFRNT = 0
      DO 40 IEQ = 1,NSUP
         I = PERMSV(IEQ)
         IW(I) = IEQ
         JSTRT = ICPTR(I)
         JSTOP = ICPTR(I+1) - 1
         DO 30 J = JSTRT,JSTOP
            MFR = IRN(J)
            IW(MFR) = IEQ
   30    CONTINUE
   40 CONTINUE
      DO 90 IEQ = 1,NSUP
         I = PERMSV(IEQ)
         JSTRT = ICPTR(I)
         JSTOP = ICPTR(I+1) - 1
         IF (JSTRT.GT.JSTOP) GO TO 90
         NV = VARS(I)
         DO 50 K = 1,NV
            KFRNT = KFRNT + 1
            RINFO(3) = RINFO(3) + DBLE(KFRNT)**2
   50    CONTINUE
         RINFO(1) = MAX(RINFO(1),DBLE(KFRNT))
         IF (IW(I).GE.0) THEN
            LFRNT = LFRNT + NV
            IW(I) = -IW(I)
         END IF
         DO 60 J = JSTRT,JSTOP
            MFR = IRN(J)
            IF (IW(MFR).GE.0) THEN
               NV = VARS(MFR)
               LFRNT = LFRNT + NV
               IW(MFR) = -IW(MFR)
            END IF
   60    CONTINUE
         RINFO(2) = MAX(RINFO(2),DBLE(LFRNT))
         IF (-IW(I).EQ.IEQ) THEN
            IW(I) = 0
            NV = VARS(I)
            DO 65 K = 1,NV
               RINFO(4) = RINFO(4) + DBLE(LFRNT)**2
               LFRNT = LFRNT - 1
               KFRNT = KFRNT - 1
   65       CONTINUE
         END IF
         DO 80 J = JSTRT,JSTOP
            MFR = IRN(J)
            IF (-IW(MFR).EQ.IEQ) THEN
               NV = VARS(MFR)
               DO 70 K = 1,NV
                  RINFO(4) = RINFO(4) + DBLE(LFRNT)**2
                  LFRNT = LFRNT - 1
                  KFRNT = KFRNT - 1
   70          CONTINUE
            END IF
   80    CONTINUE
   90 CONTINUE
      RINFO(3) = SQRT(RINFO(3)/DBLE(N))
      RINFO(4) = SQRT(RINFO(4)/DBLE(N))
      END
      SUBROUTINE MC60HD(N,NSUP,LIRN,IRN,ICPTR,VARS,MASK,LS,XLS,LIST,
     +                 INFO)
      INTEGER N,NSUP,LIRN
      INTEGER ICPTR(NSUP+1),IRN(LIRN),LIST(NSUP),LS(NSUP),
     +        MASK(NSUP),VARS(NSUP),XLS(NSUP+1),INFO(6)
      INTEGER DEGREE,I,J,LSIZE,LWIDTH,MAIN,MAXDEP,
     +        MINDEG,MINWID,NLSIZE,NLVL,NODE,NODES,NSTOP,NSTRT,NVARS
      EXTERNAL MC60LD
      MINDEG = N+1
      INFO(5) = 0
      DO 10 I = 1,NSUP
          IF (MASK(I).EQ.1) THEN
            INFO(5) = INFO(5) + VARS(I)
            DEGREE = ICPTR(I+1) - ICPTR(I)
            IF (DEGREE.LE.MINDEG) THEN
              IF (DEGREE.LT.MINDEG) THEN
                  NSTRT = I
                  MINDEG = DEGREE
              END IF
            END IF
          END IF
   10 CONTINUE
      CALL MC60LD(NSTRT,N,NSUP,LIRN,IRN,ICPTR,VARS,MASK,LS,XLS,
     +            MAXDEP,LWIDTH,NVARS)
      NODES = XLS(MAXDEP+1) - 1
      NSTOP = 0
      DO 70 MAIN = 1, NODES
        INFO(4) = LWIDTH
        LSIZE = 0
        DO 30 I = XLS(MAXDEP),XLS(MAXDEP+1) - 1
          NODE = LS(I)
          LSIZE = LSIZE + 1
          LIST(LSIZE) = NODE
          XLS(NODE) = ICPTR(NODE+1) - ICPTR(NODE)
   30   CONTINUE
        DO 50 NLSIZE = 1,5
           MINDEG = N+1
           DO 41 I = NLSIZE,LSIZE
             IF(XLS(LIST(I)).LT.MINDEG)THEN
                J = I
                MINDEG = XLS(LIST(I))
             END IF
   41     CONTINUE
          IF(MINDEG.EQ.N+1) GO TO 55
          NODE = LIST(J)
          LIST(J) = LIST(NLSIZE)
          LIST(NLSIZE) = NODE
          DO 42 I = ICPTR(NODE), ICPTR(NODE+1)-1
             XLS(IRN(I)) = N+1
   42     CONTINUE
   50   CONTINUE
   55   NLSIZE = NLSIZE-1
        MINWID = N
        DO 60 I = 1,NLSIZE
          NODE = LIST(I)
          CALL MC60LD(NODE,MINWID,NSUP,LIRN,IRN,ICPTR,VARS,MASK,LS,XLS,
     +               NLVL,LWIDTH,NVARS)
          IF (LWIDTH.LT.MINWID) THEN
            IF (NLVL.GT.MAXDEP) THEN
              NSTRT = NODE
              MAXDEP = NLVL
              GO TO 70
            ELSE
              NSTOP = NODE
              MINWID = LWIDTH
            END IF
          END IF
   60   CONTINUE
        GO TO 80
   70 CONTINUE
   80 IF (INFO(4) .LT. MINWID) THEN
         INFO(1) = NSTRT
         NSTRT = NSTOP
         NSTOP = INFO(1)
      END IF
      IF(NSTOP.NE.NODE) CALL MC60LD(NSTOP,N,NSUP,LIRN,IRN,ICPTR,VARS,
     +               MASK,LS,XLS,NLVL,LWIDTH,NVARS)
      INFO(1) = NSTRT
      INFO(2) = NSTOP
      INFO(3) = MAXDEP
      INFO(4) = LWIDTH
      INFO(5) = NVARS
      INFO(6) = NODES
      END
      SUBROUTINE MC60JD(NSUP,LIRN,NODES,NSTRT,LSTNUM,IRN,ICPTR,VARS,
     +                  STATUS,WEIGHT,NLIST,QUEUE,DEG,PRIOR)
      INTEGER LIRN,LSTNUM,NSUP,NODES,NSTRT
      INTEGER DEG(NSUP),ICPTR(NSUP+1),IRN(LIRN),NLIST(NSUP),
     +        QUEUE(0:NODES-1),STATUS(NSUP),VARS(NSUP)
      DOUBLE PRECISION WEIGHT(2),PRIOR(NSUP)
      INTEGER ADDRES,DEGREE,FATHER,I,ISTOP,ISTRT,J,JSTOP,JSTRT,J1,J2,
     +        K,L,NABOR,NBR,NEXT,NODE,NQ,QNODE,SON,THRESH
      PARAMETER (THRESH=100)
      DOUBLE PRECISION MAXPRT,PNODE,PRTY
      DO 10 I = 1,NODES
         NODE = NLIST(I)
         DEGREE = VARS(NODE)
         K = DEGREE
         VARS(NODE) = 0
         DO 7 J = ICPTR(NODE),ICPTR(NODE+1) - 1
            DEGREE = DEGREE + VARS(IRN(J))
    7    CONTINUE
         VARS(NODE) = K
         PRIOR(NODE) = -WEIGHT(1)*DEGREE-WEIGHT(2)*STATUS(NODE)
         STATUS(NODE) = 2
         DEG(NODE) = DEGREE
   10 CONTINUE
      NQ = 1
      QUEUE(NQ) = NSTRT
      QUEUE(0) = NSTRT
      STATUS(NSTRT) = 1
      PRIOR(NSTRT) = 1.0E30
      DO 70 L = 1, NODES
         IF(NQ.GT.THRESH)GO TO 100
         ADDRES = 1
         MAXPRT = PRIOR(QUEUE(1))
         DO 30 I = 2,NQ
            PRTY = PRIOR(QUEUE(I))
            IF (PRTY.GT.MAXPRT) THEN
               ADDRES = I
               MAXPRT = PRTY
            END IF
   30    CONTINUE
         NEXT = QUEUE(ADDRES)
         QUEUE(ADDRES) = QUEUE(NQ)
         NQ = NQ - 1
         ISTRT = ICPTR(NEXT)
         ISTOP = ICPTR(NEXT+1) - 1
         IF (STATUS(NEXT).EQ.1) THEN
            DO 40 I = ISTRT,ISTOP
               NBR = IRN(I)
               PRIOR(NBR) = PRIOR(NBR) + WEIGHT(1)*VARS(NEXT)
               DEG(NBR) = DEG(NBR) - VARS(NEXT)
               IF(DEG(NBR).EQ.0) PRIOR(NBR) = 1.0E30
               IF (STATUS(NBR).EQ.2) THEN
                  NQ = NQ + 1
                  QUEUE(NQ) = NBR
                  STATUS(NBR) = 1
                  PRIOR(NBR) = PRIOR(NBR)
               END IF
   40       CONTINUE
         END IF
         LSTNUM = LSTNUM + 1
         STATUS(NEXT) = -LSTNUM
         DO 60 I = ISTRT,ISTOP
            NBR = IRN(I)
            IF (STATUS(NBR).NE.1) GO TO 60
            PRIOR(NBR) = PRIOR(NBR) + WEIGHT(1)*VARS(NBR)
            STATUS(NBR) = -1
            DEG(NBR) = DEG(NBR) - VARS(NBR)
            IF(DEG(NBR).EQ.0) PRIOR(NBR) = 1.0E30
            JSTRT = ICPTR(NBR)
            JSTOP = ICPTR(NBR+1) - 1
            DO 50 J = JSTRT,JSTOP
               NABOR = IRN(J)
               IF (STATUS(NABOR).LT.0) GO TO 50
               PRIOR(NABOR) = PRIOR(NABOR) + WEIGHT(1)*VARS(NBR)
               DEG(NABOR) = DEG(NABOR) - VARS(NBR)
               IF(DEG(NABOR).EQ.0) PRIOR(NABOR) = 1.0E30
               IF (STATUS(NABOR).EQ.2) THEN
                  NQ = NQ + 1
                  QUEUE(NQ) = NABOR
                  STATUS(NABOR) = 1
               END IF
   50       CONTINUE
            STATUS(NBR) = 0
   60    CONTINUE
   70 CONTINUE
      RETURN
  100  DO 120 I = 1, NQ
         NBR = QUEUE(I)
         PNODE = PRIOR(NBR)
         J1 = I
         DO 116 K = 1, J1
            J2 = J1/2
            FATHER = QUEUE(J2)
            IF(PRIOR(FATHER).GE.PNODE) GO TO 118
            QUEUE(J1) = FATHER
            NLIST(FATHER) = J1
            J1 = J2
  116   CONTINUE
  118   QUEUE(J1) = NBR
        NLIST(NBR) = J1
  120 CONTINUE
      I = L
      DO 170 L =I, NODES
         NEXT = QUEUE(1)
         QNODE = QUEUE(NQ)
         PNODE = PRIOR(QNODE)
         NQ = NQ - 1
         J = 2
         J2 = 1
         IF(NQ.GT.1) QUEUE(NQ+1) = QUEUE(NQ)
         DO 125 I = 2, NQ
            IF(J.GT.NQ) GO TO 130
            IF( PRIOR(QUEUE(J)).LT.PRIOR(QUEUE(J+1)) ) J=J+1
            SON = QUEUE(J)
            IF(PNODE.GE.PRIOR(SON)) GO TO 130
            QUEUE(J2) = SON
            NLIST(SON) = J2
            J2 = J
            J = J*2
  125    CONTINUE
  130    QUEUE(J2) = QNODE
         NLIST(QNODE) = J2
         ISTRT = ICPTR(NEXT)
         ISTOP = ICPTR(NEXT+1) - 1
         IF (STATUS(NEXT).EQ.1) THEN
            DO 140 I = ISTRT,ISTOP
               NBR = IRN(I)
               IF (NBR.EQ.NEXT) GO TO 140
               PRIOR(NBR) = PRIOR(NBR) + WEIGHT(1)*VARS(NEXT)
               DEG(NBR) = DEG(NBR) - VARS(NEXT)
               IF(DEG(NBR).EQ.0) PRIOR(NBR) = 1.0E30
               IF (STATUS(NBR).EQ.2) THEN
                  NQ = NQ + 1
                  QUEUE(NQ) = NBR
                  STATUS(NBR) = 1
                  NLIST(NBR) = NQ
               END IF
               PNODE = PRIOR(NBR)
               J = NLIST(NBR)
               DO 133 K = 1, NQ
                  J2 = J/2
                  FATHER = QUEUE(J2)
                  IF(PRIOR(FATHER).GE.PNODE) GO TO 137
                  QUEUE(J) = FATHER
                  NLIST(FATHER) = J
                  J = J2
  133          CONTINUE
  137          QUEUE(J) = NBR
               NLIST(NBR) = J
  140       CONTINUE
         END IF
         LSTNUM = LSTNUM + 1
         STATUS(NEXT) = -LSTNUM
         DO 160 I = ISTRT,ISTOP
            NBR = IRN(I)
            IF (STATUS(NBR).NE.1) GO TO 160
            PRIOR(NBR) = PRIOR(NBR) + WEIGHT(1)*VARS(NBR)
            STATUS(NBR) = -1
            DEG(NBR) = DEG(NBR) - VARS(NBR)
            IF(DEG(NBR).EQ.0) PRIOR(NBR) = 1.0E30
            PNODE = PRIOR(NBR)
            J = NLIST(NBR)
            DO 142 K = 1, NQ
               J2 = J/2
               FATHER = QUEUE(J2)
               IF(PRIOR(FATHER).GE.PNODE) GO TO 144
               QUEUE(J) = FATHER
               NLIST(FATHER) = J
               J = J2
  142       CONTINUE
  144       QUEUE(J) = NBR
            NLIST(NBR) = J
            JSTRT = ICPTR(NBR)
            JSTOP = ICPTR(NBR+1) - 1
            DO 150 J = JSTRT,JSTOP
               NABOR = IRN(J)
               IF (STATUS(NABOR).LT.0) GO TO 150
               PRIOR(NABOR) = PRIOR(NABOR) + WEIGHT(1)*VARS(NBR)
               DEG(NABOR) = DEG(NABOR) - VARS(NBR)
               IF(DEG(NABOR).EQ.0) PRIOR(NABOR) = 1.0E30
               IF (STATUS(NABOR).EQ.2) THEN
                  NQ = NQ + 1
                  QUEUE(NQ) = NABOR
                  STATUS(NABOR) = 1
                  NLIST(NABOR) = NQ
               END IF
               PNODE = PRIOR(NABOR)
               J1 = NLIST(NABOR)
               J2 = J1/2
               FATHER = QUEUE(J2)
               IF(PRIOR(FATHER).GE.PNODE) GO TO 148
               QUEUE(J1) = FATHER
               NLIST(FATHER) = J1
               J1 = J2
               DO 146 K = 2, NQ
                  J2 = J1/2
                  FATHER = QUEUE(J2)
                  IF(PRIOR(FATHER).GE.PNODE) GO TO 148
                  QUEUE(J1) = FATHER
                  NLIST(FATHER) = J1
                  J1 = J2
  146          CONTINUE
  148          QUEUE(J1) = NABOR
               NLIST(NABOR) = J1
  150       CONTINUE
            STATUS(NBR) = 0
  160    CONTINUE
  170 CONTINUE
      END
      SUBROUTINE MC60LD(ROOT,MAXWID,NSUP,LIRN,IRN,ICPTR,VARS,MASK,LS,
     +                 XLS,NLVL,LWIDTH,NVARS)
      INTEGER LWIDTH,MAXWID,NSUP,NLVL,LIRN,ROOT,NVARS
      INTEGER ICPTR(NSUP+1),IRN(LIRN),LS(NSUP),MASK(NSUP),
     *        VARS(NSUP),XLS(NSUP+1)
      INTEGER I,J,LBEGIN,LNBR,LVLEND,LW,NBR,NODE
      INTRINSIC MAX
      MASK(ROOT) = -MASK(ROOT)
      LS(1) = ROOT
      LVLEND = 0
      NVARS = 0
      LNBR = 1
      LWIDTH = VARS(ROOT)
      DO 35 NLVL = 1,NSUP
          LBEGIN = LVLEND + 1
          LVLEND = LNBR
          XLS(NLVL) = LBEGIN
          LW = 0
          DO 30 I = LBEGIN,LVLEND
              NODE = LS(I)
              DO 20 J = ICPTR(NODE), ICPTR(NODE+1) - 1
                  NBR = IRN(J)
                  IF (MASK(NBR).GT.0) THEN
                      LNBR = LNBR + 1
                      LS(LNBR) = NBR
                      MASK(NBR) = -MASK(NBR)
                      LW = LW + VARS(NBR)
                  END IF
   20         CONTINUE
   30     CONTINUE
          LWIDTH = MAX(LW, LWIDTH)
          NVARS = NVARS + LW
          IF (LNBR.EQ.LVLEND) GO TO 40
          IF (LWIDTH.GE.MAXWID) GO TO 40
   35 CONTINUE
   40 XLS(NLVL+1) = LVLEND + 1
      DO 50 I = 1,LNBR
          MASK(LS(I)) = ABS(MASK(LS(I)))
   50 CONTINUE
      END
      SUBROUTINE MC60OD(N,NC,LIRN,IRN,ICPTR,SVAR,NSUP,NEW,VARS,FLAG)
      INTEGER N
      INTEGER NC
      INTEGER LIRN
      INTEGER IRN(LIRN)
      INTEGER ICPTR(NC+1)
      INTEGER SVAR(N)
      INTEGER NSUP
      INTEGER NEW(N)
      INTEGER VARS(N)
      INTEGER FLAG(N)
      INTEGER I,IS,J,JS,K,K1,K2
      DO 10 I = 1,N
         SVAR(I) = 1
   10 CONTINUE
      VARS(1) = N
      FLAG(1) = 0
      NSUP = 1
      DO 40 J = 1,NC
         K1 = ICPTR(J)
         K2 = ICPTR(J+1) - 1
         DO 20 K = K1,K2
            IS = SVAR(IRN(K))
            VARS(IS) = VARS(IS) - 1
   20    CONTINUE
         DO 30 K = K1,K2
            I = IRN(K)
            IS = SVAR(I)
            IF (FLAG(IS).LT.J) THEN
               FLAG(IS) = J
               IF (VARS(IS).GT.0) THEN
                  NSUP = NSUP + 1
                  VARS(NSUP) = 1
                  FLAG(NSUP) = J
                  NEW(IS) = NSUP
                  SVAR(I) = NSUP
               ELSE
                  VARS(IS) = 1
                  NEW(IS) = IS
                  SVAR(I) = IS
               END IF
            ELSE
               JS = NEW(IS)
               VARS(JS) = VARS(JS) + 1
               SVAR(I) = JS
            END IF
   30    CONTINUE
   40 CONTINUE
      END
      SUBROUTINE MC60PD(N,NSUP,LIRN,IRN,ICPTR,SVAR,VARS,VAR,FLAG)
      INTEGER N
      INTEGER LIRN
      INTEGER IRN(LIRN)
      INTEGER NSUP
      INTEGER ICPTR(N+1)
      INTEGER SVAR(N)
      INTEGER VARS(NSUP)
      INTEGER VAR(NSUP)
      INTEGER FLAG(NSUP)
      INTEGER I,IS,J,JS,K,K1,L
      DO 10 IS = 1,NSUP
         FLAG(IS) = -1
         VARS(IS) = 1
   10 CONTINUE
      L = 1
      DO 20 I = 1,N
         IS = SVAR(I)
         JS = FLAG(IS)
         IF(JS.GT.0)THEN
            SVAR(I) = JS
            VARS(JS) = VARS(JS) + 1
         ELSE IF(JS.LT.0)THEN
            FLAG(IS) = L
            VAR(L) = I
            SVAR(I) = L
            L = L + 1
         END IF
   20 CONTINUE
      DO 30 IS = 1,NSUP
         FLAG(IS) = 0
   30 CONTINUE
      L = 1
      K1 = 1
      DO 60 JS = 1,NSUP
         J = VAR(JS)
         K1 = ICPTR(J)
         ICPTR(JS) = L
         DO 50 K = K1, ICPTR(J+1)-1
            IS = SVAR(IRN(K))
            IF(FLAG(IS).NE.JS)THEN
               FLAG(IS) = JS
               IRN(L) = IS
               L = L+1
            END IF
   50    CONTINUE
   60 CONTINUE
      ICPTR(JS) = L
      END
* *******************************************************************
* COPYRIGHT (c) 1998 Council for the Central Laboratory
*                    of the Research Councils
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
* Licence, see http://hsl.rl.ac.uk/acuk/cou.html
*
* Please note that for a UK ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between AEA
*    Technology plc and the Licensee on suitable terms and conditions,
*    which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor AEA Technology plc shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
*######DATE 11 February 1998

      SUBROUTINE MC63ID(ICNTL)
      INTEGER ICNTL(10)
      INTEGER I
      ICNTL(1) = 6
      ICNTL(2) = 6
      DO 10 I = 3,10
         ICNTL(I) = 0
   10 CONTINUE
      END
C**********************************************************************
      SUBROUTINE MC63AD(DIRECT,N,NELT,NE,ELTVAR,ELTPTR,ORDER,PERM,NSUP,
     +                 VARS,SVAR,WT,LIW,IW,LW,W,ICNTL,INFO,RINFO)
C*********************************************************************
      DOUBLE PRECISION ONE,TWO,ZERO
      PARAMETER (ONE=1.0D0,TWO=2.0D0,ZERO=0.0D0)
      INTEGER LIW,LW,N,NE,NELT,NSUP
      LOGICAL DIRECT
      DOUBLE PRECISION RINFO(6),W(LW),WT(3)
      INTEGER ELTPTR(NELT+1),ELTVAR(NE),ICNTL(10),INFO(15),IW(LIW),
     +        ORDER(NELT),PERM(*),SVAR(N),VARS(N)
      INTEGER FLAG,I,IERR,IOUT,IPERM,IPTR,IRN,IW1,IWORK,J,JCNTL,K,KEY,
     +        KSUP,LENIP,LENIRN,LP,M,MP,NEED,NEW,NNE,NODES,NZ,PAIR,
     +        SELTPR,SELTVR,USED
      DOUBLE PRECISION WEIGHT(2)
      INTEGER ICNL60(2),INFO60(4),JNFO(4)
      EXTERNAL MC60CD,MC60OD,MC63BD,MC63CD,MC63DD,MC63FD,MC63GD
      INTRINSIC INT,MAX
      DO 10 I = 1,15
         INFO(I) = 0
   10 CONTINUE
      DO 20 I = 1,6
         RINFO(I) = ZERO
   20 CONTINUE
      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (N.LT.1) THEN
         INFO(1) = -1
         IF (LP.GT.0) THEN
            WRITE (LP,FMT='(/A,I2)')
     +    ' Error message from MC63A/AD: INFO(1) = ',INFO(1)
            WRITE (LP,FMT='(A,I8)') ' Value of N out-of-range.  N = '
     +           ,N
         END IF
         GO TO 120
      END IF
      IF (NELT.LT.1) THEN
         INFO(1) = -1
         IF (LP.GT.0) THEN
            WRITE (LP,FMT='(/A,I2)')
     +    ' Error message from MC63A/AD: INFO(1) = ',INFO(1)
            WRITE (LP,FMT='(A,I8)')
     +    ' Value of NELT out-of-range.  NELT = ',NELT
         END IF
         GO TO 120
      END IF
      IF (NE.LT.ELTPTR(NELT+1)-1) THEN
         INFO(1) = -1
         IF (LP.GT.0) THEN
            WRITE (LP,FMT='(/A,I2)')
     +    ' Error message from MC63A/AD: INFO(1) = ',INFO(1)
            WRITE (LP,FMT='(A)')
     +    ' Value of NE is less than ELTPTR(NELT+1)-1'
         END IF
         GO TO 120
      END IF
      IF (DIRECT) THEN
         INFO(6) = NELT
         IF (LW.LT.NELT) INFO(1) = -2
      ELSE
         INFO(6) = N
         IF (ICNTL(3).EQ.0 .AND. LW.LT.1) INFO(1) = -2
         IF (ICNTL(3).EQ.1 .AND. LW.LT.N) INFO(1) = -2
      END IF
      IF (INFO(1).EQ.-2) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT='(/A,I2)')
     +        ' Error message from MC63A/AD: INFO(1) = ',INFO(1)
            WRITE (LP,FMT='(A,I8)') ' Sufficient value for LW is ',
     +        INFO(6)
         END IF
         GO TO 120
      END IF
      INFO(5) = N + NELT
      IF (ICNTL(3).EQ.0) THEN
         INFO(5) = MAX(2*N,INFO(5))
         IF (DIRECT) INFO(5) = MAX(NE+3*NELT+3,INFO(5))
         IF (.NOT.DIRECT) INFO(5) = MAX(NE+NELT+5,INFO(5))
      ELSE
         IF (DIRECT) INFO(5) = 2* (NELT+N+1) + MAX(NE,4*NELT)
         IF (.NOT.DIRECT) INFO(5) = 3*N + 2 + NELT + MAX(NE,3*N)
      END IF
      IF (LIW.LT.INFO(5)) THEN
         INFO(1) = -4
         IF (DIRECT) INFO(5) = MAX(NE,4*NELT) + NELT* (1+NELT) +
     +                         2*MAX(NELT,N) + 3
         IF (.NOT.DIRECT) THEN
            NODES = 0
            DO 30 K = 1,NELT
               NODES = MAX(NODES,ELTPTR(K+1)-ELTPTR(K))
   30       CONTINUE
            INFO(5) = 3*N + 2 + NELT*NODES*NODES + MAX(3*N,NE)
         END IF
         IF (LP.GE.0) THEN
            WRITE (LP,FMT='(/A,I2)')
     +        ' Error message from MC63A/AD: INFO(1) = ',INFO(1)
            WRITE (LP,FMT='(A,A)')
     +        ' The workspace provided by the user is less than ',
     +        'the minimum required'
            WRITE (LP,FMT='(A,I8,A,I8)') ' LIW =  ',LIW,
     +        ' is too small. Sufficient value ',INFO(5)
         END IF
         GO TO 120
      END IF
      DO 40 I = 1,NELT
         IW(I) = I
   40 CONTINUE
      IF (ICNTL(5).EQ.1) THEN
         JCNTL = -1
      ELSE
         JCNTL = 1
      END IF
      CALL MC63BD(JCNTL,N,N,NELT,NE,ELTVAR,ELTPTR,VARS,IW,IW(NELT+1),
     +           INFO,RINFO)
      IOUT = INT(INFO(2)+INFO(3))
      IF (ICNTL(5).EQ.1 .AND. IOUT.GT.0) GO TO 110
      IF (IOUT.GT.0) INFO(7) = 1
      IF (ICNTL(3).EQ.0) THEN
         NEED = 2*N
         NEW = 1
         FLAG = NEW + N
         CALL MC60OD(N,NELT,NE,ELTVAR,ELTPTR,SVAR,NSUP,IW(NEW),VARS,
     +              IW(FLAG))
         IF (.NOT.DIRECT) THEN
            INFO(6) = NSUP
            IF (LW.LT.NSUP) INFO(1) = -2
            IF (INFO(1).EQ.-2) THEN
               IF (LP.GT.0) THEN
                  WRITE (LP,FMT='(/A,I2)')
     +              ' Error message from MC63A/AD: INFO(1) = ',INFO(1)
                  WRITE (LP,FMT='(A,I8)') ' Minimum value for LW is ',
     +              INFO(6)
               END IF
               GO TO 120
            END IF
         END IF
      ELSE
         NSUP = N
      END IF
      IF (DIRECT) THEN
         LENIP = NELT + 1
         M = MAX(NELT,NSUP)
      ELSE
         LENIP = NSUP + 1
         M = NSUP
      END IF
      NEED = NE + LENIP + NSUP + M + 1
      INFO(5) = MAX(NEED,INFO(5))
      IF (LIW.LT.INFO(5)) THEN
         IERR = -1
         IF (DIRECT) THEN
            NNE = NELT*NELT
         ELSE
            NODES = 0
            DO 50 K = 1,NELT
               NODES = MAX(NODES,ELTPTR(K+1)-ELTPTR(K))
   50       CONTINUE
            NNE = NELT*NODES*NODES
         END IF
         GO TO 60
      END IF
      LENIRN = LIW - NEED
      IPTR = 1
      IRN = IPTR + LENIP
      SELTVR = IRN + LENIRN
      SELTPR = SELTVR + NE
      IWORK = SELTPR + NSUP + 1
      CALL MC63CD(DIRECT,N,NELT,NE,NSUP,ELTVAR,ELTPTR,IW(IRN),LENIRN,
     +           IW(IPTR),LENIP,SVAR,IW(SELTVR),IW(SELTPR),IW(IWORK),M,
     +           ICNTL,IERR,NNE)
      IF (ICNTL(3).EQ.0) INFO(7) = 1
      INFO(5) = MAX(INFO(5),NEED+NNE)
   60 USED = LENIP + NNE
      IF (DIRECT) THEN
         NEED = USED + 2*NSUP + 4*NELT + 2
         INFO(5) = MAX(NEED,INFO(5))
         IF (IERR.LT.0 .OR. LIW.LT.INFO(5)) GO TO 100
         IWORK = IRN + NNE
         IW1 = IWORK + 4*NELT + 2
         CALL MC63DD(NELT,NSUP,NNE,NE,ELTVAR,ELTPTR,VARS,IW(IRN),
     +              IW(IPTR),ORDER,IW(IWORK),IW(IW1),W,ICNTL,WT)
         DO 70 I = 1,NELT
            IW(I) = -ORDER(I)
   70    CONTINUE
         DO 80 I = 1,NELT
            K = IW(I)
            ORDER(K) = I
   80    CONTINUE
      ELSE
         NEED = USED + 5*NSUP + 1
         INFO(5) = MAX(NEED,INFO(5))
         IF (IERR.LT.0 .OR. LIW.LT.INFO(5)) GO TO 100
         IPERM = IRN + NNE
         IWORK = IPERM + NSUP
         PAIR = IWORK + 3*NSUP + 1
         ICNL60(1) = 0
         ICNL60(2) = 0
         IF (ICNTL(4).EQ.1) THEN
            ICNL60(2) = 2
            DO 90 K = 1,N
               KSUP = SVAR(K)
               IF (KSUP.EQ.0) GO TO 90
               J = PERM(K)
               IW(IPERM-1+KSUP) = J
   90       CONTINUE
         END IF
         IF (ICNTL(6).EQ.0) THEN
            IF (ICNTL(4).EQ.0) THEN
               WT(1) = TWO
               WT(2) = ONE
            ELSE
               WT(1) = ONE
               WT(2) = TWO
            END IF
         END IF
         WEIGHT(1) = WT(1)
         WEIGHT(2) = WT(2)
         CALL MC60CD(N,NSUP,NNE,IW(IRN),IW(IPTR),VARS,ICNL60,IW(IPERM),
     +              WEIGHT,IW(PAIR),INFO60,IW(IWORK),W)
         KEY = IPERM + NSUP
         CALL MC63FD(NSUP,NELT,NE,ELTVAR,ELTPTR,IW(IPERM),IW(KEY),ORDER)
         CALL MC63GD(NELT,ORDER,IW(KEY))
      END IF
      NZ = ELTPTR(NELT+1) - 1
      CALL MC63BD(0,N,NSUP,NELT,NZ,ELTVAR,ELTPTR,VARS,ORDER,IW,JNFO,
     +           RINFO(4))
      IF (INFO(2).GT.0) INFO(1) = 1
      IF (INFO(3).GT.0) INFO(1) = 1
      IF (RINFO(1).LE.RINFO(4)) INFO(1) = 1
      IF (RINFO(2).LE.RINFO(5)) INFO(1) = 1
      IF (RINFO(3).LE.RINFO(6)) INFO(1) = 1
      IF (MP.GT.0) THEN
         IF (INFO(1).GT.0) THEN
            WRITE (MP,FMT='(/A,I2)')
     +        '   Warning from MC63A/AD: INFO(1) = ',INFO(1)
            IF (INFO(2).GT.0) WRITE (MP,FMT='(/3X,I8,A)') INFO(2),
     +          ' duplicate entries found'
            IF (INFO(3).GT.0) WRITE (MP,FMT='(/3X,I8,A)') INFO(3),
     +          ' out-of-range entries found'
            IF (RINFO(1).LE.RINFO(4)) WRITE (MP,FMT='(A)')
     +          '   Reordering has not reduced maximum front size'
            IF (RINFO(2).LE.RINFO(5)) WRITE (MP,FMT='(A)')
     +          '   Reordering has not reduced r.m.s. front size'
            IF (RINFO(3).LE.RINFO(6)) WRITE (MP,
     +          FMT='(A)') '   Reordering has not reduced profile'
         END IF
      END IF
      GO TO 120
  100 INFO(1) = -4
      IF (LP.GE.0) THEN
         WRITE (LP,FMT='(/A,I2)')
     +     ' Error message from MC63A/AD: INFO(1) = ',INFO(1)
         WRITE (LP,FMT='(A,I8,A,I8)') ' LIW =  ',LIW,
     +     ' is too small. Sufficient value ',INFO(5)
      END IF
      GO TO 120
  110 INFO(1) = -3
      IF (LP.GE.0) THEN
         WRITE (LP,FMT='(/A,I2)')
     +     '   Error message from MC63A/AD: INFO(1) = ',INFO(1)
         IF (INFO(2).GT.0) WRITE (LP,FMT='(/3X,I8,A)') INFO(2),
     +       ' duplicate entries found'
         IF (INFO(3).GT.0) WRITE (LP,FMT='(/3X,I8,A)') INFO(3),
     +       ' out-of-range entries found'
      END IF
  120 RETURN
      END
C***********************************************************************
      SUBROUTINE MC63BD(ICNTL,N,NSUP,NELT,NE,ELTVAR,ELTPTR,VARS,ORDER,
     +                 IW,INFO,RINFO)
C*******************************************************************
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER ICNTL,N,NE,NELT,NSUP
      DOUBLE PRECISION RINFO(3)
      INTEGER ELTPTR(NELT+1),ELTVAR(NE),INFO(4),IW(NSUP),ORDER(NELT),
     +        VARS(NSUP)
      INTEGER I,ICOUNT,IDUP,IELT,IOUT,J,JJ,JSTOP,JSTRT,JVAR,KFRNT
      INTRINSIC MAX,DBLE,SQRT
      DO 10 I = 1,4
         INFO(I) = 0
   10 CONTINUE
      DO 20 I = 1,3
         RINFO(I) = ZERO
   20 CONTINUE
      DO 30 I = 1,NSUP
         IW(I) = 0
   30 CONTINUE
      IF (N.EQ.NSUP) THEN
         DO 40 I = 1,N
            VARS(I) = 0
   40    CONTINUE
      END IF
      IF (ICNTL.EQ.0) THEN
         IF (N.EQ.NSUP) THEN
            DO 60 I = 1,NELT
               IELT = ORDER(I)
               JSTRT = ELTPTR(IELT)
               JSTOP = ELTPTR(IELT+1) - 1
               DO 50 J = JSTRT,JSTOP
                  JVAR = ELTVAR(J)
                  IW(JVAR) = I
                  VARS(JVAR) = 1
   50          CONTINUE
   60       CONTINUE
         ELSE
            DO 80 I = 1,NELT
               IELT = ORDER(I)
               JSTRT = ELTPTR(IELT)
               JSTOP = ELTPTR(IELT+1) - 1
               DO 70 J = JSTRT,JSTOP
                  JVAR = ELTVAR(J)
                  IW(JVAR) = I
   70          CONTINUE
   80       CONTINUE
         END IF
      ELSE IF (ICNTL.EQ.-1) THEN
         IOUT = 0
         IDUP = 0
         IF (N.EQ.NSUP) THEN
            DO 100 I = 1,NELT
               IELT = ORDER(I)
               JSTRT = ELTPTR(IELT)
               JSTOP = ELTPTR(IELT+1) - 1
               DO 90 J = JSTRT,JSTOP
                  JVAR = ELTVAR(J)
                  IF (JVAR.LT.1 .OR. JVAR.GT.NSUP) THEN
                     IOUT = IOUT + 1
                  ELSE IF (IW(JVAR).LT.I) THEN
                     IW(JVAR) = I
                     VARS(JVAR) = 1
                  ELSE
                     IDUP = IDUP + 1
                  END IF
   90          CONTINUE
  100       CONTINUE
         ELSE
            DO 120 I = 1,NELT
               IELT = ORDER(I)
               JSTRT = ELTPTR(IELT)
               JSTOP = ELTPTR(IELT+1) - 1
               DO 110 J = JSTRT,JSTOP
                  JVAR = ELTVAR(J)
                  IF (JVAR.LT.1 .OR. JVAR.GT.NSUP) THEN
                     IOUT = IOUT + 1
                  ELSE IF (IW(JVAR).LT.I) THEN
                     IW(JVAR) = I
                  ELSE
                     IDUP = IDUP + 1
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
         IF (IOUT+IDUP.GT.0) THEN
            INFO(1) = -1
            INFO(2) = IDUP
            INFO(3) = IOUT
            RETURN
         END IF
      ELSE
         IOUT = 0
         IDUP = 0
         IF (N.EQ.NSUP) THEN
            DO 140 I = 1,NELT
               IELT = ORDER(I)
               JSTRT = ELTPTR(IELT)
               JSTOP = ELTPTR(IELT+1) - 1
               DO 130 J = JSTRT,JSTOP
                  JVAR = ELTVAR(J)
                  IF (JVAR.LT.1 .OR. JVAR.GT.NSUP) THEN
                     ELTVAR(J) = 0
                     IOUT = IOUT + 1
                  ELSE IF (IW(JVAR).LT.I) THEN
                     IW(JVAR) = I
                     VARS(JVAR) = 1
                  ELSE
                     ELTVAR(J) = 0
                     IDUP = IDUP + 1
                  END IF
  130          CONTINUE
  140       CONTINUE
         ELSE
            DO 160 I = 1,NELT
               IELT = ORDER(I)
               JSTRT = ELTPTR(IELT)
               JSTOP = ELTPTR(IELT+1) - 1
               DO 150 J = JSTRT,JSTOP
                  JVAR = ELTVAR(J)
                  IF (JVAR.LT.1 .OR. JVAR.GT.NSUP) THEN
                     ELTVAR(J) = 0
                     IOUT = IOUT + 1
                  ELSE IF (IW(JVAR).LT.I) THEN
                     IW(JVAR) = I
                  ELSE
                     ELTVAR(J) = 0
                     IDUP = IDUP + 1
                  END IF
  150          CONTINUE
  160       CONTINUE
         END IF
         IF (IOUT+IDUP.GT.0) THEN
            INFO(1) = 1
            INFO(2) = IDUP
            INFO(3) = IOUT
            ICOUNT = 1
            DO 180 IELT = 1,NELT
               JSTRT = ELTPTR(IELT)
               JSTOP = ELTPTR(IELT+1) - 1
               ELTPTR(IELT) = ICOUNT
               DO 170 J = JSTRT,JSTOP
                  JVAR = ELTVAR(J)
                  IF (JVAR.NE.0) THEN
                     ELTVAR(ICOUNT) = JVAR
                     ICOUNT = ICOUNT + 1
                  END IF
  170          CONTINUE
  180       CONTINUE
            ELTPTR(NELT+1) = ICOUNT
         END IF
      END IF
      KFRNT = 0
      DO 220 I = 1,NELT
         IELT = ORDER(I)
         JSTRT = ELTPTR(IELT)
         JSTOP = ELTPTR(IELT+1) - 1
         DO 190 J = JSTRT,JSTOP
            JVAR = ELTVAR(J)
            IF (IW(JVAR).GE.0) THEN
               IW(JVAR) = -IW(JVAR)
               KFRNT = KFRNT + VARS(JVAR)
            END IF
  190    CONTINUE
         RINFO(1) = MAX(RINFO(1),DBLE(KFRNT))
         DO 210 J = JSTRT,JSTOP
            JVAR = ELTVAR(J)
            IF (-IW(JVAR).EQ.I) THEN
               DO 200 JJ = 1,VARS(JVAR)
                  RINFO(2) = RINFO(2) + DBLE(KFRNT)**2
                  RINFO(3) = RINFO(3) + DBLE(KFRNT)
                  KFRNT = KFRNT - 1
  200          CONTINUE
               INFO(4) = INFO(4) + VARS(JVAR)
            END IF
  210    CONTINUE
  220 CONTINUE
      RINFO(2) = SQRT(RINFO(2)/DBLE(INFO(4)))
      RETURN
      END
C*****************************************************************
      SUBROUTINE MC63CD(DIRECT,N,NELT,NZ,NSUP,ELTVAR,ELTPTR,IRN,LENIRN,
     +                 IP,LENIP,SVAR,SELTVR,SELTPR,IW,M,ICNTL,IERR,NNE)
      INTEGER IERR,LENIP,LENIRN,M,N,NELT,NNE,NSUP,NZ
      LOGICAL DIRECT
      INTEGER ELTPTR(NELT+1),ELTVAR(NZ),ICNTL(10),IP(LENIP),IRN(LENIRN),
     +        IW(M),SELTPR(NSUP+1),SELTVR(NZ),SVAR(N)
      INTEGER I,ICOUNT,IELT,J,J1,J2,JJ,JSUP,JVAR,K,K1,K2,KADJ,KELT,KSUP
      INTRINSIC MAX
      IERR = 0
      DO 10 JSUP = 1,NSUP
         IW(JSUP) = 0
         SELTPR(JSUP) = 0
   10 CONTINUE
      IF (ICNTL(3).EQ.0) THEN
         ICOUNT = 1
         DO 30 IELT = 1,NELT
            J1 = ELTPTR(IELT)
            J2 = ELTPTR(IELT+1) - 1
            ELTPTR(IELT) = ICOUNT
            DO 20 J = J1,J2
               JJ = ELTVAR(J)
               JSUP = SVAR(JJ)
               IF (IW(JSUP).LT.IELT) THEN
                  ELTVAR(ICOUNT) = JSUP
                  SELTPR(JSUP) = SELTPR(JSUP) + 1
                  ICOUNT = ICOUNT + 1
                  IW(JSUP) = IELT
               END IF
   20       CONTINUE
   30    CONTINUE
         ELTPTR(NELT+1) = ICOUNT
      ELSE
         DO 40 I = 1,N
            SVAR(I) = 0
   40    CONTINUE
         DO 60 IELT = 1,NELT
            J1 = ELTPTR(IELT)
            J2 = ELTPTR(IELT+1) - 1
            DO 50 J = J1,J2
               JVAR = ELTVAR(J)
               IF (IW(JVAR).LT.IELT) THEN
                  SELTPR(JVAR) = SELTPR(JVAR) + 1
                  IW(JVAR) = IELT
                  SVAR(JVAR) = JVAR
               END IF
   50       CONTINUE
   60    CONTINUE
      END IF
      SELTPR(1) = 1 + SELTPR(1)
      DO 70 JSUP = 2,NSUP
         SELTPR(JSUP) = SELTPR(JSUP-1) + SELTPR(JSUP)
   70 CONTINUE
      SELTPR(NSUP+1) = SELTPR(NSUP)
      DO 90 IELT = 1,NELT
         J1 = ELTPTR(IELT)
         J2 = ELTPTR(IELT+1) - 1
         DO 80 J = J1,J2
            JSUP = ELTVAR(J)
            K = SELTPR(JSUP) - 1
            SELTVR(K) = IELT
            SELTPR(JSUP) = K
   80    CONTINUE
   90 CONTINUE
      IF (.NOT.DIRECT) GO TO 150
      DO 100 IELT = 1,NELT
         IW(IELT) = 0
  100 CONTINUE
      KADJ = 1
      DO 140 IELT = 1,NELT
         IP(IELT) = KADJ
         J1 = ELTPTR(IELT)
         J2 = ELTPTR(IELT+1) - 1
         DO 130 J = J1,J2
            JSUP = ELTVAR(J)
            K1 = SELTPR(JSUP)
            K2 = SELTPR(JSUP+1) - 1
            IF (KADJ+K2-K1+1.GT.LENIRN) THEN
               DO 110 K = K1,K2
                  KELT = SELTVR(K)
                  IF (KELT.NE.IELT .AND. IW(KELT).NE.IELT) THEN
                     IF (KADJ.LE.LENIRN) IRN(KADJ) = KELT
                     KADJ = KADJ + 1
                     IW(KELT) = IELT
                  END IF
  110          CONTINUE
            ELSE
               DO 120 K = K1,K2
                  KELT = SELTVR(K)
                  IF (KELT.NE.IELT .AND. IW(KELT).NE.IELT) THEN
                     IRN(KADJ) = KELT
                     KADJ = KADJ + 1
                     IW(KELT) = IELT
                  END IF
  120          CONTINUE
            END IF
  130    CONTINUE
  140 CONTINUE
      IP(NELT+1) = KADJ
      GO TO 210
  150 CONTINUE
      DO 160 JSUP = 1,NSUP
         IW(JSUP) = 0
  160 CONTINUE
      KADJ = 1
      DO 200 JSUP = 1,NSUP
         IP(JSUP) = KADJ
         J1 = SELTPR(JSUP)
         J2 = SELTPR(JSUP+1) - 1
         DO 190 J = J1,J2
            IELT = SELTVR(J)
            K1 = ELTPTR(IELT)
            K2 = ELTPTR(IELT+1) - 1
            IF (KADJ+K2-K1+1.GT.LENIRN) THEN
               DO 170 K = K1,K2
                  KSUP = ELTVAR(K)
                  IF (KSUP.NE.JSUP .AND. IW(KSUP).NE.JSUP) THEN
                     IF (KADJ.LE.LENIRN) IRN(KADJ) = KSUP
                     KADJ = KADJ + 1
                     IW(KSUP) = JSUP
                  END IF
  170          CONTINUE
            ELSE
               DO 180 K = K1,K2
                  KSUP = ELTVAR(K)
                  IF (KSUP.NE.JSUP .AND. IW(KSUP).NE.JSUP) THEN
                     IRN(KADJ) = KSUP
                     KADJ = KADJ + 1
                     IW(KSUP) = JSUP
                  END IF
  180          CONTINUE
            END IF
  190    CONTINUE
  200 CONTINUE
      IP(NSUP+1) = KADJ
  210 CONTINUE
      IF (KADJ-1.GT.LENIRN) IERR = -4
      NNE = MAX(1,KADJ-1)
      RETURN
      END
C******************************************************************
      SUBROUTINE MC63DD(NELT,NSUP,NNE,NE,ELTVAR,ELTPTR,VARS,IRN,ICPTR,
     +                 ORDER,IW,IW1,W,ICNTL,WT)
C***********************************************************************
      DOUBLE PRECISION FIVE,ONE,TEN,TWO,ZERO
      PARAMETER (FIVE=5.0D0,ONE=1.0D0,TEN=10.0D0,TWO=2.0D0,ZERO=0.0D0)
      INTEGER NE,NELT,NNE,NSUP
      DOUBLE PRECISION W(NELT),WT(3)
      INTEGER ELTPTR(NELT+1),ELTVAR(NE),ICNTL(10),ICPTR(NELT+1),
     +        IRN(NNE),IW(4*NELT+2),IW1(2*NSUP),ORDER(NELT),VARS(NSUP)
      INTEGER DEGREE,DIST,I,IELIM,IFRNT,IL,IVARS,J,LIST,LS,LSTNUM,
     +        LWIDTH,MAXPSV,MINPSV,NLIST,NLVL,NODES,NSTRT,NVARS,QUEUE,
     +        UNNUM,XLS
      EXTERNAL MC60HD,MC60LD,MC63ED
      DOUBLE PRECISION WEIGHT(3)
      INTEGER HINFO(6)
      IF (ICNTL(6).EQ.0) THEN
         IF (ICNTL(4).EQ.0) THEN
            WT(1) = TEN
            WT(2) = FIVE
            WT(3) = ONE
         ELSE
            WT(1) = ONE
            WT(2) = TWO
            WT(3) = ZERO
         END IF
      END IF
      LSTNUM = 0
      WEIGHT(1) = WT(1)
      WEIGHT(2) = WT(2)
      WEIGHT(3) = WT(3)
      IF (ICNTL(4).NE.1) THEN
         DO 10 I = 1,NELT
            ORDER(I) = 1
   10    CONTINUE
      END IF
      DO 20 I = 1,NELT
         DEGREE = ICPTR(I+1) - ICPTR(I)
         IF (DEGREE.EQ.0) THEN
            LSTNUM = LSTNUM + 1
            ORDER(I) = -LSTNUM
         END IF
   20 CONTINUE
      LS = 1
      XLS = LS + NELT
      LIST = XLS + NELT + 1
      IVARS = LIST + NELT
      DO 30 I = 1,NELT
         IW(IVARS+I-1) = 1
   30 CONTINUE
      DO 70 I = 1,NELT
         IF (LSTNUM.GE.NELT) GO TO 80
         IF (ICNTL(4).NE.1) THEN
            CALL MC60HD(NELT,NELT,NNE,IRN,ICPTR,IW(IVARS),ORDER,IW(LS),
     +                 IW(XLS),IW(LIST),HINFO)
            NSTRT = HINFO(1)
            NLVL = HINFO(3)
            NODES = HINFO(6)
            DO 50 IL = 1,NLVL
               DO 40 J = IW(XLS+IL-1),IW(XLS+IL) - 1
                  ORDER(IW(J)) = NLVL - IL
   40          CONTINUE
   50       CONTINUE
         ELSE
            MAXPSV = 0
            DO 60 J = 1,NELT
               IF (ORDER(J).GT.0) THEN
                  IF (ORDER(J).GT.MAXPSV) THEN
                     IF (MAXPSV.EQ.0) THEN
                        MINPSV = ORDER(J)
                        NSTRT = J
                     END IF
                     MAXPSV = ORDER(J)
                  END IF
                  IF (ORDER(J).LT.MINPSV) THEN
                     MINPSV = ORDER(J)
                     NSTRT = J
                  END IF
               END IF
   60       CONTINUE
            CALL MC60LD(NSTRT,NELT,NELT,NNE,IRN,ICPTR,IW(IVARS),ORDER,
     +                 IW(LS),IW(XLS),NLVL,LWIDTH,NVARS)
            NODES = IW(XLS+NLVL) - 1
            IF (MAXPSV.NE.MINPSV) WEIGHT(2) = (WEIGHT(2)* (NLVL-1))/
     +          (MAXPSV-MINPSV)
         END IF
         NLIST = 1
         QUEUE = NLIST + NELT
         DIST = QUEUE + NELT
         UNNUM = DIST + NELT
         IFRNT = 1
         IELIM = IFRNT + NSUP
         CALL MC63ED(NELT,NNE,NODES,NSUP,NSTRT,LSTNUM,VARS,NE,ELTVAR,
     +              ELTPTR,IRN,ICPTR,ORDER,IW(NLIST),IW(QUEUE),W,
     +              IW1(IFRNT),IW1(IELIM),IW(DIST),IW(UNNUM),WEIGHT)
   70 CONTINUE
   80 RETURN
      END
C***********************************************************************
      SUBROUTINE MC63ED(NELT,NNE,NODES,NSUP,NSTRT,LSTNUM,VARS,NE,ELTVAR,
     +                 ELTPTR,IRN,ICPTR,STATUS,NLIST,QUEUE,PRIOR,IFRNT,
     +                 IELIM,DIST,UNNUM,WEIGHT)
C***********************************************************************
      INTEGER THRESH
      PARAMETER (THRESH=100)
      INTEGER LSTNUM,NE,NELT,NNE,NODES,NSTRT,NSUP
      DOUBLE PRECISION PRIOR(NELT),WEIGHT(3)
      INTEGER DIST(NELT),ELTPTR(NELT+1),ELTVAR(NE),ICPTR(NELT+1),
     +        IELIM(NSUP),IFRNT(NSUP),IRN(NNE),NLIST(NELT),
     +        QUEUE(0:NODES-1),STATUS(NELT),UNNUM(NELT),VARS(NSUP)
      DOUBLE PRECISION MAXPRT,PNODE,PRTY,W1,W2,W3
      INTEGER ADDRES,FATHER,I,ISTOP,ISTRT,J,J1,J2,K,L,L1,LFRNT,NBR,NEXT,
     +        NFRNT,NGAIN,NODE,NQ,QNODE,SON
      W1 = WEIGHT(1)
      W2 = WEIGHT(2)
      W3 = WEIGHT(3)
      DO 10 K = 1,NSUP
         IELIM(K) = 0
         IFRNT(K) = 0
   10 CONTINUE
      DO 30 I = 1,NODES
         NODE = NLIST(I)
         DIST(NODE) = STATUS(NODE)
         UNNUM(NODE) = ICPTR(NODE+1) - ICPTR(NODE)
         STATUS(NODE) = 2
         J1 = ELTPTR(NODE)
         J2 = ELTPTR(NODE+1) - 1
         DO 20 J = J1,J2
            K = ELTVAR(J)
            IELIM(K) = IELIM(K) + 1
   20    CONTINUE
   30 CONTINUE
      NQ = 1
      QUEUE(NQ) = NSTRT
      PRIOR(NSTRT) = 1.0D30
      QUEUE(0) = NSTRT
      DO 80 L = 1,NODES
         IF (NQ.GT.THRESH) GO TO 90
         ADDRES = 1
         MAXPRT = PRIOR(QUEUE(1))
         DO 40 I = 2,NQ
            PRTY = PRIOR(QUEUE(I))
            IF (PRTY.GT.MAXPRT) THEN
               ADDRES = I
               MAXPRT = PRTY
            END IF
   40    CONTINUE
         NEXT = QUEUE(ADDRES)
         QUEUE(ADDRES) = QUEUE(NQ)
         NQ = NQ - 1
         LSTNUM = LSTNUM + 1
         STATUS(NEXT) = -LSTNUM
         J1 = ELTPTR(NEXT)
         J2 = ELTPTR(NEXT+1) - 1
         DO 50 J = J1,J2
            K = ELTVAR(J)
            IELIM(K) = IELIM(K) - 1
            IF (IFRNT(K).EQ.0) IFRNT(K) = LSTNUM
   50    CONTINUE
         ISTRT = ICPTR(NEXT)
         ISTOP = ICPTR(NEXT+1) - 1
         DO 70 I = ISTRT,ISTOP
            NBR = IRN(I)
            IF (STATUS(NBR).EQ.2) THEN
               NQ = NQ + 1
               QUEUE(NQ) = NBR
               STATUS(NBR) = 1
            END IF
            IF (STATUS(NBR).EQ.1) THEN
               NFRNT = 0
               LFRNT = 0
               J1 = ELTPTR(NBR)
               J2 = ELTPTR(NBR+1) - 1
               DO 60 J = J1,J2
                  K = ELTVAR(J)
                  IF (IFRNT(K).EQ.0) NFRNT = NFRNT + VARS(K)
                  IF (IELIM(K).EQ.1) LFRNT = LFRNT + VARS(K)
   60          CONTINUE
               NGAIN = NFRNT - LFRNT
               UNNUM(NBR) = UNNUM(NBR) - 1
               PRIOR(NBR) = -W1*NGAIN - W2*DIST(NBR) - W3*UNNUM(NBR)
               IF (NFRNT.LE.0) PRIOR(NBR) = 1.0D30
            END IF
   70    CONTINUE
   80 CONTINUE
      RETURN
   90 DO 120 I = 1,NQ
         NBR = QUEUE(I)
         PNODE = PRIOR(NBR)
         J1 = I
         DO 100 K = 1,I
            J2 = J1/2
            FATHER = QUEUE(J2)
            IF (PRIOR(FATHER).GE.PNODE) GO TO 110
            QUEUE(J1) = FATHER
            NLIST(FATHER) = J1
            J1 = J2
  100    CONTINUE
  110    QUEUE(J1) = NBR
         NLIST(NBR) = J1
  120 CONTINUE
      L1 = L
      DO 200 L = L1,NODES
         NEXT = QUEUE(1)
         QNODE = QUEUE(NQ)
         PNODE = PRIOR(QNODE)
         NQ = NQ - 1
         J = 2
         J2 = 1
         IF (NQ.GT.1) QUEUE(NQ+1) = QUEUE(NQ)
         DO 130 I = 2,NQ
            IF (J.GT.NQ) GO TO 140
            IF (PRIOR(QUEUE(J)).LT.PRIOR(QUEUE(J+1))) J = J + 1
            SON = QUEUE(J)
            IF (PNODE.GE.PRIOR(SON)) GO TO 140
            QUEUE(J2) = SON
            NLIST(SON) = J2
            J2 = J
            J = J*2
  130    CONTINUE
  140    QUEUE(J2) = QNODE
         NLIST(QNODE) = J2
         LSTNUM = LSTNUM + 1
         STATUS(NEXT) = -LSTNUM
         J1 = ELTPTR(NEXT)
         J2 = ELTPTR(NEXT+1) - 1
         DO 150 J = J1,J2
            K = ELTVAR(J)
            IELIM(K) = IELIM(K) - 1
            IF (IFRNT(K).EQ.0) IFRNT(K) = LSTNUM
  150    CONTINUE
         ISTRT = ICPTR(NEXT)
         ISTOP = ICPTR(NEXT+1) - 1
         DO 190 I = ISTRT,ISTOP
            NBR = IRN(I)
            IF (STATUS(NBR).EQ.2) THEN
               NQ = NQ + 1
               QUEUE(NQ) = NBR
               STATUS(NBR) = 1
               NLIST(NBR) = NQ
            END IF
            IF (STATUS(NBR).EQ.1) THEN
               NFRNT = 0
               LFRNT = 0
               J1 = ELTPTR(NBR)
               J2 = ELTPTR(NBR+1) - 1
               DO 160 J = J1,J2
                  K = ELTVAR(J)
                  IF (IFRNT(K).EQ.0) NFRNT = NFRNT + VARS(K)
                  IF (IELIM(K).EQ.1) LFRNT = LFRNT + VARS(K)
  160          CONTINUE
               NGAIN = NFRNT - LFRNT
               UNNUM(NBR) = UNNUM(NBR) - 1
               PRIOR(NBR) = -W1*NGAIN - W2*DIST(NBR) - W3*UNNUM(NBR)
               IF (NFRNT.LE.0) PRIOR(NBR) = 1.0D30
               PNODE = PRIOR(NBR)
               J = NLIST(NBR)
               DO 170 K = 1,NQ
                  IF (J.EQ.1) GO TO 180
                  J2 = J/2
                  FATHER = QUEUE(J2)
                  IF (PRIOR(FATHER).GE.PNODE) GO TO 180
                  QUEUE(J) = FATHER
                  NLIST(FATHER) = J
                  J = J2
  170          CONTINUE
  180          QUEUE(J) = NBR
               NLIST(NBR) = J
            END IF
  190    CONTINUE
  200 CONTINUE
      RETURN
      END
C***********************************************************
      SUBROUTINE MC63FD(NSUP,NELT,NZ,SELTVR,SELTPR,IPERM,KEY,LIST)
C***********************************************************
      INTEGER NELT,NSUP,NZ
      INTEGER IPERM(NSUP),KEY(NELT),LIST(NELT),SELTPR(NELT+1),SELTVR(NZ)
      INTEGER I,IELL,J,JSTOP,JSTRT,L
      DO 10 I = 1,NELT
         LIST(I) = I
         KEY(I) = 0
   10 CONTINUE
      DO 30 IELL = 1,NELT
         JSTRT = SELTPR(IELL)
         JSTOP = SELTPR(IELL+1) - 1
         J = SELTVR(JSTRT)
         KEY(IELL) = IPERM(J)
         DO 20 I = JSTRT + 1,JSTOP
            J = SELTVR(I)
            L = IPERM(J)
            IF (L.LT.KEY(IELL)) KEY(IELL) = L
   20    CONTINUE
   30 CONTINUE
      RETURN
      END
C********************************************************************
      SUBROUTINE MC63GD(N,LIST,KEY)
C***********************************************************************
      INTEGER N
      INTEGER KEY(N),LIST(N)
      INTEGER GUESS,I,LL,LM,LR,LTEMP,NL,NR,STKTOP
      INTEGER IW(64)
      DO 10 I = 1,64
         IW(64) = 0
   10 CONTINUE
      LL = 1
      LR = N
      STKTOP = 0
   20 CONTINUE
      IF (LL.LT.LR) THEN
         NL = LL
         NR = LR
         LM = (LL+LR)/2
         GUESS = KEY(LIST(LM))
   30    CONTINUE
         IF (KEY(LIST(NL)).LT.GUESS) THEN
            NL = NL + 1
            GO TO 30
         END IF
   40    CONTINUE
         IF (GUESS.LT.KEY(LIST(NR))) THEN
            NR = NR - 1
            GO TO 40
         END IF
         IF (NL.LT. (NR-1)) THEN
            LTEMP = LIST(NL)
            LIST(NL) = LIST(NR)
            LIST(NR) = LTEMP
            NL = NL + 1
            NR = NR - 1
            GO TO 30
         END IF
         IF (NL.LE.NR) THEN
            IF (NL.LT.NR) THEN
               LTEMP = LIST(NL)
               LIST(NL) = LIST(NR)
               LIST(NR) = LTEMP
            END IF
            NL = NL + 1
            NR = NR - 1
         END IF
         STKTOP = STKTOP + 1
         IF (NR.LT.LM) THEN
            IW(STKTOP) = NL
            IW(32+STKTOP) = LR
            LR = NR
         ELSE
            IW(STKTOP) = LL
            IW(32+STKTOP) = NR
            LL = NL
         END IF
         GO TO 20
      END IF
      IF (STKTOP.NE.0) THEN
         LL = IW(STKTOP)
         LR = IW(32+STKTOP)
         STKTOP = STKTOP - 1
         GO TO 20
      END IF
      RETURN
      END





