!***************************************************************************
!====      read a DCD file and perform analysis of collective motion      ====
!====        Wensheng Xu, CIAC, May 2025                                  ====
!****************************************************************************

!====SetCOM: define global variables and arrays
MODULE SetCOM

	INTEGER, PARAMETER:: T0Max=100
	
	INTEGER NPart, NBead, NPoly,NBond
	INTEGER IT0, TMax
	
	REAL, ALLOCATABLE:: AX(:), AY(:), AZ(:)
	
	DOUBLE PRECISION, PARAMETER:: PI=4.0*ATAN(1.0)
	
	DOUBLE PRECISION Box, Vol, Rho, Delt
	DOUBLE PRECISION, ALLOCATABLE:: RX(:), RY(:), RZ(:)
	DOUBLE PRECISION, ALLOCATABLE:: PX(:), PY(:), PZ(:)
	CHARACTER*1 FType
	
END MODULE

!====main program
PROGRAM MAIN
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER NFrame, I, K, Foutput
	
	DOUBLE PRECISION TimeStart, TimeFinal
!=================================================================
	OPEN (UNIT=12, FILE='LogDynL')
	DO K = 1, 3
		CALL CPU_TIME(TimeStart)
		
		IF (K==1) FType = 'L'
		IF (K==2) FType = 'M'
		IF (K==3) FType = 'S'
	
!		====open DCD file
		OPEN (UNIT=30, FILE='../Conf/Coord'//FType//'.dcd', STATUS='old', FORM='unformatted')
	
!		====set basic information
		CALL INIT(NFrame)
	
		Foutput = INT(ANINT(DBLE(NFrame)/10))

!		====start of calculations
		CALL Get_Dynamic(0)
	
		WRITE (12, "(/, '---->>calculation begins')")
		DO I = 1, NFrame
!			====get confiruration
			CALL Get_Conf
			IF (MOD(I, 10)==0) WRITE (*, *) I, NFrame
!			====dynamic properties
			CALL Get_Dynamic(1)
		
!			====write intermediate results
			IF (MOD(I, Foutput)==0) THEN
				CALL Get_Dynamic(2)
				WRITE (12, "('  #frames completed: ', I8)") I
			END IF
		END DO
		CALL Get_Dynamic(2)
		CALL Get_Dynamic(3)
		WRITE (12, "('<<----calculation completed', /)")
	
!		====free memory
		CALL FreeMemory
	
		CALL CPU_TIME(TimeFinal)
		WRITE (12, "('Total time is ', F8.2, ' minutes')") &
			(TimeFinal - TimeStart)/60.0
	
		CLOSE (30)
	END DO
	CLOSE (12)

END

!====INIT: read input parameters and set initial information
SUBROUTINE INIT(NFrame)
	
	USE SetCOM
	IMPLICIT NONE
	
	INTEGER NFrame, Natom

	DOUBLE PRECISION DeltL, TProdL, DeltM, TProdM, DeltS, TProdS
!================================================================
!	====read system information
	OPEN (UNIT=15, FILE='../Conf/Basic.txt')
	READ (15, *)
	READ (15, *) NPart, NBead, NPoly
	READ (15, *)
	READ (15, *) Box
	READ (15, *)
	READ (15, *) DeltL, TProdL, DeltM, TProdM, DeltS, TProdS
	CLOSE (15)
	
	Vol = Box**3
	Rho = NPart/Vol
	NBond = NPoly*(NBead-1)
	
	IF (FType=='L') THEN
		Delt = DeltL
		IT0 = INT(ANINT(0.5*TProdL/(T0Max*Delt)))
	END IF
	IF (FType=='M') THEN
		Delt = DeltM
		IT0 = INT(ANINT(0.5*TProdM/(T0Max*Delt)))
	END IF
	IF (FType=='S') THEN
		Delt = DeltS
		IT0 = INT(ANINT(0.5*TProdS/(T0Max*Delt)))
	END IF
	
!	====allocate arrays
	CALL SetMemory
	
!	====print basic information
	WRITE (12, "('---->>Basic information')")
	WRITE (12, "('  Total particles        : ', I10)") NPart
	WRITE (12, "('  Chain length           : ', I10)") NBead
	WRITE (12, "('  Chain number           : ', I10)") NPoly
	WRITE (12, "('  Box dimension          : ', F10.4)") Box
	WRITE (12, "('  Number density         : ', F10.4, /)") Rho
	WRITE (12, "('  Time interval for TCF  : ', F10.3)") Delt
	WRITE (12, "('  Maximum time  for TCF  : ', F10.2, /)") Delt*IT0*T0Max

!	====read header in the DCD file
	CALL ReadHeader(NFrame, Natom)
	
	IF (Natom .NE. NPart) STOP 'Natom from DCD file is inconsistent with NPart!'
	
	WRITE (12, "('  Number of frames in dcd: ', I10)") NFrame
	WRITE (12, "('  Number of atoms  in dcd: ', I10, /)") Natom
	
	TMax = NFrame
	
RETURN
END

!====Get_Dynamic: call subroutines for computing dynamic properties
SUBROUTINE Get_Dynamic(Switch)

	IMPLICIT NONE
	
	INTEGER Switch
!=================================================================

	CALL Get_ASL(Switch)

RETURN
END

!====Subroutines
INCLUDE "SetConf.f90"
INCLUDE "DynaASL.f90"
