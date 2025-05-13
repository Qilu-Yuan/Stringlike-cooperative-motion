!====ReadHeader: read header in the DCD file
SUBROUTINE ReadHeader(NFrame, NAtom)
	
	USE SetCOM
	IMPLICIT NONE
	
	INTEGER NFrame, NAtom, I, dummyi
	
	REAL dummyr
	
	CHARACTER*4 dummyc
!================================================================
	READ (30) dummyc, NFrame, (dummyi, I=1, 8), dummyc, (dummyi, I=1, 9)
	READ (30) dummyi, dummyr
	READ (30) NAtom
	
RETURN
END

!====Get_Conf: get configuration
SUBROUTINE Get_Conf
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER I,J, P, B
	
	DOUBLE PRECISION D
	DOUBLE PRECISION L2
!================================================================
!	====read particle coordinates
	READ (30) (D, I = 1, 6)
	READ (30) (AX(I), I = 1, NPart)
	READ (30) (AY(I), I = 1, NPart)
	READ (30) (AZ(I), I = 1, NPart)
		
	DO I = 1, NPart
		RX(I) = DBLE(AX(I))
		RY(I) = DBLE(AY(I))
		RZ(I) = DBLE(AZ(I))
	END DO

!	====get chain centers of mass
	DO P = 1, NPoly
		PX(P) = 0
		PY(P) = 0
		PZ(P) = 0
		DO B = 1, NBead
			I = (P - 1)*NBead + B
			PX(P) = PX(P) + RX(I)
			PY(P) = PY(P) + RY(I)
			PZ(P) = PZ(P) + RZ(I)
		END DO
		PX(P) = PX(P)/NBead
		PY(P) = PY(P)/NBead
		PZ(P) = PZ(P)/NBead
	END DO
	
RETURN
END

!====SetMemory: allocate memory for all variables
SUBROUTINE SetMemory
	
	USE SetCOM
	IMPLICIT NONE
	
	ALLOCATE(AX(NPart), AY(NPart), AZ(NPart))
	ALLOCATE(RX(NPart), RY(NPart), RZ(NPart))
	ALLOCATE(PX(NPoly), PY(NPoly), PZ(NPoly))
RETURN
END

!====FreeMemory: free memory for all variables
SUBROUTINE FreeMemory
	
	USE SetCOM
	IMPLICIT NONE
	
	DEALLOCATE(AX, AY, AZ)
	DEALLOCATE(RX, RY, RZ)
	DEALLOCATE(PX, PY, PZ)
RETURN
END
