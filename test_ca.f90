PROGRAM test_ca
IMPLICIT NONE

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:,:)::COMPOSITON,Melt_CAO,MELT_SIO2,MELT_FRACTION
INTEGER::timesteps, T,I,J,K,BASE_X,BASE_Y,BASE_Z
CHARACTER::DUMMY_CHAR



OPEN(800,FILE='MELT_CAO',form='unformatted')
OPEN(801,FILE='MELT_SIO2',form='unformatted')
OPEN(802,FILE='MELT_FRACTION',form='unformatted')


OPEN(1000,FILE='INPUT')

OPEN(1001,FILE='ZERO_CA')


! Allocate
REWIND(1000)

READ(1000,400) DUMMY_CHAR
READ(1000,200) BASE_X
READ(1000,400) DUMMY_CHAR
READ(1000,200) BASE_Y
READ(1000,400) DUMMY_CHAR
READ(1000,200) BASE_Z

timesteps=200


ALLOCATE(MELT_CAO(BASE_X,BASE_Y,BASE_Z,timesteps))
ALLOCATE(MELT_SIO2(BASE_X,BASE_Y,BASE_Z,timesteps))
ALLOCATE(MELT_FRACTION(BASE_X,BASE_Y,BASE_Z,timesteps))


DO T=1,timesteps
 READ(800) MELT_CAO(:,:,:,T)
 READ(801) MELT_SIO2(:,:,:,T)
 READ(802) MELT_FRACTION(:,:,:,T)
write(1001,*) 'TIME',T

DO I=1,BASE_X
 DO J=1,BASE_Y
  DO K=1,BASE_Z
   IF (MELT_FRACTION(I,J,K,T)>0.0001) THEN
     IF (MELT_CAO(I,J,K,T)<.0001) THEN
      write(1001,100) I,J,K,T, MELT_CAO(I,J,K,T),MELT_SIO2(I,J,K,T),MELT_FRACTION(I,J,K,T)
     END IF
   END IF
 END DO
END DO
END DO
END DO

100 FORMAT(7F20.10)
200 FORMAT(I8)
300 FORMAT(1F40.15)
301 FORMAT(4F40.15)
302 FORMAT(4F40.20)
400 FORMAT(A30)
END PROGRAM









