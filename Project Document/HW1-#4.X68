*-----------------------------------------------------------
* Title      :Represent a negative integer with two's complement
* Written by :Dawit Abera
* Date       :Oct 12/2017
* Description:Two's complement 
*-----------------------------------------------------------
CR           EQU     $0D                  ASCII code for Carriage Return
LF           EQU     $0A                  ASCII code for Line Feed
FirstInput   EQU     $FF9A                -102 in HEX
SecondInput  EQU     $FFA9                -87 in HEX
StoreResult  EQU     $5000                Address location for the result output
*-----------------------------------------------------------

START:       ORG         $1000
ADD         MOVE.W      #FirstInput,D0               Load D0 with first input
            ADDI.L      #SecondInput,D0              Load the second input 
            MOVEA.W     #StoreResult,A0              Load the address regi for the result
            BVS         ERROR                        Branch when overflow occurs
            MOVE.W      D0,(A0)                      Load D0 and D1 
            MOVE.W      (A0),D1                      Move the data from A0 to D1
            
            MOVE.B      #3,D0                        Show the value in D1
            TRAP        #15
            MOVE.B      #9,D0
            TRAP        #15
            
            
            
ERROR       LEA          MESSAGE,A1                  Loads MESSAGE into address registor A1
            MOVE.B       #14,D0
            TRAP         #15
        
MESSAGE     DC.B         'Their is an overflow',0    Display MESSAGE if overflow occoures 
            END          START


*~Font name~Courier New~
*~Font size~10~
*~Tab type~1~
*~Tab size~4~
