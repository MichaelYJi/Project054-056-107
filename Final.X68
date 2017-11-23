*-----------------------------------------------------------
* Title      :
* Written by :
* Date       :
* Description:
*-----------------------------------------------------------
    ORG    $1000
START:          

    LEA filename,A1
    MOVE.B #51,D0
    TRAP #15
    
    MOVE.L D3,D2
    MOVE.B #55,D0
    TRAP #15
    
    LEA buffer,A1
    MOVE.W #datasize,D2
    MOVE.B #53,D0
    TRAP #15
    MOVEA.L A1,A3
    MOVE.B D0,D7
     
loop1 MOVE.B (A1),D3
      CMP #$FF,D3
      BEQ done
      JSR convert
      ADDQ #1,A1
      BRA loop1

done CLR D3
    MOVE.B (A3)+,D3
    
loop2
    JSR checkop
    CMP #$CC,D5 *CC represents invalid opcode
    BEQ clr
    JSR printdata 
    JSR print
clr CMP #$DD,D3 *DD represents new line
    BNE continue
    MOVE.B (A3)+,D3
    BRA clr
continue CMP #$FF,D3 *FF represents end of stream
    BNE loop2
        
    SIMHALT

* Put variables and constants here
filename dc.b 'test1.txt', 0
datasize equ 2000
buffer ds.b 2000
data dc.b 'DATA ',0
source dc.b 0
dest dc.b 0
sourcem dc.b 0
destm dc.b 0

pmove_b dc.b 'MOVE.B ',0

dregister dc.b 'D',0

comma dc.b ',',0

convert CMP #$FF,D3
       BEQ exit
       CMP #$30,D3
       BLT fake
       CMP #$46,D3
       BGT fake
       CMP #$39,D3
       BLE not
       CMP #$41,D3
       BGE not
       BRA fake
*       CMP #$0D,D5
*       BEQ fake
*       CMP #$0A,D5
*       BEQ fake
*       CMP #$A0,D5
*       BEQ fake
*       CMP #$F0,D5
*       BEQ fake
*       BRA not
fake   MOVE.B #$DD,D3
       BRA exit
not    CMP #$40,D3
       BGT let
       SUBI.L #$30,D3
       BRA exit
let    SUBI.L #$37,D3
exit   MOVE.B D3,(A1)
       RTS
       
newline MOVE #$0D,D1
        MOVE #6,D0
        TRAP #15
        MOVE #$0A,D1
        MOVE #6,D0
        TRAP #15
        RTS
        
print 
loop    CMP #$DD,D3
        BEQ done1
        CMP #$FF,D3
        BEQ done2
        MOVE.B D3,D1
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
        MOVE.B (A3)+,D3
        BRA loop
done1   JSR newline
done2   RTS

        
printdata LEA data,A1
        MOVE.B #14,D0
        TRAP #15
        RTS 
        
checkop MOVEA.W A3,A4

        MOVE.B -(A4),D4
        
        CMP #$01,D4
        BEQ MOVE_B
        
        *First level jump table goes here
        RTS
        
MOVE_B  CLR D5 *Crashes program if checks aren't in
        CLR D6
        CLR D7
        
        *Isolate and get source mode
        ADDQ #2,A4
        MOVE.B (A4),D5
        ADDQ #1,A4
        MOVE.B (A4),D6
        MOVE.B #30,D7
        ASL.L D7,D5
        ASR.L D7,D5
        ASR.B #$3,D6
        ADD.B D6,D5
        MOVE.B D5,sourcem
        
        *Isolate and get source register
        MOVEA.W A3,A4
        ADDQ #3,A4
        MOVE.B (A4),D5
        MOVE.B #29,D7
        ASL.L D7,D5
        ASR.L D7,D5
        MOVE.B D5,source
        
        *Isolate and get destination mode
        MOVEA.W A3,A4
        ADDQ #1,A4
        MOVE.B (A4),D5
        ADDQ #1,A4
        MOVE.B (A4),D6
        MOVE.B #31,D7
        ASL.L D7,D5
        SUBI.B #2,D7
        ASR.L D7,D5
        ASR.B #$2,D6
        ADD.B D5,D6
        MOVE.B D6,destm
        
        *Isolate and get destination register
        MOVEA.W A3,A4
        ADDQ #1, A4
        MOVE.B (A4),D5
        ASR.B #$1,D5
        MOVE.B D5,dest
        
        *Check valid opcode
        
        *Print
        LEA pmove_b,A1
        MOVE.B #14,D0
        TRAP #15
        MOVEA.W EA,A4
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR 0(A4,D6)
        LEA comma,A1
        MOVE.B #14,D0
        MOVE.B (dest),D5
        MOVE.B (sourcem),D6
        JSR 0(A4,D6)
        
        CLR D5
        MOVE.B #$CC,D5
        ADDQ #4,A3
mbexit  RTS
                
EA      JMP EA_000
        JMP EA_001
        JMP EA_010
        JMP EA_011
        JMP EA_100
        JMP EA_101
        JMP EA_110
        JMP EA_111

EA_000  LEA dregister,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B D5,D1
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
        RTS
        
EA_001  NOP
        RTS

EA_010  NOP
        RTS

EA_011  NOP
        RTS 

EA_100  NOP
        RTS    

EA_101  NOP
        RTS  
        
EA_110  NOP
        RTS

EA_111  NOP
        RTS        
    END    START        ; last line of source











*~Font name~Courier New~
*~Font size~10~
*~Tab type~1~
*~Tab size~4~
