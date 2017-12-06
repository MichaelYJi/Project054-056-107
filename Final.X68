*-----------------------------------------------------------
* Title      :
* Written by :
* Date       :
* Description:
*-----------------------------------------------------------
    ORG    $1000
START:          
start
    LEA filename,A1
    MOVE.B #51,D0
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
    MOVEA.L #0,A6
    MOVE.B #11,D0
    MOVE.W #$FF00,D1
    TRAP #15
innerloop
    ADDQ #1,A6
    MOVE.L A3,D1
    MOVE.B #16,D2
    MOVE.B #15,D0
    TRAP #15
    MOVE.L #$0,D1
    LEA tab,A1
    MOVE.B #14,D0
    TRAP #15
    CLR D5
    CLR D6
    MOVE.B #$0,badop
    JSR checkop
    MOVE.B badop,D6
    CMP #$EE,D6
    BNE inlpn1
    LEA invalid,A1
    MOVE.B #14,D0
    TRAP #15
    JSR newline
    BRA clr
inlpn1 CMP #$CC,D5 *CC represents invalid opcode
    BNE notvalid
    JSR newline
    BRA clr
notvalid JSR printdata 
    JSR print
clr CMP #$DD,D3 *DD represents new line
    BNE continue
    MOVE.B (A3)+,D3
    BRA clr
continue CMP #$FF,D3 *FF represents end of stream
    BEQ complete
    CMPA #29,A6
    BNE innerloop
requery LEA nextpage,A1
    MOVE.B #18,D0
    TRAP #15
    CMP #1,D1
    BEQ loop2
    BRA requery
    
complete MOVE.B #56,D0
    TRAP #15
    
reprompt LEA restart,A1
    MOVE.B #18,D0
    TRAP #15
    CMP #1,D1
    BNE nextquest
    MOVE.B #11,D0
    MOVE.W #$FF00,D1
    TRAP #15
    BRA start
nextquest CMP #2,D1
    BNE reprompt
    MOVE.B #9,D0
    TRAP #15
    SIMHALT

* Put variables and constants here
filename dc.b 'test1.txt', 0
datasize equ 10000
buffer ds.b 10000
data dc.b 'DATA ',0
invalid dc.b HT,'INVALID EA',0
source dc.b 0
dest dc.b 0
sourcem dc.b 0
destm dc.b 0
counter dc.b 0
badop dc.b 0
size dc.b 0
dir dc.b 0
mparam dc.b 0
rparam dc.b 0
dirparam dc.b 0

CR EQU $0D
LF EQU $0A
HT EQU $09
restart dc.b 'Press 1 to restart or 2 to exit...',CR,LF,0
nextpage dc.b 'Press 1 to continue simulation...',CR,LF,0
tab dc.b '',HT,0

pmove_b dc.b 'MOVE.B ',0
pmove_w dc.b 'MOVE.W ',0
pmove_l dc.b 'MOVE.L ',0
pmovea_w dc.b 'MOVEA.W ',0
pmovea_l dc.b 'MOVEA.L ',0
pmoveq dc.b 'MOVEQ ',0
pmovem dc.b 'MOVEM',0
pnop dc.b 'NOP ',0
padd dc.b 'ADD',0
padda dc.b 'ADDA',0
paddi dc.b 'ADDI',0
paddq dc.b 'ADDQ',0
pmuls dc.b 'MULS.W ',0
plea dc.b 'LEA ',0
pdivu dc.b 'DIVU.W ',0
psub dc.b 'SUB',0
pclr dc.b 'CLR',0

dregister dc.b 'D',0
aregister dc.b 'A',0
hashtag dc.b '#',0
dollarsign dc.b '$',0
comma dc.b ',',0
dotb dc.b '.B ',0
dotw dc.b '.W ',0
dotl dc.b '.L ',0
leftp dc.b '(',0
rightp dc.b ')',0
plus dc.b   '+',0
minus dc.b '-',0
slash dc.b '/',0

pas dc.b 'AS',0
pls dc.b 'LS',0
pro dc.b 'RO',0
pright dc.b 'R',0
pleft dc.b 'L',0
space dc.b ' ',0

pjsr dc.b 'JSR ',0
prts dc.b 'RTS',0

pbcc dc.b 'BCC ',0
pbgt  dc.b 'BGT ',0
pble dc.b 'BLE ',0

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

*First level jump table        
checkop MOVEA.W A3,A4

        MOVE.B -(A4),D4
        
        CMP #$00,D4
        BEQ ADDI_
        
        CMP #$01,D4
        BEQ MOVE_B
        
        CMP #$02,D4
        BEQ MOVE_L
        
        CMP #$03,D4
        BEQ MOVE_W
        
        CMP #$04,D4
        BEQ code0100
        
        CMP #$05,D4
        BEQ ADDQ_
        
        CMP #$06,D4
        BEQ BXX
        
        CMP #$07,D4
        BEQ MOVEQ_
        
        CMP #$08,D4
        BEQ DIVU_
        
        CMP #$09,D4
        BEQ SUB_
        
        CMP #$0C,D4
        BEQ code1100
        
        CMP #$0D,D4
        BEQ ADD_
        
        CMP #$0E,D4
        BEQ rl_check

        RTS
        
code0100 JSR IsolateRegisters

        *CLR check
        MOVE.B (A4),D4
        CMP #$2,D4
        BEQ CLR_
        
        *LEA check
        MOVE.B (destm),D4
        CMP #$7,D4
        BEQ LEA_
        
        *JSR check
        CLR D5
        MOVEA A3,A4
        MOVE.B -(A4),D4 *resets to 15->12
        MOVE.B (destm),D5
        CMP #$2,D5
        BEQ JSR_

        *NOP/RTS check
        MOVEA.W A3,A4
        ADDQ #1,A4
        MOVE.B (A4)+,D4
        CMP #$07,D4
        BNE next0100
        MOVE.B (A4),D4
        CMP #$01,D4
        BEQ NOP_
        CMP #$05,D4
        BEQ RTS_   
        
        *Else always check MOVEM
next0100 BRA MOVEM_
        
        RTS
        
code1100
        *MULS check
        JSR IsolateRegisters
        MOVE.B (destm),D4
        CMP #$7,D4
        BEQ MULS_
        
        *Branch to AND otherwise

        RTS
        
CLR_    CLR D5
        JSR IsolateRegisters 
        
        MOVE.B (sourcem),D5
        CMP #$1,D5
        BEQ clrexit
        CMP #$7,D5
        BNE _GET_SIZE
        MOVE.B (source),D6
        CMP #$4,D6
        BEQ clrexit
      
_GET_SIZE
        LEA pclr,A1
        MOVE.B #14,D0
        TRAP #15
         
        LSR.B       #6, D4              * SSXXXXXX -> 000000SS WHERE SS :-> SIZE
        CMP.B       #0, D4              * IF SIZE IS 0 (BYTE)
        BEQ         CLR_BYTE            * SIZE IS O (BYTE) BRANCH TO CLR_BYTE
        CMP.B       #1, D4              * IF SIZE IS 1 (WORD)
        BEQ         CLR_WORD            * SIZE IS 1 (WORD) BRANCH TO CLR_WORD
        BRA         CLR_LONG            * SIZE IS 2 (LONG) GO TO CLR_LONG
        
CLR_BYTE
        LEA dotb,A1
        MOVE.B #14,D0
        TRAP #15
        BRA         CLR_MODE            * GO TO CLR_MODE   
CLR_WORD
        LEA dotw,A1
        MOVE.B #14,D0
        TRAP #15
        BRA         CLR_MODE            * GO TO CLR_MODE    
CLR_LONG
        LEA dotl,A1
        MOVE.B #14,D0
        TRAP #15
       BRA         CLR_MODE
CLR_MODE       
        MOVE.B sourcem,D6
        MOVE.B source,D5
        JSR EA

        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
clrexit RTS

ADDI_   JSR IsolateRegisters

        *Check valid opcode
        MOVE.B (sourcem),D5
        CMP #$1,D5
        BEQ addiexit
        MOVE.B (sourcem),D5
        CMP #$7,D5
        BNE addin0
        MOVE.B (source),D5
        CMP #$4,D5
        BEQ addiexit
addin0  MOVEA.W A3,A4
        MOVE.B (A4)+,D5
        CMP #$6,D5
        BNE addiexit
        
        *Isolate size
        MOVE.B (A4),D5
        MOVE.B #28,D7
        LSL.L D7,D5
        ADDQ #2,D7
        LSR.L D7,D5
        
        *Check size
        CMP #$0,D5
        BNE addin1
        LEA dotb,A5
        MOVE.B #2,D4
        BRA prnaddi
addin1  CMP #$1,D5
        BNE addin2
        LEA dotw,A5
        MOVE.B #4,D4
        BRA prnaddi
addin2  CMP #$2,D5
        BNE addiexit
        LEA dotl,A5
        MOVE.B #8,D4
        
        *Print
prnaddi LEA paddi,A1
        MOVE.B #14,D0
        TRAP #15
        LEA (A5),A1
        TRAP #15
        MOVE.B #$7,D6
        MOVE.B #$4,D5
        JSR EA
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR EA

        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        *JSR newline
addiexit RTS
        
MOVE_B  JSR IsolateRegisters
                
        *Check valid opcode
        MOVE.B #$1,D7    
        MOVE.B (destm),D6
        CMP D6,D7
        BEQ mbexit
        
        *Print
        LEA pmove_b,A1
        MOVE.B #14,D0
        TRAP #15
        CLR D5
        CLR D6
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR EA
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (dest),D5
        MOVE.B (destm),D6
        JSR EA
        
        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        *JSR newline
mbexit  RTS

MOVE_W  JSR IsolateRegisters
        
        *Check movea or move
        MOVE.B #$1,D7
        MOVE.B (destm),D6
        CMP D6,D7
        BNE printnaw
        LEA pmovea_w,A1
        MOVE.B #14,D0
        TRAP #15
        BRA printmw
        
        *Print
printnaw LEA pmove_w,A1
        MOVE.B #14,D0
        TRAP #15
printmw CLR D5
        CLR D6
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR EA
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (dest),D5
        MOVE.B (destm),D6
        JSR EA
        
        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        *JSR newline
mwexit  RTS

MOVE_L  JSR IsolateRegisters 
        
        *Check movea or move
        MOVE.B #$1,D7
        MOVE.B (destm),D6
        CMP D6,D7
        BNE printnal
        LEA pmovea_l,A1
        MOVE.B #14,D0
        TRAP #15
        BRA printml
        
        *Print
printnal LEA pmove_l,A1
        MOVE.B #14,D0
        TRAP #15
printml CLR D5
        CLR D6
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR EA
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (dest),D5
        MOVE.B (destm),D6
        JSR EA
        
        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        *JSR newline
mlexit  RTS

ADDQ_   JSR IsolateRegisters
        
        *Check valid source mode
        MOVE.B (sourcem),D5
        CMP #$7,D5
        BNE addqn0
        MOVE.B (source),D5
        CMP #$4,D5
        BEQ addqexit

        *Check 8th bit for 0
addqn0  MOVEA.W A3,A4
        MOVE.B (A4)+,D5
        MOVE.L #32,D7
        ASL.L D7,D5
        ASR.L D7,D5
        CMP #$1,D5
        BEQ addqexit

        *Isolate size
        MOVE.B (A4),D5
        MOVE.B #28,D7
        LSL.L D7,D5
        ADDQ #2,D7
        LSR.L D7,D5
        
        *Check size
        CMP #$0,D5
        BNE addq1
        LEA dotb,A5
        MOVE.B #2,D4
        BRA prnaddq
addq1  CMP #$1,D5
        BNE addq2
        LEA dotw,A5
        MOVE.B #4,D4
        BRA prnaddq
addq2  CMP #$2,D5
        BNE addqexit
        LEA dotl,A5
        MOVE.B #8,D4
        
        *Print
prnaddq LEA paddq,A1
        MOVE.B #14,D0
        TRAP #15
        LEA (A5),A1
        TRAP #15
        LEA hashtag,A1
        TRAP #15
        LEA dollarsign,A1
        TRAP #15
        MOVE.B dest,D1
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR EA

        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        *JSR newline
addqexit RTS

BXX
        CLR D4
        ADDQ #1,A4  *STORE 11->8 
        MOVE.B (A4),D4
        
        CMP #$4,D4
        BEQ BCC_
        
        CMP #$E,D4
        BEQ BGT_

        CMP #$F,D4
        BEQ BLE_

        BRA notBXX
        
BCC_
        LEA pbcc,A1
        MOVE.B #14,D0
        TRAP #15
        BRA BXX_2

BGT_
        LEA pbgt,A1
        MOVE.B #14,D0
        TRAP #15
        BRA BXX_2

BLE_
        LEA pble,A1
        MOVE.B #14,D0
        TRAP #15
        BRA BXX_2

BXX_2
        ADDQ #1,A4  *STORE 7->4
        MOVE.B (A4),D4
        CLR.L D2
        MOVE.B D4,D2
        LSL.B #4,D2
        
        CLR.L D4
        ADDQ #1,A4  *STORE 3->0
        MOVE.B (A4),D4
        ADD.B D2,D4
        NOT.B D4    *FLIP BITS
        ADD.B #$1,D4 *NOW 2 COMPLEMENT
        
        CMP #127,D4
        BGT BXX_NEG
        
        LEA plus,A1
        MOVE.B #14,D0
        TRAP #15
        
        MOVE.B D4,D1
        MOVE.B #3,D0
        TRAP #15
        
        MOVE.B #$CC,D5  *DONE
	    ADDQ #4,A3
	    MOVE.B (A3),D3
	    RTS
        
BXX_NEG    
        LEA minus,A1
        MOVE.B #14,D0
        TRAP #15

        CLR.L D2
        MOVE.W #256,D2
        SUB.L D2,D4
        
        MOVE.B D4,D1
        MOVE.B #3,D0
        TRAP #15
        
        MOVE.B #$CC,D5  *DONE
	    ADDQ #4,A3
	    MOVE.B (A3),D3
	    RTS

notBXX  RTS

MOVEQ_  CLR D5
        CLR D6
        CLR D7

        *Check 8th bit for 0
        ADDQ #1,A4
        MOVE.B (A4),D5
        MOVE.L #31,D7
        ASL.L D7,D5
        ASR.L D7,D5
        CMP #$1,D5
        BEQ mqexit
        
        *Get data register
        MOVE.B (A4)+,D5
        MOVE.B #28,D7
        LSL.L D7,D5
        ADDQ #1,D7
        LSR.L D7,D5
        MOVE.B D5,dest
        
        *Get and print value
        LEA pmoveq,A1
        MOVE.B #14,D0
        TRAP #15
        LEA hashtag,A1
        TRAP #15
        LEA dollarsign,A1
        TRAP #15
        MOVE.B (A4)+,D1
        CMP #$0,D1
        BEQ mqpnext
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
mqpnext MOVE.B (A4),D1
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        *Print data register
        CLR D5
        MOVE.B #$0,D6
        MOVE.B (dest),D5
        JSR EA
        
        CLR D5
        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
mqexit  RTS

MOVEM_  MOVE.B (destm),D5
        CMP #$2,D5
        BEQ mvmn1
        CMP #$3,D5
        BNE mvmexit
        
        *More checks to confirm valid opcode
mvmn1   MOVE.B (dest),D5
        CMP #$4,D5
        BEQ mvmn2
        CMP #$6,D5
        BNE mvmexit
        
        *Isolate size
mvmn2   MOVE.B (destm),D5
        MOVE.B #$31,D7
        LSL.L D7,D5
        LSR.L D7,D5
        MOVE.B D5,size
        
        *Store size letter
        CMP #$0,D5
        BNE mvmns1
        MOVEA.W dotw,A5
        BRA mvmns2
mvmns1  MOVEA.W dotl,A5
        
        *Isolate direction
mvmns2  MOVE.B (dest),D5
        MOVE.B #$30,D7
        LSL.L D7,D5
        ADDQ #1,D7
        LSR.L D7,D5
        MOVE.B D5,dir
        
        *Checks to confirm valid source EA
        MOVE.B (sourcem),D5
        CMP #$2,D5
        BEQ mvmpsrt
        CMP #$4,D5
        BNE mvmn3
        MOVE.B (dir),D6
        CMP #$0,D6
        BNE mvmexit
        BRA mvmpsrt
mvmn3   CMP #$3,D5
        BNE mvmn4
        MOVE.B (dir),D6
        CMP #$1,D6
        BNE mvmexit
        BRA mvmpsrt
mvmn4   CMP #$7,D5
        BNE mvmexit
        
        *Print start
mvmpsrt LEA pmovem,A1
        MOVE.B #14,D0
        TRAP #15
*        LEA (A5),A1
*        MOVE.B #14,D0
*        TRAP #15
        
        *If ea->list, print ea first
        MOVE.B (dir),D5
        CMP #$1,D5
        BNE mvmn5
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR EA
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        
        MOVEA.W A3,A4
        ADDQ #3,A4
        MOVE.B #$0,rparam
        MOVE.B #$0,mparam
        MOVE.B #$0,dirparam
        JSR D0_D3
*        LEA slash,A1
*        MOVE.B #14,D0
*        TRAP #15
        ADDQ #1,A4
        MOVE.B #$4,rparam
        JSR D0_D3
*        LEA slash,A1
*        MOVE.B #14,D0
*        TRAP #15
        ADDQ #1,A4
        MOVE.B #$0,rparam
        MOVE.B #$1,mparam
        JSR D0_D3
*        LEA slash,A1
*        MOVE.B #14,D0
*        TRAP #15
        ADDQ #1,A4
        MOVE.B #$4,rparam
        JSR D0_D3
        
        BRA mvmdone
        
        *If list->ea, print list first
mvmn5   MOVEA.W A3,A4
        ADDQ #6,A4
        MOVE.B #$7,rparam
        MOVE.B #$1,mparam
        MOVE.B #$1,dirparam
        JSR D0_D3
*        LEA slash,A1
*        MOVE.B #14,D0
*        TRAP #15
        SUBQ #1,A4
        MOVE.B #$3,rparam
        JSR D0_D3
*        LEA slash,A1
*        MOVE.B #14,D0
*        TRAP #15
        SUBQ #1,A4
        MOVE.B #$7,rparam
        MOVE.B #$0,mparam
        JSR D0_D3
*        LEA slash,A1
*        MOVE.B #14,D0
*        TRAP #15
        SUBQ #1,A4
        MOVE.B #$3,rparam
        JSR D0_D3
        
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR EA

*        *Loop through registers D0-A7
*mvmn5   MOVEA.W A3,A4
*        ADDQ #3,A4
*        CLR D2
*mvmlp1  MOVE.B (A4),D5
*        *Shift left
*        *Add to shift
*        *Shift right
*        *
        

mvmdone MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
mvmexit RTS

ApplyDirection
        MOVE.B rparam,D6
        MOVE.B (dirparam),D5
        CMP #$0,D5
        BNE appdirn1
        ADDQ #1,D6
        BRA appdirexit
appdirn1 SUBQ #1,D6

appdirexit MOVE.B D6,rparam
        RTS

D0_D3   LEA slash,A1
        MOVE.B #14,D0
        TRAP #15
        
        MOVE.B (A4),D5
        LSR.B #3,D5
        CMP #1,D5
        BNE d03n1
        LEA dregister,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B mparam,D6
        MOVE.B mparam,D5
        JSR EA
        JSR ApplyDirection
        LEA slash,A1
        MOVE.B #14,D0
        TRAP #15
d03n1   MOVE.B (A4),D5
        MOVE.B #29,D7
        LSL.L D7,D5
        ADDQ #2,D7
        LSR.L D7,D5
        CMP #1,D5
        BNE d03n2
        MOVE.B mparam,D6
        MOVE.B rparam,D5
        JSR EA
        JSR ApplyDirection
        LEA slash,A1
        MOVE.B #14,D0
        TRAP #15
d03n2   MOVE.B (A4),D5
        MOVE.B #30,D7
        LSL.L D7,D5
        ADDQ #1,D5
        LSR.L D7,D5
        CMP #1,D5
        BNE d03n3
        MOVE.B mparam,D6
        MOVE.B rparam,D5
        JSR EA
        JSR ApplyDirection
        LEA slash,A1
        MOVE.B #14,D0
        TRAP #15
d03n3   MOVE.B (A4),D5
        MOVE.B #31,D7
        LSL.L D7,D5
        LSR.L D7,D5
        CMP #1,D5
        BNE d03exit
        MOVE.B mparam,D6
        MOVE.B rparam,D5
        JSR EA
        LEA slash,A1
        MOVE.B #14,D0
        TRAP #15
        
d03exit RTS        

MULS_   *Check for invalid opcode
        MOVE.B (sourcem),D5
        CMP #$1,D5
        BEQ mulsexit

        *Print MULS
        LEA pmuls,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (sourcem),D6
        MOVE.B (source),D5
        JSR EA
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B #$0,D6
        MOVE.B (dest),D5
        JSR EA

        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        *JSR newline
mulsexit RTS

DIVU_   JSR IsolateRegisters

        *Check for valid opcode
        MOVE.B (destm),D5
        CMP #$3,D5
        BNE divuexit
        
        *Check for valid source mode
        MOVE.B (sourcem),D5
        CMP #$1,D5
        BEQ divuexit

        *Print MULS
        LEA pdivu,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (sourcem),D6
        MOVE.B (source),D5
        JSR EA
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B #$0,D6
        MOVE.B (dest),D5
        JSR EA

        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        *JSR newline
divuexit RTS


ADD_    JSR IsolateRegisters
        
        *EA + Dn to Dn
        MOVE.B (destm),D5
        CMP #$0,D5
        BNE addn1
        *Check for address register with size byte
        MOVE.B (sourcem),D6
        CMP #$1,D6
        BEQ addexit
        LEA dotb,A5
        MOVE.B #$0,destm
        BRA printadd
addn1   CMP #$1,D5
        BNE addn2
        LEA dotw,A5
        MOVE.B #$0,destm
        BRA printadd
addn2   CMP #$2,D5
        BNE addn3
        LEA dotl,A5
        MOVE.B #$0,destm
        BRA printadd
        
        *Dn + EA to EA
addn3   CMP #$4,D5
        BNE addn4
        MOVE.B (sourcem),D6
        CMP #$0,D6
        BEQ addexit
        CMP #$1,D6
        BEQ addexit
        CMP #$7,D6
        BNE addcnt1
        MOVE.B (source),D6
        CMP #$4,D6
        BEQ addexit
addcnt1 LEA dotb,A5
        MOVE.B #$0,destm
        BRA printadd
addn4   CMP #$5,D5
        BNE addn5
        MOVE.B (sourcem),D6
        CMP #$0,D6
        BEQ addexit
        CMP #$1,D6
        BEQ addexit
        CMP #$7,D6
        BNE addcnt2
        MOVE.B (source),D6
        CMP #$4,D6
        BEQ addexit
addcnt2 LEA dotw,A5
        MOVE.B #$0,destm
        BRA printadd
addn5   CMP #$6,D5
        BNE addn6
        MOVE.B (sourcem),D6
        CMP #$0,D6
        BEQ addexit
        CMP #$1,D6
        BEQ addexit
        CMP #$7,D6
        BNE addcnt3
        MOVE.B (source),D6
        CMP #$4,D6
        BEQ addexit
addcnt3 LEA dotl,A5
        MOVE.B #$0,destm
        BRA printadd
        
        *ADDA
addn6   CMP #$3,D5
        BNE addn7
        LEA dotw,A5
        MOVE.B #$1,destm
        BRA printadda
addn7   LEA dotl,A5
        MOVE.B #$1,destm
        BRA printadda
        
printadda LEA padda,A1
        BRA addrest
printadd LEA padd,A1
addrest MOVE.B #14,D0
        TRAP #15
        LEA (A5),A1
        MOVE.B #14,D0
        TRAP #15
        CLR D5
        CLR D6
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR EA
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (dest),D5
        MOVE.B (destm),D6
        JSR EA

        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        *JSR newline
addexit RTS

SUB_    JSR IsolateRegisters
        MOVE.B (destm),D5
        
        *Dn - ea to Dn
        CMP #$0,D5
        BNE subn1
        LEA dotb,A5
        MOVE.B #$0,destm
        BRA printsub
subn1   CMP #$1,D5
        BNE subn2
        LEA dotw,A5
        MOVE.B #$0,destm
        BRA printsub
subn2   CMP #$2,D5
        BNE subn3
        LEA dotl,A5
        MOVE.B #$0,destm
        BRA printsub
        
        *EA - Dn to EA
subn3   CMP #$4,D5
        BNE subn4
        MOVE.B (sourcem),D6
        CMP #$0,D6
        BEQ subexit
        CMP #$1,D6
        BEQ subexit
        CMP #$7,D6
        BNE subcnt1
        MOVE.B (source),D6
        CMP #$4,D6
        BEQ subexit
subcnt1 LEA dotb,A5
        MOVE.B #$0,destm
        BRA printsub
subn4   CMP #$5,D5
        BNE subn5
        MOVE.B (sourcem),D6
        CMP #$0,D6
        BEQ subexit
        CMP #$1,D6
        BEQ subexit
        CMP #$7,D6
        BNE subcnt2
        MOVE.B (source),D6
        CMP #$4,D6
        BEQ subexit
subcnt2 LEA dotw,A5
        MOVE.B #$0,destm
        BRA printsub
subn5   CMP #$6,D5
        BNE subexit
        MOVE.B (sourcem),D6
        CMP #$0,D6
        BEQ subexit
        CMP #$1,D6
        BEQ subexit
        CMP #$7,D6
        BNE subcnt3
        MOVE.B (source),D6
        CMP #$4,D6
        BEQ subexit
subcnt3 LEA dotl,A5
        MOVE.B #$0,destm
        
printsub LEA psub,A1
subrest MOVE.B #14,D0
        TRAP #15
        LEA (A5),A1
        MOVE.B #14,D0
        TRAP #15
        CLR D5
        CLR D6
        MOVE.B (source),D5
        MOVE.B (sourcem),D6
        JSR EA
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (dest),D5
        MOVE.B (destm),D6
        JSR EA

        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        *JSR newline
subexit RTS


LEA_    *Check for valid source mode
        CLR D5
        MOVE.B (sourcem),D5
        CMP #$2,D5
        BEQ leaprnt
        CMP #$7,D5
        BNE leaexit
        MOVE.B (source),D5
        CMP #$0,D5
        BEQ leaprnt
        CMP #$1,D5
        BNE leaexit

        *Print
leaprnt LEA plea,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B (sourcem),D6
        MOVE.B (source),D5
        JSR EA
        LEA comma,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B #$1,D6
        MOVE.B (dest),D5
        JSR EA
        
        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        *JSR newline
leaexit RTS

NOP_    MOVEA.W A3,A4
        MOVE.B -(A4),D4

        MOVE.B (A4)+,D4
        CMP #$04,D4
        BEQ n1NOP
        BRA notNOP
n1NOP   MOVE.B (A4)+,D4
        CMP #$0E,D4
        BEQ n2NOP
        BRA notNOP
n2NOP   MOVE.B (A4)+,D4
        CMP #$07,D4
        BEQ n3NOP
        BRA notNOP
n3NOP   MOVE.B (A4),D4
        CMP #$01,D4
        BNE notNOP

        LEA pnop,A1
        MOVE.B #14,D0
        TRAP #15
        
        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        *JSR newline
notNOP  RTS

RTS_    MOVEA.W A3,A4
        MOVE.B -(A4),D4

        MOVE.B (A4)+,D4
        CMP #$04,D4
        BEQ r1RTS
        BRA notRTS
r1RTS   MOVE.B (A4)+,D4
        CMP #$0E,D4
        BEQ r2RTS
        BRA notRTS
r2RTS   MOVE.B (A4)+,D4
        CMP #$07,D4
        BEQ r3RTS
        BRA notRTS
r3RTS   MOVE.B (A4),D4
        CMP #$05,D4
        BNE notRTS

        LEA prts,A1
        MOVE.B #14,D0
        TRAP #15

        MOVE.B #$CC,D5
        ADDQ #4,A3
        MOVE.B (A3),D3
        *JSR newline
notRTS  RTS

JSR_
    CLR.L D5
    *D4 ALREADY == 15->12
    ADDQ #2,A4
    MOVE.B (A4),D4  *SETS D4 == 7->4
    MOVE.B #30,D7   *isolate 5->4
    LSL.L D7,D4
    LSR.L D7,D4
    
    LSL.L #1,D4 *PREP TO ADD SPOT 3
    CLR.L D2
    MOVE.L D4,D2    *STORE 5->4 ON D2
    
    CLR.L D4
    ADDQ #1,A4  *MOVE TO RANGE 3->0
    MOVE.B (A4),D4
    LSR.L #3,D4 *ISOLATE 3
    
    ADD.L D4,D2 *D2 == 5->3
    CLR.L D4
    
    *checks 5->3
    CMP #$0,D2
    BEQ notJSR
    CMP #$1,D2
    BEQ notJSR
    CMP #$3,D2
    BEQ notJSR
    CMP #$4,D2
    BEQ notJSR
    CMP #$5,D2
    BEQ notJSR
    CMP #$6,D2
    BEQ notJSR
    
    CMP #$2,D2
    BEQ JSR_pANp
    
    CMP #$7,D2
    BEQ JSR_xxx
    
    BRA notJSR
    
JSR_pAnp
    MOVE.B (A4),D4  *GETS 3->0
    CLR.L D7
    MOVE.B #29,D7
    LSL.L D7,D4 *ISOLATE 2->0
    LSR.L D7,D4
    MOVE.B D2,D1
    
    LEA pjsr,A1
    MOVE.B #14,D0
    TRAP #15
    
    LEA leftp,A1
    MOVE.B #14,D0
    TRAP #15
    
    LEA aregister,A1
    MOVE.B #14,D0
    TRAP #15
    
    MOVE.B #3,D0
    TRAP #15
    
    LEA rightp,A1
    MOVE.B #14,D0
    TRAP #15
    
    MOVE.B #$CC,D5
    ADDQ #4, A3
    MOVE.B (A3),D3
    *BRA newline
    RTS

JSR_xxx
    MOVE.B (A4),D4  *GETS 3->0
    CLR.L D7
    MOVE.B #29,D7
    LSL.L D7,D4 *ISOLATE 2->0
    LSR.L D7,D4
    MOVE.B D4,D1
    
    CLR.L D4
    ADDQ #2,A4  * spot -1 first bit after 3->0
    MOVE.B (A4),D4
    CMP #$FF,D4
	BEQ notJSR
    CMP #$DD,D4
    BEQ notJSR
    
    SUBQ #1,A4  *gets range 3->0
	
	CLR.L D4
    
    CMP #$0,D1
    BEQ JSRea
    
    CMP #$1,D1
    BEQ JSRea
    
    BRA notJSR
    
JSRea
    LEA pjsr,A1
    MOVE.B #14,D0
    TRAP #15
    
    LEA dollarsign,A1
	MOVE.B #14,D0
	TRAP #15

    *CLR.L D4
    *ADDQ #1,A4  *first bite after last opcode $
    *MOVE.B (A4),D4
    *CMP #$FF,D4
	*BEQ notJSR
    *CMP #$DD,D4
    *BEQ notJSR
    
JSRprintea

    CLR.L D4
    MOVE.B (A4),D4
	CLR.L D2
	MOVE.B #16,D2
	CLR.L D1
	MOVE.L D4,D1	*print extra ea from D4
	MOVE.B #15,D0
	TRAP #15
	
	ADDQ #1,A4
	MOVE.B (A4),D4
	
	ADDQ #1,A3
	
	CMP #$FF,D4
	BEQ JSRdone
	CMP #$DD,D4
	BNE JSRprintea

	MOVE.B #$CC,D5  *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	*BRA newline
	RTS

JSRdone
	MOVE.B #$CC,D5  *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	*BRA newline
	RTS
    
notJSR  RTS

*============================ASL, ASR, LSL, LSR, ROL, ROR================================

rl_check    *for registers
    CLR.L D4
	ADDQ #$2,A4	*gets digits 7->4 into D4
	MOVE.B (A4),D4
	
	LSR.L #2,D4	*isolates the 7->6 spot, if it == 11 branch
	CMP #$3,D4
    BEQ rl_check2   *branch to check2 for ea
    
    CLR.L D4
	MOVE.B (A4),D4	*7->4 back into D4
	MOVE.B #31,D7
	LSL.L D7,D4	*isolates 4
	LSR.L D7,D4
	LSL.L #1,D4	*preps to get added
	CLR.L D2
	MOVE.B D4,D2	*stores 4 in D4

	ADDQ #$1,A4	*SHIFT gets digits 3->0 into D4
	MOVE.B (A4),D4
	LSR.L #3,D4	*isolate 3
	
	ADD.B D2,D4	*gets 4->3 in D4
	
	CMP #$2,D4  *if 4->3 == 2, bad
    BEQ not_rl_shift 
    
    BRA rl_shift
    
rl_check2   *ea check
    CLR.L D4
    SUBQ #1,A4  *shift 11->8
    MOVE.B (A4),D4  *get 11->8
    
    LSR.L #1,D4 *shift right to isolate 11->9
    
	CMP #$2,D4  *cant be 010
    BEQ not_rl_shift
    
    CMP #$4,D4  cant be 100
    BEQ not_rl_shift
    
    CMP #$5,D4
    BEQ not_rl_shift 
    
    CMP #$6,D4
    BEQ not_rl_shift 
    
    CMP #$7,D4
    BEQ not_rl_shift
    
    CLR.L D4
	ADDQ #$1,A4	*gets digits 7->4 into D4
	MOVE.B (A4),D4

    MOVE.B #30,D7
	LSL.L D7,D4	*isolate 5->4
	LSR.L D7,D4

	LSL.L #1,D4	*prep to add 3 onto it
	CLR.L D2
	MOVE.L D4,D2	*store 5->4 on D2 to read

	CLR.L D4
	ADDQ #1,A4	*gets 3->0
	MOVE.B (A4),D4
	
	LSR.L #3,D4	*isolate 3

	ADD.L D4,D2	*add shifted spot 5->4 to spot 3, D2 = 5->3
	
	CMP #$0,D2
	BEQ not_rl_shift
	
	CMP #$1,D2
	BEQ not_rl_shift
	
	CMP #$5,D2
	BEQ not_rl_shift
	
	CMP #$6,D2
	BEQ not_rl_shift
	
	CMP #$2,D2
	BEQ rl_shift
	
	CMP #$3,D2
	BEQ rl_shift
	
	CMP #$4,D2
	BEQ rl_shift
	
	*D4 equals 111 CONTINUE
	
    CLR.L D4
    ADDQ #1,A4  * spot -1 first bit after 3->0
    MOVE.B (A4),D4
    CMP #$FF,D4
	BEQ not_rl_shift
    CMP #$DD,D4
    BEQ not_rl_shift
    
    SUBQ #1,A4  *gets range 3->0
	
	CLR.L D4
	*ADDQ #3,A4  *get range 3->0
	MOVE.B (A4),D4  *load 3->0 onto D4
	
	CLR.L D7
	MOVE.B #29,D7
	LSL.L D7,D4     *isolate 2->0 on D4
	LSR.L D7,D4

    *SUBQ #3,A4  *reset counter back to 15->12
    
    CMP #$0,D4
    BEQ rl_shift
    
    CMP #$1,D4
    BEQ rl_shift
    
    BRA not_rl_shift

rl_shift	*is it memory or register shift?
    SUBQ #3,A4  *reset counter back to 15->12
	CLR.L D4
	ADDQ #$2,A4	*gets digits 7->4 into D4
	MOVE.B (A4),D4

	LSR.L #2,D4	*isolates the 7->6 spot
	CMP #$3,D4
	BNE rl_shift_b	*branch to part b(register)if necessary

*continue with ea shift

	CLR.L D4
	SUBQ #1,A4	*gets 11->8 number
	MOVE.B (A4),D4

	LSR.L #1,D4	*isolate 11->9

	
	CMP #$0,D4	*if 11->9 == 000, jump to print AS
	BEQ print_as_a

	CMP #$1,D4	*if 11->9 == 001, jump to print LS
	BEQ print_ls_a

	CMP #$2,D4	*if 11->9 == 010, bad op 
	BEQ not_rl_shift 

	CMP #$3,D4	*if 11->9 == 011. jump to print RO
	BEQ print_ro_a

	*could probably just BRA not_rl_shift
	
	CMP #$4,D4	*if 11->9 == 100, bad op 
	BEQ not_rl_shift

	CMP #$5,D4	*if 11->9 == 101, bad op 
	BEQ not_rl_shift 

	CMP #$6,D4	*if 11->9 == 110, bad op 
	BEQ not_rl_shift 

	CMP #$7,D4	*if 11->9 == 111, bad op 
	BEQ not_rl_shift

	BRA not_rl_shift

print_as_a  *print AS

	LEA pas,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_a1

print_ls_a  *print LS
	LEA pls,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_a1

print_ro_a  *print RO
	LEA pro,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_a1

rl_shift_a1	*ea
	CLR.L D4
	MOVE.B (A4),D4	*gets 11->8 into D4
	MOVE.B #31,D7
	LSL.L D7,D4	*isolate 8
	LSR.L D7,D4
	
	CMP #$0,D4	*if spot 8 == 0,go to print R
	BEQ print_r_a
	
	CMP #$1,D4	*if spot 8 == 1, go to print L
	BEQ print_l_a

	BRA not_rl_shift

print_r_a   *print R
	LEA pright,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_a2

print_l_a   *print L
	LEA pleft,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_a2

rl_shift_a2	*register shift
	CLR.L D4
	ADDQ #$1,A4	*gets digits 7->4 into D4
	MOVE.B (A4),D4

    MOVE.B #30,D7
	LSL.L D7,D4	*isolate 5->4
	LSR.L D7,D4

	LSL.L #1,D4	*prep to add 3 onto it
	CLR.L D2
	MOVE.L D4,D2	*store 5->4 on D2 to read

	CLR.L D4
	ADDQ #1,A4	*gets 3->0
	MOVE.B (A4),D4
	
	LSR.L #3,D4	*isolate 3

	ADD.L D4,D2	*add shifted spot 5->4 to spot 3, D2 = 5->3

	CLR.L D4

	MOVE.B (A4),D4	*load 3->0 on D4
	MOVE.B #29,D7
	LSL.L D7,D4	*isolate 2->0
	LSR.L D7,D4
	CLR.L D1	*load 2->0 into D1
	MOVE.B D4,D1

	LEA space,A1
	MOVE.B #14,D0
	TRAP #15

	CMP #$2,D2	*if == 010, print (An)
	BEQ print_pAnp_rl_a

	CMP #$3,D2	*if == 011, print (An)+
	BEQ print_pAnpplus_rl_a

	CMP #$4,D2	*if == 100, print -(An)
	BEQ print_minuspAnp_rl_a

	CMP #$7,D2	*if == 111, print ea
	BEQ print_ea_rl_a
	
	BRA not_rl_shift

print_pAnp_rl_a     *print (An)
	LEA leftp,A1
	MOVE.B #14,D0
	TRAP #15

	LEA aregister,A1
	MOVE.B #14,D0
	TRAP #15
	
	MOVE.B #3,D0
	TRAP #15
    
	LEA rightp,A1
	MOVE.B #14,D0
	TRAP #15

	MOVE.B #$CC,D5      *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	*BRA newline
	RTS

print_pAnpplus_rl_a *print (An)+
	LEA leftp,A1
	MOVE.B #14,D0
	TRAP #15

	LEA aregister,A1
	MOVE.B #14,D0
	TRAP #15
	
	MOVE.B #3,D0
	TRAP #15

	LEA rightp,A1
	MOVE.B #14,D0
	TRAP #15

	LEA plus,A1
	MOVE.B #14,D0
	TRAP #15

	MOVE.B #$CC,D5      *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	*BRA newline
	RTS

print_minuspAnp_rl_a    *print -(An)

	LEA minus,A1
	MOVE.B #14,D0
	TRAP #15

	LEA leftp,A1
	MOVE.B #14,D0
	TRAP #15

	LEA aregister,A1
	MOVE.B #14,D0
	TRAP #15
	
	MOVE.B #3,D0
	TRAP #15

	LEA rightp,A1
	MOVE.B #14,D0
	TRAP #15

	MOVE.B #$CC,D5      *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	*BRA newline
	RTS
	
print_ea_rl_a   *print ea 
	CMP #$2,D1
	BEQ not_rl_shift

	CMP #$3,D1
	BEQ not_rl_shift

	CMP #$4,D1
	BEQ not_rl_shift

	CMP #$5,D1
	BEQ not_rl_shift

	CMP #$6,D1
	BEQ not_rl_shift

	CMP #$7,D1
	BEQ not_rl_shift

	LEA dollarsign,A1
	MOVE.B #14,D0
	TRAP #15

    CLR.L D4
    ADDQ #1,A4
    MOVE.B (A4),D4
    CMP #$FF,D4
	BEQ not_rl_shift
    CMP #$DD,D4
    BEQ not_rl_shift
    
    BRA print_rl_ea
    
print_rl_ea
    CLR.L D4
    MOVE.B (A4),D4
	CLR.L D2
	MOVE.B #16,D2
	CLR.L D1
	MOVE.L D4,D1	*print extra ea from D4
	MOVE.B #15,D0
	TRAP #15
	
	ADDQ #1,A4
	MOVE.B (A4),D4
	
	ADDQ #1,A3
	
	CMP #$FF,D4
	BEQ rl_shift_done
	CMP #$DD,D4
	BNE print_rl_ea

	MOVE.B #$CC,D5  *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	*BRA newline
	RTS
	
rl_shift_done
	MOVE.B #$CC,D5  *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	*BRA newline
	RTS

rl_shift_b  *=====================PART B FOR MEMORY LOCATIONS==================
    CLR.L D4
	MOVE.B (A4),D4	*7->4 back into D4
	MOVE.B #31,D7
	LSL.L D7,D4	*isolates 4
	LSR.L D7,D4
	LSL.L #1,D4	*preps to get added
	CLR.L D2
	MOVE.B D4,D2	*stores 4 in D4

	ADDQ #$1,A4	*SHIFT gets digits 3->0 into D4
	MOVE.B (A4),D4
	LSR.L #3,D4	*isolate 3
	
	ADD.B D2,D4	*gets 4->3
	
	CMP #$0,D4
	BEQ print_as_b

	CMP #$1,D4
	BEQ print_ls_b

	CMP #$3,D4
	BEQ print_ro_b

	BRA not_rl_shift
	
print_as_b      *print AS
	LEA pas,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b1

print_ls_b      *print LS
	LEA pls,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b1

print_ro_b  *print RO
	LEA pro,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b1

rl_shift_b1
	SUBQ #2,A4      *SHIFT
	MOVE.B (A4),D4	*gets 11->8 into D4
	MOVE.B #31,D7
	LSL.L D7,D4	    *isolate 8
	LSR.L D7,D4

	CMP #$0,D4	    *if spot 8 == 0, print R
	BEQ print_r_b

	CMP #$1,D4	    *if spot 8 == 1, print L
	BEQ print_l_b

	BRA not_rl_shift

print_r_b   *print R
	LEA pright,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b2

print_l_b   *print L
	LEA pleft,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b2

rl_shift_b2
	ADDQ #1,A4      *SHIFT
	MOVE.B (A4),D4	*gets 7->4 into D4
	LSR.L #2,D4	    *isolate 7->6

	CMP #$0,D4
	BEQ print_byte_rl_b

	CMP #$1,D4
	BEQ print_word_rl_b

	CMP #$2,D4
	BEQ print_long_rl_b
	
	BRA not_rl_shift

print_byte_rl_b *print .B
	LEA dotb,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b3

print_word_rl_b *print .W
	LEA dotw,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b3

print_long_rl_b *print .L
	LEA dotl,A1
	MOVE.B #14,D0
	TRAP #15
	BRA rl_shift_b3

rl_shift_b3
    CLR.L D4
	MOVE.B (A4),D4	*gets 7->4 into D4
	LSR.L #1,D4	    *isolate 5
	MOVE.B #31,D7
	LSL.L D7,D4
	LSR.L D7,D4

	CMP #$0,D4
	BEQ print_number_rl_b

	CMP #$1,D4
	BEQ print_d_rl_b

	BRA not_rl_shift

print_number_rl_b   *prints #number
	SUBQ #1,A4	    *SHIFT gets 11->8 into D4
	MOVE.B (A4),D4

	LSR.L #1,D4	    *isolates 11->9
	CLR.L D1
	MOVE.L D4,D1
	
	LEA hashtag,A1  *print #
	MOVE.B #14,D0
	TRAP #15

	MOVE.B #3,D0    *print number in D1
	TRAP #15
	
	LEA comma,A1    *print ,
	MOVE.B #14,D0
	TRAP #15

	LEA dregister,A1	*print D
	MOVE.B #14,D0
	TRAP #15

	ADDQ #2,A4	    *SHIFT get 3->0 into D4
	MOVE.B (A4),D4

	CLR.L D1
	MOVE.B #29,D7
	LSL.L D7,D4	    *isolate 2->0
	LSR.L D7,D4
	MOVE.L D4,D1

	MOVE.B #3,D0	*print number in D1, dregister
	TRAP #15

	MOVE.B #$CC,D5  *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	*BRA newline
	RTS

print_d_rl_b
	SUBQ #1,A4	    *SHIFT gets 11->8 into D4
	MOVE.B (A4),D4

	LSR.L #1,D4	    *isolates 11->9
	CLR.L D1
	MOVE.L D4,D1
	
	LEA dregister,A1	*print D
	MOVE.B #14,D0
	Trap #15

	MOVE.B #3,D0	*print number in D1
	TRAP #15

	LEA comma,A1	*print ,
	MOVE.B #14,D0
	Trap #15

	LEA dregister,A1	*print D
	MOVE.B #14,D0
	TRAP #15

	ADDQ #2,A4  	*get 3->0 into D4
	MOVE.B (A4),D4

	CLR.L D1
	MOVE.B #29,D7
	LSL.L D7,D4	    *isolate 2->0
	LSR.L D7,D4
	MOVE.L D4,D1

	MOVE.B #3,D0	*print number in D1, dregister
	TRAP #15

	MOVE.B #$CC,D5  *DONE
	ADDQ #4,A3
	MOVE.B (A3),D3
	*BRA  newline
	RTS

not_rl_shift
	RTS
*===================================================================================


*Pre-conditions:    
*D6 needs to contain addressing mode
*D5 needs to contain EA/register number          
EA      CMP #$0,D6
        BEQ EA_000
        
        CMP #$1,D6
        BEQ EA_001
        
        CMP #$2,D6
        BEQ EA_010
        
        CMP #$3,D6
        BEQ EA_011
        
        CMP #$4,D6
        BEQ EA_100
        
        CMP #$7,D6
        BEQ EA_111
        
        MOVE.B #$EE,badop
        RTS

EA_000  LEA dregister,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B D5,D1
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
        RTS
        
EA_001  LEA aregister,A1
        MOVE.B #14,D0
        TRAP #15
        MOVE.B D5,D1
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
        RTS

EA_010  LEA leftp,A1
        MOVE.B #14,D0
        TRAP #15
        LEA aregister,A1
        TRAP #15
        MOVE.B D5,D1
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
        LEA rightp,A1
        MOVE.B #14,D0
        TRAP #15
        RTS

EA_011  LEA leftp,A1
        MOVE.B #14,D0
        TRAP #15
        LEA aregister,A1
        TRAP #15
        MOVE.B D5,D1
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
        LEA rightp,A1
        MOVE.B #14,D0
        TRAP #15
        LEA plus,A1
        TRAP #15
        RTS 

EA_100  LEA minus,A1
        MOVE.B #14,D0
        TRAP #15
        LEA leftp,A1
        TRAP #15
        LEA aregister,A1
        TRAP #15
        MOVE.B D5,D1
        MOVE.B #16,D2
        MOVE.B #15,D0
        TRAP #15
        LEA rightp,A1
        MOVE.B #14,D0
        TRAP #15
        RTS    

EA_111  CMP #$0,D5
        BEQ ABSword
        
        CMP #$1,D5
        BEQ ABSlong
        
        CMP #$4,D5
        BEQ Immediate
        
        MOVE.B #$EE,badop
        RTS 
        
ABSword LEA dollarsign,A1
        MOVE.B #14,D0
        TRAP #15
        MOVEA.W A3,A4
        ADDQ #3,A4
        CLR D4
        CLR D5
        MOVE.B #4,D4
        MOVE.B #16,D2
        MOVE.B #15,D0
abswlp  ADDQ #1,D5
        MOVE.B (A4)+,D1
        TRAP #15
        CMP D4,D5
        BNE abswlp
        ADDA.L D4,A3
        RTS
        
ABSlong LEA dollarsign,A1
        MOVE.B #14,D0
        TRAP #15
        MOVEA.W A3,A4
        ADDQ #3,A4
        CLR D4
        CLR D5
        MOVE.B #8,D4
        MOVE.B #16,D2
        MOVE.B #15,D0
absllp  ADDQ #1,D5
        MOVE.B (A4)+,D1
        TRAP #15
        CMP D4,D5
        BNE absllp
        ADDA.L D4,A3
        RTS
        
Immediate LEA hashtag,A1
        MOVE.B #14,D0
        TRAP #15
        LEA dollarsign,A1
        TRAP #15
        MOVEA.W A3,A4
        ADDQ #3,A4
        CLR D6
        CLR D5
        MOVE.B #8,D4
        MOVE.B #16,D2
        MOVE.B #15,D0
absilp  ADDQ #1,D5
        MOVE.B (A4)+,D1
        TRAP #15
        MOVE.B (A4),D6
        CMP #$DD,D6
        BNE absin1
        BRA absilpext
absin1  CMP #$FF,D6
        BNE absilp
absilpext ADDA.L D5,A3
        RTS

IsolateRegisters
        CLR D5
        CLR D6
        CLR D7
        MOVEA.W A3,A4
        SUBQ #1,A4
        
        *Isolate and get source mode
        ADDQ #2,A4
        MOVE.B (A4),D5
        ADDQ #1,A4
        MOVE.B (A4),D6
        MOVE.B #30,D7
        LSL.L D7,D5
        SUBI.B #1,D7
        LSR.L D7,D5
        LSR.B #$3,D6
        ADD.B D6,D5
        MOVE.B D5,sourcem
        
        *Isolate and get source register
        MOVEA.W A3,A4
        ADDQ #2,A4
        MOVE.B (A4),D5
        MOVE.B #29,D7
        LSL.L D7,D5
        LSR.L D7,D5
        MOVE.B D5,source
        
        *Isolate and get destination mode
        MOVEA.W A3,A4
        MOVE.B (A4),D5
        ADDQ #1,A4
        MOVE.B (A4),D6
        MOVE.B #31,D7
        LSL.L D7,D5
        SUBI.B #2,D7
        LSR.L D7,D5
        LSR.B #$2,D6
        ADD.B D5,D6
        MOVE.B D6,destm
        
        *Isolate and get destination register
        MOVEA.W A3,A4
        MOVE.B (A4),D5
        LSR.B #$1,D5
        MOVE.B D5,dest

        RTS 
      
    END    START        ; last line of source





































*~Font name~Courier New~
*~Font size~10~
*~Tab type~1~
*~Tab size~4~
