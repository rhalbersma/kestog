
#ifndef KALLISTO_H
#define KALLISTO_H

//#include <windows.h>

#ifdef __cplusplus
extern "C" {
#endif

#define KALLISTOAPI WINAPI

int         WINAPI getmove (int board[8][8], int color, double maxtime, char str[1024],
             int *playnow, int info, int unused, struct CBmove *move);
int         WINAPI enginecommand (char command[256], char reply[1024]);
int         WINAPI islegal (int b[8][8], int color, int from, int to,struct  CBmove *move);

void KALLISTOAPI EI_MakeMove(char *move);
char * KALLISTOAPI  EI_Think();
void KALLISTOAPI EI_Ponder();
char * KALLISTOAPI EI_PonderHit(char *move);
void KALLISTOAPI EI_Initialization(PF_SearchInfo si, int mem_lim);
void KALLISTOAPI EI_Stop();
void KALLISTOAPI EI_SetupBoard(char *p);
void KALLISTOAPI EI_NewGame();
void KALLISTOAPI EI_SetTimeControl(int time, int inc);
void KALLISTOAPI EI_SetTime(int time, int otime);
char *KALLISTOAPI EI_GetName();
void KALLISTOAPI EI_OnExit();
void KALLISTOAPI EI_Analyse();
void KALLISTOAPI EI_EGDB(EdAccess *eda);

#ifdef __cplusplus
   }
#endif

#endif
