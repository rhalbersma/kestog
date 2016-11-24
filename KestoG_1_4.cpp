//************************************************
//*                                               *
//*               Kallisto support                *
//*                                               *
//*************************************************
// функция обратного вызова для отображения информации о ходе вычислений
// pv - лучший вариант
// cm - ход анализируемый в данный момент
typedef void (__stdcall *PF_SearchInfo)(int score, int depth, int speed, char *pv, char *cm);
PF_SearchInfo pfSearchInfo = 0;
// This function needed to show search info

#include <fstream>
using namespace std;

//*************************************************
//*                                               *
//*               KestoG sources                  *
//*                                               *
//*************************************************

/* includes */
#define WIN32_LEAN_AND_MEAN // exclude rarely used stuff from headers
#include <stdio.h>
#include <stdlib.h>  // for malloc()
#include <string.h>  // for memset()
#include <windows.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "EdAccess.h"
#include "KALLISTO.H" // f-ion prototypes


/* board values */
#define OCCUPIED 16
#define WHITE 1
#define BLACK 2
#define MAN 4
#define KING 8
#define FREE 0
#define BLK_MAN (BLACK|MAN)
#define WHT_MAN (WHITE|MAN)
#define BLK_KNG  (BLACK|KING)
#define WHT_KNG (WHITE|KING)

#define SHIFT1 (BLK_MAN << 8)
#define SHIFT2 (WHT_MAN << 8)
#define SHIFT3 (BLK_KNG << 8)
#define SHIFT4 (WHT_KNG << 8)

#define CHANGECOLOR 3
#define MAXDEPTH 60
#define MAXMOVES 60 // ??
#define MAXHIST 16384

#define UPPER 1
#define LOWER 2
#define EXACT 3
#define MATE 10000   // de facto accepted value
#define HASHMATE 9744
#define ED_WIN 30000
#define REHASH 4
#define ETCDEPTH 4  // if depth >= ETCDEPTH do ETC
#define WINDOW 40
//#define KALLISTOAPI WINAPI
//#undef KALLISTO // uncomment if compile to CheckerBoard
#define KALLISTO // uncomment if compile to KallistoGUI

/*----------> compile options  */
// not used options
//#define PERFT
#undef PERFT
#undef MUTE
#undef VERBOSE
#undef STATISTICS

/* getmove return values */
#define DRAW 0
#define WIN 1
#define LOSS 2
#define UNKNOWN 3

#define TICKS CLK_TCK
#define SWAPINT(a,b) {int t;t=a;a=b;b=t;}

/* definitions of new types   */
typedef unsigned __int64 U64;
typedef unsigned __int32 U32;

typedef struct {
    U32 m_lock;
    int m_value:16;
    int m_depth:8;
    unsigned int m_valuetype:2;
    signed int m_best_from:8;
    signed int m_best_to:8;
    int m_age;
   } TEntry;

#define NODE_OPP(type) (-(type))
#define ISWINSCORE(score) ((score) >= HASHMATE)
#define ISLOSSSCORE(score) ((score) <= -HASHMATE)

struct coor             /* coordinate structure for board coordinates */
 {
   unsigned int x;
   unsigned int y;
 };

 struct CBmove    // all the information there is about a move
  {
  int jumps;               // how many jumps are there in this move?
  unsigned int newpiece;            // what type of piece appears on to
  unsigned int oldpiece;            // what disappears on from
  struct  coor from, to;            // coordinates of the piece in 8x8 notation!
  struct coor path[12];            // intermediate path coordinates
  struct  coor del[12];            // squares whose pieces are deleted
  unsigned int delpiece[12];            // what is on these squares
  } GCBmove;

struct move2
  {
    unsigned short m[12];
    char path[11]; // was short path[12];
    char l;// move's length
   };

/* function prototypes */

static int PVSearch(int b[46],int depth,int fracDepth,int trunc,int alpha,int beta,int color,int node_type,int iid_search,int xcapture);
static int rootsearch(int b[46], int alpha, int beta, int depth,int color,int search_type);
static int   compute( int b[46],int color, int time, char output[256]);
static int   eval(int b[46], int color, int alpha, int beta, bool in_pv);
static void                      domove(int b[46],struct move2 *move,int stm);
static void                      domove2(int b[46],struct move2 *move,int stm);
static void                      undomove(int b[46],struct move2 *move,int stm);
static void                      update_hash(struct move2 *move);
static int            Gen_Captures(int b[46],struct move2 movelist[MAXMOVES],int color );
static int            Gen_Moves(int b[46],struct move2 movelist[MAXMOVES],int color );
static int            Gen_Proms(int b[46],struct move2 movelist[10],int color );
static void black_king_capture(int b[46], int *n, struct move2 movelist[MAXMOVES], int j,int in_dir);
static void black_man_capture(int b[46], int *n, struct move2 movelist[MAXMOVES], int j,int in_dir);
static void white_king_capture(int b[46], int *n, struct move2 movelist[MAXMOVES], int j,int in_dir);
static void white_man_capture(int b[46], int *n, struct move2 movelist[MAXMOVES], int j,int in_dir);
static int            Test_Capture(int b[46], int color);
static int            Test_From_pb( int b[46],int mloc,int dir);
static int            Test_From_cb( int b[46], int mloc,int dir);
static int            Test_From_pw( int b[46],int mloc,int dir);
static int            Test_From_cw( int b[46], int mloc,int dir);

static void            setbestmove(struct move2 move);
static void            MoveToStr(move2 m, char *s);
static struct          coor numbertocoor(int n);
static void            movetonotation(struct move2 move,char str[80]);

static U64             rand64(void);
static U64             Position_to_Hashnumber( int b[46] , int color );
static void            Create_HashFunction(void);
static void            TTableInit( unsigned int size);
static void            retrievepv( int b[46], char *pv, int color);
static int hashretrieve(int depth,int *value,int *alpha,int *beta,int *best_from,int *best_to,int *try_mcp,bool in_pv);
static void                hashstore(int value,int depth,int alpha,int beta,int best_from,int best_to);
static U64                 getUniqueZobrist(void);
static BOOL                isZobristUnique(U64 key);

static __inline int           weight();
static void                   ClearHistory(void);
static void QuickSort( int SortVals[MAXMOVES],struct move2 movelist[MAXMOVES], int inf, int sup);
static void                     Perft(int b[46],int color,int depth);
static void                     good_move( int from,int to,int ply);
static __inline void    history_bad( int from,int to );
static __inline void    history_good( int from,int to );
static __inline bool ok_to_reduce( int from,int to );
static __inline bool ok_to_prune( int from,int to,int depth );
static          bool    move_is_dangerous( int b[46],struct move2 *move );
static int __inline            is_promotion(struct move2 *move);
static int                     matval(int b[46],int color);
static __inline int            get_phase(void);
static int       pick_next_move( int *marker,int SortVals[MAXMOVES],int n );
static void      Sort( int start,int num, int SortVals[],struct move2 movelist[] );
static void      init( int b[46] );
int              EdProbe(int c[46],int color);

/*----------> globals  */
CRITICAL_SECTION AnalyseSection;
int *play;
U64 ZobristNumbers[47][17];
U64 HashSTM; // random number - side to move
U64 HASH_KEY; // global variable HASH_KEY
U32 MASK;

unsigned int size = 8; // default size 8 MB
int realdepth;
int nodes; // better double nodes;?
struct move2 bestrootmove;
unsigned int History[46][46]; // array used by history heuristics
unsigned int HistHit[46][46];  // good moves which not failed-low
unsigned int HistTot[46][46]; // all played moves

// killer slots
int killersf1[MAXDEPTH+1];
int killerst1[MAXDEPTH+1];

int killersf2[MAXDEPTH+1];
int killerst2[MAXDEPTH+1];

int killersf3[MAXDEPTH+1];
int killerst3[MAXDEPTH+1];

int killer_scores1[MAXDEPTH+1];
int killer_scores2[MAXDEPTH+1];
int killer_scores3[MAXDEPTH+1];

double PerftNodes;
static const int directions[4] = { 0x4,0x5,-0x4,-0x5 };
// static const int g_dirs[4] = { 0x8,0x4,0x2,0x1 };
static const int NodeAll = -1;
static const int NodePV  =  0;
static const int NodeCut = +1;
//static int HistoryValue = 10650; // 10650 == 65%  9830 == 60% 9010==55% 11469==70%
int root_depth;
enum game_phase {OPENING = 0,MIDGAME = 1,ENDGAME = 2};
static const int MVALUE[11] = {0,0,0,0,0,100,100,0,0,300,300};
static int g_pieces[11]; // global array
double start,t,maxtime; // time variables
const int SearchNormal = 0;
const int SearchShort  = 1;
// margins for reduction based on remaining depth
int red_margin[] =
{
   0, 0, 0, 400, 500,  600,  700,  800,  900,  1000, 1100, 1200,
   1300, 1400,   1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200,
   2300, 2400,   2500, 2600, 2700, 2800, 2900, 3000, 3100, 4000
};
//int red[] = { 0,110,200,300};

static unsigned int indices[41]; // indexes into p_list for a given square
static unsigned int p_list[3][16]; // p_list contains the location of all the pieces of either color
static unsigned int captured_number[MAXDEPTH][12];

int num_wpieces;
int num_bpieces;

int searches_performed_in_game;
TEntry *ttable=NULL;

EdAccess *ED = NULL;
unsigned int EdPieces = 0;
bool Reversible[MAXDEPTH];
int EdRoot[3];
bool EdNocaptures = false;

#pragma warn -8057 // disable warning #8057 when compiling
#pragma warn -8004 // disable warning #8004 when compiling

/*    dll stuff                  */
BOOL WINAPI
DllEntryPoint (HANDLE hDLL, DWORD dwReason, LPVOID lpReserved)
{
  /* in a dll you used to have LibMain instead of WinMain in
     windows programs, or main in normal C programs win32
     replaces LibMain with DllEntryPoint. */

  switch (dwReason)
    {
    case DLL_PROCESS_ATTACH:
      /* dll loaded. put initializations here! */
      break;
    case DLL_PROCESS_DETACH:
      /* program is unloading dll. put clean up here! */
      break;
    case DLL_THREAD_ATTACH:
      break;
    case DLL_THREAD_DETACH:
      break;
    default:
      break;
    }
  return TRUE;
}


/* CheckerBoard API: enginecommand(), islegal(), getmove() */

int WINAPI
enginecommand (char str[256], char reply[1024])
{
  /* answer to commands sent by CheckerBoard.  This does not
   * answer to some of the commands, eg it has no engine
   * options. */

  char command[256], param1[256], param2[256];
  char *stopstring;
  sscanf (str, "%s %s %s", command, param1, param2);

  // check for command keywords

  if (strcmp (command, "name") == 0)
    {
      sprintf (reply, "KestoG engine v1.4");
      return 1;
    }

  if (strcmp (command, "about") == 0)
    {
     sprintf (reply,"KestoG engine v1.4\n by Kestutis Gasaitis\nEngine which plays Russian checkers\n2007\nE-mail to:kestasl8t@delfi.lt with any comments.");
     return 1;
    }

  if (strcmp (command, "help") == 0)
    {
      sprintf (reply, "missing.htm");
      return 0;
    }

  if (strcmp (command, "set") == 0)
    {
      if (strcmp (param1, "hashsize") == 0)
         {
              size = strtol( param2, &stopstring, 10 );
              if ( size < 1) return 0;
              if ( size > 128) size = 128;
              return 1;
          }
     if (strcmp (param1, "book") == 0)
   {
     return 0;
   }
    }

  if (strcmp (command, "get") == 0)
    {
      if (strcmp (param1, "hashsize") == 0)
   {
     return 0;
   }
      if (strcmp (param1, "book") == 0)
   {
     return 0;
   }
      if (strcmp (param1, "protocolversion") == 0)
   {
     sprintf (reply, "2");
     return 1;
   }
      if (strcmp (param1, "gametype") == 0)
   {
     sprintf (reply, "25"); // 25 stands for Russian checkers
     return 1;
   }
    }
    strcpy (reply, "?");
  return 0;
}


int WINAPI
islegal (int b[8][8], int color, int from, int to,struct CBmove *move)
{
  /* islegal tells CheckerBoard if a move the user wants to
   * make is legal or not. to check this, we generate a
   * movelist and compare the moves in the movelist to the
   * move the user wants to make with from & to */

    int n,i,found=0,Lfrom,Lto;
    struct move2 movelist[MAXMOVES];
    int board[46];
    int capture;
    char Lstr[80];

   /* initialize board */
   for(i=0;i<46;i++)
     board[i]=OCCUPIED;
   for(i=5;i<=40;i++)
     board[i]=FREE;
       board[5]=b[0][0];board[6]=b[2][0];board[7]=b[4][0];board[8]=b[6][0];
       board[10]=b[1][1];board[11]=b[3][1];board[12]=b[5][1];board[13]=b[7][1];
       board[14]=b[0][2];board[15]=b[2][2];board[16]=b[4][2];board[17]=b[6][2];
       board[19]=b[1][3];board[20]=b[3][3];board[21]=b[5][3];board[22]=b[7][3];
       board[23]=b[0][4];board[24]=b[2][4];board[25]=b[4][4];board[26]=b[6][4];
       board[28]=b[1][5];board[29]=b[3][5];board[30]=b[5][5];board[31]=b[7][5];
       board[32]=b[0][6];board[33]=b[2][6];board[34]=b[4][6];board[35]=b[6][6];
       board[37]=b[1][7];board[38]=b[3][7];board[39]=b[5][7];board[40]=b[7][7];
   for(i=5;i<=40;i++)
     if(board[i] == 0) board[i]=FREE;
   for(i=9;i<=36;i+=9)
     board[i]=OCCUPIED;
     init(board);

     /* board initialized */

     n = Gen_Captures( board,movelist,color );
     capture=n;

     if (!n)
     n = Gen_Moves( board,movelist,color );
     if (!n) return 0;

     /* now we have a movelist - check if from and to are the same */
    for(i=0;i<n;i++)
      {
                     movetonotation(movelist[i],Lstr);
                     if ( capture )
                        sscanf(Lstr,"%i%*c%i",&Lfrom,&Lto);
                    else
                        sscanf(Lstr,"%i%*c%i",&Lfrom,&Lto);

                    if(from==Lfrom && to==Lto)
         {
         found=1;
         break;
         }
      }
      if(found){
      /* sets GCBmove to movelist[i] */
       setbestmove(movelist[i]);
      *move=GCBmove;
        }

    return found;
   }


int WINAPI
getmove (int board[8][8],int color,double maxtime,char str[1024],int *playnow,int info,int unused,struct CBmove *move)
{
  /* getmove is what checkerboard calls. you get the parameters:

     - b[8][8]
     is the current position. the values in the array are
     determined by the #defined values of BLACK, WHITE, KING,
     MAN. a black king for instance is represented by BLACK|KING.

     - color
     is the side to make a move. BLACK or WHITE.

     - maxtime
     is the time your program should use to make a move. this
     is what you specify as level in checkerboard. so if you
     exceed this time it's
not too bad - just
don't exceed it
     too much...

     - str
     is a pointer to the output string of the checkerboard status bar.
     you can use sprintf(str,"information"); to print any information you
     want into the status bar.

     - *playnow
     is a pointer to the playnow variable of checkerboard. if
     the user would like your engine to play immediately, this
     value is nonzero, else zero. you should respond to a
     nonzero value of *playnow by interrupting your search
     IMMEDIATELY.

     - CBmove
     tells checkerboard what your move is, see above.
   */

   int i;
   int value;
   int b[46];
   int time = (int)(maxtime * 1000.0);
   double t, elapsed;
   // Create_HashFunction();

   /* initialize board */
   for(i=0;i<46;i++)
     b[i]=OCCUPIED;
   for(i=5;i<=40;i++)
     b[i]=FREE;
          /*    (white)
                37  38  39  40
              32  33  34  35
                28  29  30  31
              23  24  25  26
                19  20  21  22
              14  15  16  17
                10  11  12  13
               5   6   7   8
         (black)   */
     b[5]=board[0][0];b[6]=board[2][0];b[7]=board[4][0];b[8]=board[6][0];
     b[10]=board[1][1];b[11]=board[3][1];b[12]=board[5][1];b[13]=board[7][1];
     b[14]=board[0][2];b[15]=board[2][2];b[16]=board[4][2];b[17]=board[6][2];
     b[19]=board[1][3];b[20]=board[3][3];b[21]=board[5][3];b[22]=board[7][3];
     b[23]=board[0][4];b[24]=board[2][4];b[25]=board[4][4];b[26]=board[6][4];
     b[28]=board[1][5];b[29]=board[3][5];b[30]=board[5][5];b[31]=board[7][5];
     b[32]=board[0][6];b[33]=board[2][6];b[34]=board[4][6];b[35]=board[6][6];
     b[37]=board[1][7];b[38]=board[3][7];b[39]=board[5][7];b[40]=board[7][7];
   for(i=5;i<=40;i++)
       if ( b[i] == 0 ) b[i]=FREE;
   for(i=9;i<=36;i+=9)
       b[i]=OCCUPIED;
       play=playnow;

#ifdef PERFT
       t=clock();
       init(b);
       PerftNodes = 0;
       realdepth = 0;
       Perft(b,color,10);
       elapsed = (clock()-t)/(double)CLK_TCK;
       sprintf(str,"[done][time %.2fs][PerftNodes %.0f]" ,elapsed,PerftNodes);
#else
      init(b);
      if ( info&1 ){
      TTableInit(size);
      searches_performed_in_game = 0;
      Create_HashFunction();
                    }
      value = compute(b, color,time, str);
#endif

   for(i=5;i<=40;i++)
       if( b[i] == FREE) b[i]=0;
     /* return the board */
    board[0][0]=b[5];board[2][0]=b[6];board[4][0]=b[7];board[6][0]=b[8];
    board[1][1]=b[10];board[3][1]=b[11];board[5][1]=b[12];board[7][1]=b[13];
    board[0][2]=b[14];board[2][2]=b[15];board[4][2]=b[16];board[6][2]=b[17];
    board[1][3]=b[19];board[3][3]=b[20];board[5][3]=b[21];board[7][3]=b[22];
    board[0][4]=b[23];board[2][4]=b[24];board[4][4]=b[25];board[6][4]=b[26];
    board[1][5]=b[28];board[3][5]=b[29];board[5][5]=b[30];board[7][5]=b[31];
    board[0][6]=b[32];board[2][6]=b[33];board[4][6]=b[34];board[6][6]=b[35];
    board[1][7]=b[37];board[3][7]=b[38];board[5][7]=b[39];board[7][7]=b[40];

    /* set the move */
   *move=GCBmove;
    if(value>=HASHMATE) return WIN;
    if(value<=-HASHMATE) return LOSS;
  return UNKNOWN;
}


static void black_king_capture( int b[46],  int *n, struct move2 movelist[MAXMOVES],int j,int in_dir){
   //
   int m; // auxiliary variable
   int temp; // temporary variable
   int capsq; // captured square address
   int found_cap = 0;
   int found_pd;
   struct move2 move,orgmove;
   int i; // loop counter
   int dir;
   int next_dir;

   orgmove = movelist[*n];
               for ( i =0; i < 4 ; i++ ){      // scan all 4 directions
               if ( !(in_dir >> (( 3 - i ) & 1)) ) continue;
               dir = directions[i];             // dir = 4,5,-4,-5.
               temp = j; // from square
               do temp += dir; while ( !b[temp] );
               if ( ( b[temp] & WHITE ) != 0 ){     /*   ----===========----   */
                 temp = temp + dir;
                 if ( b[temp] ) continue;
                 capsq = temp - dir; // fix captured square address in capsq
                 found_pd = 0;
                 do{
                      if ( Test_From_pb(b,temp,dir) ){
                            // add to movelist
                            move = orgmove;
                            move.l++;
                            move.path[move.l - 2] = temp;
                            m = SHIFT3|temp;
                            move.m[1] = m;
                            m = (b[capsq]<<8)|capsq;
                            move.m[move.l - 1] = m;
                            found_pd++;
                            found_cap = 1;
                            movelist[*n] = move;
                            b[capsq]++;
     // further jumps
     switch (i){
     case 0:next_dir = ( found_pd == 1 ) ? 13:5;break; // in binary form 1101:0101
     case 1:next_dir = ( found_pd == 1 ) ? 14:10;break; // in binary form 1110:1010
     case 2:next_dir = ( found_pd == 1 ) ? 7:5;break; // in binary form 0111:0101
     case 3:next_dir = ( found_pd == 1 ) ? 11:10;break; // in binary form 1011:1010
               }
                            black_king_capture(b, n, movelist, temp,next_dir);
                            b[capsq]--;
                                                                  } // if
                            temp = temp + dir;
                            } while ( !b[temp] );

                            if ( !found_pd ){
                            if ( (b[temp] & WHITE) != 0 && !b[temp+dir] ){
                            temp = capsq + dir;
                                  // add to movelist
                                  move = orgmove;
                                  move.l++;
                                  move.path[move.l - 2] = temp;
                                  m = SHIFT3|temp;
                                  move.m[1] = m;
                                  m = (b[capsq]<<8)|capsq;
                                  move.m[move.l - 1] = m;
                                  found_cap = 1;
                                  movelist[*n] = move;
                                  b[capsq]++;
                                  // further 1 jump
                                  // next_dir = g_dirs[i]; // 8,4,2,1
                                  switch (i){
                                  case 0:next_dir = 8;break;
                                  case 1:next_dir = 4;break;
                                  case 2:next_dir = 2;break;
                                  case 3:next_dir = 1;break;
                                               }
                                  black_king_capture(b, n, movelist, temp,next_dir);
                                  b[capsq]--;
                                } // if
                                  else{
                                  temp = capsq + dir;
                                  do{
                                  // add to movelist
                                  move = orgmove;
                                  move.l++;
                                  move.path[move.l - 2] = temp;
                                  m = SHIFT3|temp;
                                  move.m[1] = m;
                                  m = (b[capsq]<<8)|capsq;
                                  move.m[move.l - 1] = m;
                                  found_cap = 1;
                                  movelist[*n] = move;
                                  (*n)++;
                                  temp = temp + dir;
                                       }while ( !b[temp] );
                                          }
                                       }
                     } //   /*   ----===========----   */
                } // for

     if ( !found_cap ) (*n)++;
 }


static void black_man_capture( int b[46],  int *n,struct move2 movelist[MAXMOVES] ,int j,int in_dir){
 //
   int m;   // auxiliary variable
   int found_cap = 0;
   struct move2 move,orgmove;
   int i;
   int dir;

   orgmove = movelist[*n];
              for ( i = 0; i < 4; i++ ){ // scan all 4 directions
                      dir = directions[i];   // dir = 4,5,-4,-5.
                      if ( dir == in_dir ) continue;
                      int sq1 = j + dir;
                      if (  ( b[sq1] & WHITE ) != 0 ){
                      int sq2 = j + (dir<<1);
                          if ( !b[sq2] ){
                                     // add to movelist
                                     move = orgmove;
                                     move.l++;
                                     move.path[move.l - 2] = sq2;
                                     m = (b[sq1]<<8)|sq1;
                                     move.m[move.l - 1] = m;
                                     if ( sq2 >= 37 ){  // promotion
                                     m = SHIFT3|sq2;
                                     move.m[1] = m;
                                     movelist[*n] = move;
                                     found_cap = 1;
                                     if ( sq2 == 37 || sq2 == 40 )
                                     (*n)++;
                                     else{
                                     dir = (dir == 4)?1:2;
                                     b[sq1]++;
                                     //
                                     black_king_capture(b, n, movelist,sq2,dir);
                                     b[sq1]--;
                                            }
                                                            }
                                     else{     // non-promotion
                                     m = SHIFT1|sq2;
                                     move.m[1] = m;
                                     found_cap = 1;
                                     movelist[*n] = move;
                                     b[sq1]++;
                                     //
                                     black_man_capture(b, n, movelist,sq2,-dir);
                                     b[sq1]--;
                                           }
                                     }
                                  }
                            } // for
        if ( !found_cap ) (*n)++;
}


static void white_king_capture( int b[46],  int *n, struct move2 movelist[MAXMOVES] ,int j,int in_dir){
   //
   int m;  // auxiliary variable
   int temp; // temporary variable
   int capsq;  // captured square address
   int found_cap = 0;
   int found_pd;
   struct move2 move,orgmove;
   int i;
   int next_dir;
   int dir;

   orgmove = movelist[*n];
              for ( i = 0; i < 4; i++ ){    // scan all 4 directions
               if ( !(in_dir >> (( 3 - i) & 1))  ) continue;
              dir = directions[i];             // dir = 4,5,-4,-5.
              temp = j; // from square
              do temp += dir;while ( !b[temp] );
              if ( ( b[temp] & BLACK ) != 0 ){ // /*-----------====================*/
              temp = temp + dir;
              if ( b[temp] ) continue;
              capsq = temp - dir; // fix captured square address in capsq
              found_pd = 0;
               do{
                    if ( Test_From_pw(b,temp,dir) ){
                     // add to movelist
                     move = orgmove;
                     move.l++;
                     move.path[move.l - 2] = temp;
                     m = SHIFT4|temp;
                     move.m[1] = m;
                     m = (b[capsq]<<8)|capsq;
                     move.m[move.l - 1] = m;
                     found_pd++;
                     found_cap = 1;
                     movelist[*n] = move;
                     b[capsq]--;
                     // further jumps
         switch (i){
         case 0:next_dir = ( found_pd == 1 ) ? 13:5;break; // in binary form 1101:0101
         case 1:next_dir = ( found_pd == 1 ) ? 14:10;break; // in binary form 1110:1010
         case 2:next_dir = ( found_pd == 1 ) ? 7:5;break; // in binary form 0111:0101
         case 3:next_dir = ( found_pd == 1 ) ? 11:10;break; // in binary form 1011:1010
                   }
                     white_king_capture(b, n, movelist, temp,next_dir);
                     b[capsq]++;
                                 } // if
                     temp = temp + dir;
                     } while ( !b[temp] );

                     if ( !found_pd ){
                     if ( (b[temp] & BLACK) != 0 && !b[temp+dir] ){
                            temp = capsq + dir;
                            // add to movelist
                            move = orgmove;
                            move.l++;
                            move.path[move.l - 2] = temp;
                            m = SHIFT4|temp;
                            move.m[1] = m;
                            m = (b[capsq]<<8)|capsq;
                            move.m[move.l - 1] = m;
                            found_cap = 1;
                            movelist[*n] = move;
                            b[capsq]--;
                             // further 1  jump
                             // next_dir = g_dirs[i]; // 8,4,2,1
                               switch (i){
                                  case 0:next_dir = 8;break;
                                  case 1:next_dir = 4;break;
                                  case 2:next_dir = 2;break;
                                  case 3:next_dir = 1;break;
                                             }
                            white_king_capture(b, n, movelist, temp,next_dir);
                             b[capsq]++;
                               } // if
                            else{
                              temp = capsq + dir;
                              do{
                              // add to movelist
                              move = orgmove;
                              move.l++;
                              move.path[move.l - 2] = temp;
                              m = SHIFT4|temp;
                              move.m[1] = m;
                              m = (b[capsq]<<8)|capsq;
                              move.m[move.l - 1] = m;
                              found_cap = 1;
                              movelist[*n] = move;
                              (*n)++;
                              temp = temp + dir;
                                  } while ( !b[temp] );
                                  }
                              }
                  } // /*-----------====================*/
        } // for

     if ( !found_cap ) (*n)++;
 }


static void white_man_capture( int b[46],  int *n, struct move2 movelist[MAXMOVES] , int j,int in_dir){
   //
   int m;  // auxiliary variable
   int found_cap = 0;
   struct move2 move,orgmove;
   int i;
   int dir;

   orgmove = movelist[*n];
                    for ( i = 0; i < 4; i++ ){ // scan all 4 directions
                             dir = directions[i]; // dir = 4,5,-4,-5.
                             if ( dir == in_dir ) continue;
                             int sq1 = j + dir;
                             if (  ( b[sq1] & BLACK ) != 0 ){
                             int sq2 = j + (dir<<1);
                             if ( !b[sq2] ){
                                    // add to movelist
                                    move = orgmove;
                                    move.l++;
                                    move.path[move.l - 2] = sq2;
                                    m = (b[sq1]<<8)|sq1;
                                    move.m[move.l - 1] = m;
                             if ( sq2 <= 8 ){  // promotion
                                   m = SHIFT4|sq2;
                                   move.m[1] = m;
                                   found_cap = 1;
                                   movelist[*n] = move;
                                   if ( sq2 == 5 || sq2 == 8 )
                                   (*n)++;
                                   else{
                                   dir = ( dir == -4 ) ? 4:8;
                                   b[sq1]--;
                                   //
                                   white_king_capture(b, n, movelist,sq2,dir);
                                   b[sq1]++;
                                         }
                                                }
                                  else{  // non-promotion
                                    m = SHIFT2|sq2;
                                    move.m[1] = m;
                                    found_cap = 1;
                                    movelist[*n] = move;
                                    b[sq1]--;
                                    //
                                    white_man_capture(b, n, movelist,sq2,-dir);
                                    b[sq1]++;
                                        }
                                   }
                               }
                          } // for
       if ( !found_cap ) (*n)++;
   }


static int Test_From_pb(int b[46], int temp, int dir){
// test if there is capture in perpendicular direction to dir from square
// for black color
int d;
int square;

       if ((dir&1) == 0)
         d = 5;
      else
         d = 4;

      square = temp;

          do{
            square = square + d;
              }while ( !b[square]  );
           if ( ( b[square] & WHITE ) != 0 )
           if ( !b[square + d] )
           return (1);

     square = temp;
     d = -d; // another perp. direction
          do{
            square = square + d;
              }while ( !b[square] );
           if ( ( b[square] & WHITE ) != 0 )
           if ( !b[square + d] )
           return (1);
    return (0);
 }


static int Test_From_cb( int b[46],int temp,int dir){
// test if there is capture in current direction dir
// for black color
          do{
            temp = temp + dir;
              }while ( !b[temp] );
           if ( ( b[temp] & WHITE ) != 0 )
           if ( !b[temp + dir] )
           return (1);
           return (0);
 }


static int Test_From_pw(int b[46],int temp,int dir){
// test if there is capture in perpendicular direction to dir from square
// for white color
int d;
int square;

       if ((dir&1) == 0)
         d = 5;
       else
         d = 4;

       square = temp;

          do{
            square = square + d;
              }while ( !b[square] );
           if ( ( b[square] & BLACK ) != 0 )
           if ( !b[square + d] )
           return (1);

       square = temp;
       d = -d; // another perp. direction

          do{
            square = square + d;
              }while ( !b[square] );
           if ( ( b[square] & BLACK ) != 0 )
           if ( !b[square + d] )
           return (1);
    return (0);
 }


static int Test_From_cw(int b[46],int temp,int dir){
// test if there is capture in current direction dir from square
// for white color
          do{
            temp = temp + dir;
              }while ( !b[temp] );
           if ( ( b[temp] & BLACK ) != 0 )
           if ( !b[temp + dir] )
           return (1);
           return (0);
 }


static int  Test_Capture(int b[46], int color){
      //
      int i;
      int square;

        if ( color == BLACK ){
        for( i = 1;i <= num_bpieces;i++){
            square = p_list[2][i];
             if ( !square )
                  continue;
            if ( (b[square] & MAN) != 0 ){
                       if( (b[square+4] & WHITE) !=0)
                       if( !b[square+8] )
                           return(1);
                       if( (b[square+5] & WHITE) !=0)
                       if( !b[square+10] )
                           return(1);
                       if( (b[square-4] & WHITE) !=0)
                       if( !b[square-8] )
                           return(1);
                       if( (b[square-5] & WHITE) !=0)
                       if( !b[square-10] )
                           return(1);
                                                          }
            else{ // KING
                  if ( Test_From_cb(b,square,4) ) return (1);
                  if ( Test_From_cb(b,square,5) ) return (1);
                  if ( Test_From_cb(b,square,-4) ) return (1);
                  if ( Test_From_cb(b,square,-5) ) return (1);
                  }
                         }
            return (0);
                                   } // if ( color == BLACK )

      if ( color == WHITE ){
        for( i = num_wpieces;i >= 1;i-- ){
           square = p_list[1][i];
           if ( !square )
                continue;
           if ( (b[square] & MAN) != 0 ){
                       if( (b[square+4] & BLACK) !=0)
                       if( !b[square+8] )
                           return(1);
                       if( (b[square+5] & BLACK) !=0)
                       if( !b[square+10] )
                           return(1);
                       if( (b[square-4] & BLACK) !=0)
                       if( !b[square-8] )
                           return(1);
                       if( (b[square-5] & BLACK) !=0)
                       if( !b[square-10] )
                           return(1);
                                                        }
         else{  // KING
                if ( Test_From_cw(b,square,4) ) return (1);
                if ( Test_From_cw(b,square,5) ) return (1);
                if ( Test_From_cw(b,square,-4) ) return (1);
                if ( Test_From_cw(b,square,-5) ) return (1);
                }
            }
            return (0);
                            } // if ( color == WHITE )
      return (0);
     }


static int Gen_Captures( int b[46], struct move2 movelist[MAXMOVES] , int color ){
/*------------------------> purpose: generate all possible captures.
   ------------------------> returns: number of captures.
*/
  int n = 0;              // move number or number of moves
  int m;                   // auxiliary variable
  int square,i;          // loop counters
  int temp;               //  temporary variable
  int capsq;             // captured square address
  int found_pd;       // found jump in perpendicular direction  to current direction
  int dir;                 // current direction
  int next_dir;         // next direction
  int j;                     // from square

          if ( color == BLACK ){
          for ( square=1; square <= num_bpieces;square++ ){
                  j = p_list[2][square];
                  if ( !j )
                         continue;

//          if ( (b[j] & BLACK ) != 0 ){
          if ( (b[j] & MAN ) != 0 ){
                b[j] = FREE;
                for ( i = 0; i < 4; i++ ){          // scan all 4 directions
                       dir = directions[i];           // dir = 4,5,-4,-5.
                        int sq1 = j + dir;
                        if ( ( b[sq1] & WHITE ) != 0 ){
                        int sq2 = j + (dir<<1);
                        if ( !b[sq2] ){
                        // add to movelist
                        movelist[n].l = 3;
                        movelist[n].path[1] = sq2;
                        m = SHIFT1|j;
                        movelist[n].m[0] = m;
                        m = (b[sq1]<<8)|sq1;
                        movelist[n].m[2] = m;
                        if ( sq2 >= 37 ){ // promotion
                        m = SHIFT3|sq2;
                        movelist[n].m[1] = m;
                        if ( sq2 == 37 || sq2 == 40 )
                        n++;
                        else{
                        next_dir = (dir == 4)?1:2; // in binary form 0001:0010
                        b[sq1]++;
                        // assert dir != -4 dir != -5 because can't promote capturing backwards
                        black_king_capture(b, &n, movelist,sq2,next_dir);
                        b[sq1]--;
                               }
                                               }
                        else{ // non-promotion
                        m = SHIFT1|sq2;
                        movelist[n].m[1] = m;
                        b[sq1]++;
                        //next_dir = -dir;
                        black_man_capture(b, &n, movelist,sq2,-dir);
                        b[sq1]--;
                             }
                           } // if
                        } // if
                  } // for
               b[j] = BLK_MAN;
                              } // if MAN

                      else{     // b[j] is a KING
                      b[j] = FREE;
                      for ( i = 0; i < 4; i++ ){          // scan all 4 directions
                           dir = directions[i];             // dir = 4,5,-4,-5.
                           temp = j; // from square
                           do temp += dir;while ( !b[temp] );
                           if ( ( b[temp] & WHITE ) != 0 ){       /*   ----===========----   */
                           temp = temp + dir;
                           if ( b[temp] ) continue;
                           capsq = temp - dir;
                           found_pd = 0;
                           do{
                           if ( Test_From_pb(b,temp,dir) ){
                             found_pd++;
                             // add to movelist
                             movelist[n].l = 3;
                             movelist[n].path[1] = temp;
                             m = SHIFT3|temp;
                             movelist[n].m[1] = m;
                             m = SHIFT3|j;
                             movelist[n].m[0] = m;
                             m = (b[capsq]<<8)|capsq;
                             movelist[n].m[2] = m;
                             b[capsq]++;
                             // further jumps
                             switch (i){
                              case 0:next_dir = ( found_pd == 1 ) ? 13:5;break; // in binary form 1101:0101
                              case 1:next_dir = ( found_pd == 1 ) ? 14:10;break; // in binary form 1110:1010
                              case 2:next_dir = ( found_pd == 1 ) ? 7:5;break; // in binary form 0111:0101
                              case 3:next_dir = ( found_pd == 1 ) ? 11:10;break; // in binary form 1011:1010
                                           }
                             black_king_capture(b, &n, movelist,temp,next_dir);
                             b[capsq]--;
                                                                    } // if
                             temp = temp + dir;
                                } while ( !b[temp] );

                         if ( !found_pd ){
                         if ( (b[temp] & WHITE) != 0 && !b[temp + dir] ){
                             temp = capsq + dir;
                             // add to movelist
                             movelist[n].l = 3;
                             movelist[n].path[1] = temp;
                             m = SHIFT3|temp;
                             movelist[n].m[1] = m;
                             m = SHIFT3|j;
                             movelist[n].m[0] = m;
                             m = (b[capsq]<<8)|capsq;
                             movelist[n].m[2] = m;
                             b[capsq]++;
                             // further 1 jump
                             // next_dir = g_dirs[i]; // 8,4,2,1
                                switch (i){
                                  case 0:next_dir = 8;break;
                                  case 1:next_dir = 4;break;
                                  case 2:next_dir = 2;break;
                                  case 3:next_dir = 1;break;
                                               }
                             black_king_capture(b, &n, movelist,temp,next_dir);
                             b[capsq]--;
                                                                                                      } // if
                             else{
                                 temp = capsq + dir;
                                 do{
                                 // add to movelist
                                 movelist[n].l = 3;
                                 movelist[n].path[1] = temp;
                                 m = SHIFT3|temp;
                                 movelist[n].m[1] = m;
                                 m = SHIFT3|j;
                                 movelist[n].m[0] = m;
                                 m = (b[capsq]<<8)|capsq;
                                 movelist[n].m[2] = m;
                                 n++;
                                 temp = temp + dir;
                                    }while ( !b[temp] );
                                    }
                                        }
                      }     /*   ----===========----   */
               } // for
           b[j] = BLK_KNG;
        } // else
//   } // if BLACK
                 } // for
    } // if ( color == BLACK )

       if ( color == WHITE ){
                  for ( square=1; square <= num_wpieces;square++ ){
                  j = p_list[1][square];
                  if ( !j )
                       continue;

//         if ( ( b[j] & WHITE ) != 0 ){
         if ( ( b[j] & MAN ) != 0 ){
               b[j] = FREE;
               for ( i = 0; i < 4; i++ ){          // scan all 4 directions
                    dir = directions[i];             // dir = 4,5,-4,-5.
                     int sq1 = j + dir;
                     if ( ( b[sq1] & BLACK ) != 0 ){
                     int sq2 = j + (dir<<1);
                     if ( !b[sq2] ){
                         // add to movelist
                         movelist[n].l = 3;
                         movelist[n].path[1] = sq2;
                         m = SHIFT2|j;
                         movelist[n].m[0] = m;
                         m = (b[sq1]<<8)|(sq1);
                         movelist[n].m[2] = m;
                         if ( sq2 <= 8 ){ // promotion
                         m = SHIFT4|sq2;
                         movelist[n].m[1] = m;
                         if ( sq2 == 5 || sq2 == 8 )
                         n++;
                         else{
                         next_dir = (dir == -4) ? 4:8; // in binary form 0100:1000
                         b[sq1]--;
                         // assert dir != 4 && dir != 5
                         white_king_capture(b,&n,movelist,sq2,next_dir);
                         b[sq1]++;
                                }
                                              }
                         else{ // non-promotion
                         m = SHIFT2|sq2;
                         movelist[n].m[1] = m;
                         b[sq1]--;
                         //next_dir = -dir;
                         white_man_capture(b, &n, movelist,sq2,-dir);
                         b[sq1]++;
                               }
                              } // if
                            } // if
                        } // for
                        b[j] =WHT_MAN;
               }    // if MAN

              else{ // b[j] is a KING
                         b[j] = FREE;
                         for ( i = 0; i < 4; i++ ){          // scan all 4 directions
                              dir = directions[i];             // dir = 4,5,-4,-5.
                               temp = j; // from square
                              do  temp += dir;while ( !b[temp] );
                               if ( ( b[temp] & BLACK ) != 0 ){   /*   ----===========----   */
                               temp = temp + dir;
                               if ( b[temp] ) continue;
                               capsq = temp - dir;
                               found_pd = 0;
                               do{
                                if ( Test_From_pw(b,temp,dir) ){
                                   found_pd++;
                                   // add to movelist
                                   movelist[n].l = 3;
                                   movelist[n].path[1] = temp;
                                   m = SHIFT4|temp;
                                   movelist[n].m[1] = m;
                                   m = SHIFT4|j;
                                   movelist[n].m[0] = m;
                                   m = ( b[capsq]<<8 )|capsq;
                                   movelist[n].m[2] = m;
                                   b[capsq]--;
                                   // further jumps
                                   switch (i){
           case 0:next_dir = ( found_pd == 1 ) ? 13:5;break; // in binary form 1101:0101
           case 1:next_dir = ( found_pd == 1 ) ? 14:10;break; // in binary form 1110:1010
           case 2:next_dir = ( found_pd == 1 ) ? 7:5;break; // in binary form 0111:0101
           case 3:next_dir = ( found_pd == 1 ) ? 11:10;break; // in binary form 1011:1010
                                              }
                                   white_king_capture(b, &n, movelist,temp,next_dir);
                                   b[capsq]++;
                                                   } // if
                                   temp = temp + dir;
                                   } while ( !b[temp] );

                                   if ( !found_pd ){
                                   if ( (b[temp] & BLACK) != 0 && !b[temp + dir] ){
                                      temp = capsq + dir;
                                      // add to movelist
                                      movelist[n].l = 3;
                                      movelist[n].path[1] = temp;
                                      m = SHIFT4|temp;
                                      movelist[n].m[1] = m;
                                      m = SHIFT4|j;
                                      movelist[n].m[0] = m;
                                      m = (b[capsq]<<8)|capsq;
                                      movelist[n].m[2] = m;
                                      b[capsq]--;
                                      // further 1 jump
                                      // next_dir = g_dirs[i]; // 8,4,2,1
                                         switch (i){
                                         case 0:next_dir = 8;break;
                                         case 1:next_dir = 4;break;
                                         case 2:next_dir = 2;break;
                                         case 3:next_dir = 1;break;
                                                       }
                                      white_king_capture(b,&n, movelist,temp,next_dir);
                                      b[capsq]++;
                                                                                                                  } // if
                                 else{
                                 temp = capsq + dir;
                                 do{
                                 // add to movelist
                                 movelist[n].l = 3;
                                 movelist[n].path[1] = temp;
                                 m = SHIFT4|temp;
                                 movelist[n].m[1] = m;
                                 m = SHIFT4|j;
                                 movelist[n].m[0] = m;
                                 m = (b[capsq]<<8)|capsq;
                                 movelist[n].m[2] = m;
                                 n++;
                                 temp = temp + dir;
                                      }while ( !b[temp] );
                                        }
                                              }
               }  /*   ----====================----   */
           } // for
     b[j] = WHT_KNG;
           } // else
  //  } // if WHITE
                                       } // for
   } //  if ( color == WHITE )

   return (n);  // returns number of captures n
   }


static int Gen_Moves( int b[46], struct move2 movelist[MAXMOVES] , int color ){
   //
   int m; // auxiliary variable
   int n = 0;
   int i;
   int temp;
   int square;

         if ( color == BLACK ){
         for ( i = 1;i <= num_bpieces;i++){
          square = p_list[2][i];
          if ( !square )
                  continue;
                  if ( ( b[square] & MAN ) != 0 ){
                            if ( !b[square + 4] ){
                               movelist[n].l = 2;
                               if ( square >=32 ) m=BLK_KNG; else m=BLK_MAN; m = (m<<8)|(square+4);
                               movelist[n].m[1] = m;
                               m = SHIFT1|square;
                               movelist[n++].m[0] = m;
                                                           }
                            if ( !b[square + 5] ){
                               movelist[n].l = 2;
                               if ( square >=32 ) m=BLK_KNG; else m=BLK_MAN; m = (m<<8)|(square+5);
                               movelist[n].m[1] = m;
                               m = SHIFT1|square;
                               movelist[n++].m[0] = m;
                                                          }
                                                                  } // MAN
                 else{   // KING
                        temp = square + 4;
                        while  ( !b[temp] ){
                        movelist[n].l = 2;
                        m = SHIFT3|temp;
                        movelist[n].m[1] = m;
                        m = SHIFT3|square;
                        movelist[n++].m[0] = m;
                        temp = temp + 4;
                                           } // while
                        temp = square + 5;
                        while  ( !b[temp] ){
                        movelist[n].l = 2;
                        m = SHIFT3|temp;
                        movelist[n].m[1] = m;
                        m = SHIFT3|square;
                        movelist[n++].m[0] = m;
                        temp = temp + 5;
                                            } // while
                        temp = square - 4;
                        while  ( !b[temp] ){
                        movelist[n].l = 2;
                        m = SHIFT3|temp;
                        movelist[n].m[1] = m;
                        m = SHIFT3|square;
                        movelist[n++].m[0] = m;
                        temp = temp - 4;
                                           } // while
                        temp = square - 5;
                        while  ( !b[temp] ){
                        movelist[n].l = 2;
                        m = SHIFT3|temp;
                        movelist[n].m[1] = m;
                        m = SHIFT3|square;
                        movelist[n++].m[0] = m;
                        temp = temp - 5;
                                            } // while
                        }  // else
                               } // for
           } // color == BLACK

          else{ // if ( color == WHITE )
              for ( i = 1;i <=num_wpieces;i++){
                    square = p_list[1][i];
                    if ( !square )
                            continue;
                                   if ( ( b[square] & MAN ) != 0 ){
                                    if ( !b[square-4] ){
                                            movelist[n].l = 2;
                    if ( square <=13 ) m=WHT_KNG; else m=WHT_MAN; m = (m<<8)|(square-4);
                                            movelist[n].m[1] = m;
                                            m = SHIFT2|square;
                                            movelist[n++].m[0] = m;
                                                                }
                                            if ( !b[square-5] ){
                                            movelist[n].l = 2;
                   if ( square <=13 ) m=WHT_KNG; else m=WHT_MAN; m = (m<<8)|(square-5);
                                            movelist[n].m[1] = m;
                                            m = SHIFT2|square;
                                            movelist[n++].m[0] = m;
                                                                }
                                                      } // MAN
                           else{ // KING
                             temp = square + 4;
                             while  ( !b[temp] ){
                             movelist[n].l = 2;
                             m = SHIFT4|temp;
                             movelist[n].m[1] = m;
                             m = SHIFT4|square;
                             movelist[n++].m[0] = m;
                             temp = temp + 4;
                                                        } // while
                             temp = square + 5;
                             while  ( !b[temp] ){
                             movelist[n].l = 2;
                             m = SHIFT4|temp;
                             movelist[n].m[1] = m;
                             m = SHIFT4|square;
                             movelist[n++].m[0] = m;
                             temp = temp + 5;
                                                        } // while
                             temp = square - 4;
                             while  ( !b[temp] ){
                             movelist[n].l = 2;
                             m = SHIFT4|temp;
                             movelist[n].m[1] = m;
                             m = SHIFT4|square;
                             movelist[n++].m[0] = m;
                             temp = temp - 4;
                                                        } // while
                             temp = square - 5;
                             while  ( !b[temp] ){
                             movelist[n].l = 2;
                             m = SHIFT4|temp;
                             movelist[n].m[1] = m;
                             m = SHIFT4|square;
                             movelist[n++].m[0] = m;
                             temp = temp - 5;
                                                        } // while
                               } // else
                  } // for
      }
      return (n); // returns number of moves n
}


static int Gen_Proms(int b[46],struct move2 movelist[10],int color){
   // generates only promotions
   int m;
   int n = 0;
   int i;
   int dir;
   int j;

         if ( color == BLACK ){
         for ( i = 32;i <= 35;i++){
           if ( (b[i] & BLACK) != 0 ){
              if ( (b[i] & MAN) != 0 ){
                      for ( j = 0; j <= 1; j++ ){ // scan 2 directions
                           dir = directions[j];  // dir = 4,5.
                            if ( !b[i+dir] ){
                               movelist[n].l = 2;
                               m = SHIFT3|(i+dir);
                               movelist[n].m[1] = m;
                               m = SHIFT1|i;
                               movelist[n++].m[0] = m;
                                             }
                                                 } // for
                             } // MAN
                          } // BLACK
                     } // for
           } // if ( color == BLACK )
          // or else
          if ( color == WHITE ){
          for ( i = 10;i<=13;i++){
                   if ( (b[i] & WHITE) != 0 ){
                   if ( (b[i] & MAN) != 0 ){
                   for ( j = 2; j <= 3; j++ ){ // scan 2 directions
                   dir = directions[j];        // dir = -4,-5.
                   if ( !b[i+dir] ){
                   movelist[n].l = 2;
                   m = SHIFT4|(i+dir);
                   movelist[n].m[1] = m;
                   m = SHIFT2|i;
                   movelist[n++].m[0] = m;
                                             }
                                  } // for
                         } // MAN
                   } // WHITE
            } // for
   } // if ( color == WHITE )
  return (n); // returns number of promotions n
}


static void domove(int b[46],struct move2 *move,int stm)
/*----> purpose: execute move on board and update HASH_KEY */
{
   unsigned int from,target;
   unsigned int contents;
   int i;

   HASH_KEY ^= HashSTM;

   from = ((move->m[0]) & 63);
   b[from] = FREE;
   contents = ((move->m[0]) >> 8);
   Reversible[realdepth] = (contents & KING) && move->l == 2;

   HASH_KEY ^= ZobristNumbers[from][contents];
   g_pieces[contents]--;

//   plist[pos[square]] = 0;

   target = (( move->m[1]) & 63);
   contents = ((move->m[1]) >> 8);
   HASH_KEY ^= ZobristNumbers[target][contents];
   b[target] = contents;
   g_pieces[contents]++;

   indices[target] = indices[from];
   p_list[stm][indices[target]] = target;

//   plist[pos[square]] = square;

   for(i=2;i<move->l;i++){
      target = ((move->m[i]) & 63);
      b[target] = FREE;
      contents = ((move->m[i]) >> 8);
      HASH_KEY ^= ZobristNumbers[target][contents];
      g_pieces[contents]--;
      captured_number[realdepth][i] = indices[target];
      p_list[stm^3][indices[target]] = 0;
//plist[pos[square]] = 0;
                                      }
      realdepth++;
}


static void domove2(int b[46],struct move2 *move,int stm )
/*----> purpose: execute move on board without HASH_KEY updating */
{
   unsigned int from,target;
   unsigned int contents;
   int i;

   from = ((move->m[0]) & 63);
   b[from] = FREE;
   contents = ((move->m[0]) >> 8);
   Reversible[realdepth] = (contents & KING) && move->l == 2;
   g_pieces[contents]--;

   //   plist[pos[square]] = 0;

   target = ((move->m[1]) & 63);
   contents = ((move->m[1]) >> 8);
   b[target] = contents;
   g_pieces[contents]++;

//   plist[pos[square]] = square;

    indices[target] = indices[from];
    p_list[stm][indices[target]] = target;

    for(i=2;i<move->l;i++){
      target = ((move->m[i]) & 63);
      b[target] = FREE;
      contents = ((move->m[i]) >> 8);
      g_pieces[contents]--;
      captured_number[realdepth][i] = indices[target];
      p_list[stm^3][indices[target]] = 0;
//    plist[pos[square]] = 0;
                                        }
      realdepth++;
}


static void undomove(int b[46],struct move2 *move,int stm )
/*----> purpose: undo what domove did */
{
   unsigned int from,target;
   unsigned int contents;
   int i;

   realdepth--;
   target = ((move->m[1]) & 63);
   b[target] = FREE;
   contents = ((move->m[1]) >> 8);
   g_pieces[contents]--;

//   plist[pos[square]] = 0;

   from = ((move->m[0]) & 63);
   contents = ((move->m[0]) >> 8);
   b[from] = contents;
   g_pieces[contents]++;

//   plist[pos[square]] = square;

    indices[from] = indices[target];
    p_list[stm][indices[from]] = from;

   for(i=2;i<move->l;i++){
      target = ((move->m[i]) & 63);
      contents = ((move->m[i]) >> 8);
      b[target] = contents;
      g_pieces[contents]++;
      indices[target] = captured_number[realdepth][i];
      p_list[stm^3][indices[target]] = target;
//    plist[pos[square]] = square;
                                       }
 //     realdepth--;
}


static int eval(int b[46], int color, int alpha, int beta, bool in_pv){
   /*----> purpose: static evaluation of the board */
   int nbm = g_pieces[6]; // number of black men
   int nwm = g_pieces[5]; // number of white men
   int nbk = g_pieces[10]; // number of black kings
   int nwk = g_pieces[9]; // number of white kings

   if ( nbm == 0 && nbk == 0 )  return ( color == BLACK ? ( realdepth-MATE):MATE-realdepth);
   if ( nwm == 0 && nwk == 0 ) return( color == WHITE  ? ( realdepth-MATE):MATE-realdepth);

   int i,j;
   int eval;
   int v1,v2;
   unsigned int phase = get_phase(); // get game phase

   if ( phase == ENDGAME ){
           v1 = 100*nbm+300*nbk;
           v2 = 100*nwm+300*nwk;
           eval = v1-v2;
           int White = nwm + nwk; // total number of white pieces
           int Black = nbm + nbk; // total number of black pieces

           if ( nbk > 0 && ( White < (2+nbk)) && (eval < 0)) return (0);
           if ( nwk > 0 && ( Black < (2+nwk)) && (eval > 0)) return (0);

                                            int WGL = 0; // white king on a1-h8
                                            int BGL = 0; // black king on a1-h8
                                            for ( i=5;i<=40;i+=5 ){
                                            if ( b[i] ){
                                            if ( b[i] & MAN ){
                                            WGL=0;
                                            BGL=0;
                                            break;
                                                                        }
                                            if ( b[i] == WHT_KNG ) WGL = 1;
                                            else
                                            if ( b[i] == BLK_KNG ) BGL = 1;
                                                          }
                                                                                }

                    // surely winning advantage:
                    if ( White == 1 && nwm == 1 && Black >= 4 ) eval = eval + (eval>>1);
                    if ( Black == 1 && nbm == 1 && White >= 4 ) eval = eval + (eval>>1);

                    // scaling down
                    if ( nbk > 0 && eval < 0 ) eval = eval >> 1;
                    if ( nwk > 0 && eval > 0 ) eval = eval >> 1;

                    if ( nbk == 1 && BGL && !WGL && White <= 3 )
                    if ( Black <= 2 || eval < 500 )
                    return (0);

                    if ( nwk == 1 && WGL && !BGL && Black <= 3 )
                    if ( White <= 2 || eval > -500 )
                    return (0);

      static int PST_king[41] = {0,0,0,0,0,  // 0..4
                                               2,1,0,2,0, // 5..9
                                               2,1,2,2, // 10..13
                                               1,2,3,2,0, // 14..18
                                               1,3,3,0, // 19..22
                                               0,3,3,1,0, // 23..27
                                               2,3,2,1, // 28..31
                                               2,2,1,2,0, // 32..36
                                               2,0,1,2 // 37..40
                                              };

      static int be_man[41] = {0,0,0,0,0,                   // 0 .. 4
                                            0,0,0,0,0,                  // 5 .. 9
                                            2,2,2,2,                     // 10 .. 13
                                            4,4,4,4,0,               // 14 .. 18
                                            6,6,6,6,                    // 19 .. 22
                                            8,8,8,8,0,              // 23 .. 27
                                            10,10,10,4,                   //  28 .. 31
                                             4,12,12,12,0,       // 32 .. 36
                                             0,0,0,0                      // 37 .. 40
                                             };

    static int we_man[41] = {0,0,0,0,0,                        // 0 .. 4
                                           0,0,0,0,0,                  // 5 .. 9
                                          12,12,12,4,                     // 10 .. 13
                                           4,10,10,10,0,               // 14 .. 18
                                           8,8,8,8,                     // 19 .. 22
                                           6,6,6,6,0,              // 23 .. 27
                                           4,4,4,4,                         //  28 .. 31
                                           2,2,2,2,0,          // 32 .. 36
                                           0,0,0,0                   // 37 .. 40
                                            };

    static int LatticeArray[] = {0,0,0,0,0,  // 0 .. 4
                                             0,0,0,0,0,  // 5 .. 9
                                             2,2,2,0,     // 10 .. 13
                                              0,-2,-2,-2,0,   // 14 .. 18
                                              2,2,2,0,           // 19 .. 22
                                              0,-2,-2,-2,0,    // 23 .. 27
                                              2,2,2,0,            //  28 .. 31
                                              0,-2,-2,-2,0,  // 32 .. 36
                                              0,0,0,0                   // 37 .. 40
                                            };

                       int w_lattice = 0;
                       int b_lattice = 0;
                       for ( i = 5; i <= 40; i++ ){
                       if ( b[i] ){
                          switch ( b[i] ){
                          case BLK_MAN:
                          eval = eval + be_man[i];b_lattice+=LatticeArray[i];break;
                          case WHT_MAN:
                          eval = eval - we_man[i];w_lattice+=LatticeArray[i];break;
                          case BLK_KNG:
                          eval = eval + PST_king[i];break;
                          case WHT_KNG:
                          eval = eval - PST_king[i];break;
                          default:break;
                                              }
                                       }
                                                            }
    w_lattice = abs(w_lattice);
    if (w_lattice) eval += w_lattice - 2;
    b_lattice = abs(b_lattice);
    if (b_lattice) eval -= b_lattice - 2;
                /* the move */
    if ( nbk == 0 && nwk == 0 && nbm == nwm ){
    int allstones = nbm + nwm;
    int move;
    const int themoveval = 2;
    int stonesinsystem = 0;
    if ( color == BLACK )
         {
    for ( i=5; i <= 8;i++)
             {
    for ( j=0; j < 4; j++)
                    {
           if ( b[i+9*j] != FREE ) stonesinsystem++;
                     }
              }
     if ( stonesinsystem % 2 ) // the number of stones in blacks system is odd -> he has the move
     move = themoveval*(24-allstones)/6;
     else
     move = -themoveval*(24-allstones)/6;
     eval += move;
        }
     else // color is WHITE
     {
     for ( i=10; i <= 13;i++)
              {
     for ( j=0; j < 4; j++)
                  {
            if ( b[i+9*j] != FREE ) stonesinsystem++;
                   }
              }
     if ( stonesinsystem % 2 ) // the number of stones in whites system is odd -> he has the move
     move = -themoveval*(24-allstones)/6;
     else
     move = themoveval*(24-allstones)/6;
     eval += move;
     }
                    }

                // negamax formulation requires this:
                if(color == BLACK){
                eval++; // small bonus for turn
                return (eval);
                                                }
               else{
               eval--; // small bonus for turn
               return (-eval);
                     }
                                     }  // ENDGAME

            v1 = 100*nbm+250*nbk;
            v2 = 100*nwm+250*nwk;
            eval = v1-v2;
            eval += (200*(v1-v2))/(v1+v2);      /*favor exchanges if in material plus*/

            // king's balance
            if ( nbk != nwk){
            if ( nwk == 0 && nbm >= nwm-2 )
                 eval += 200;
            else
            if ( nbk == 0 && nwm >= nbm-2 )
                 eval -= 200;
                                    }
            // scaling down
            if ( nbk > 0 && eval < 0 ) eval = ((3*eval) >> 2);
            if ( nwk > 0 && eval > 0 ) eval = ((3*eval) >> 2);

             // Lazy evaluation
             // Early exit from evaluation  if eval already is extremely low or extremely high
             if ( !in_pv ){
             int teval = ( color == WHITE ) ? -eval : eval;
             if (teval - 64 >= beta) return teval;
             if (teval + 64 <= alpha) return teval;
                              }
   // back rank guard:
   static int br[32] = {0,-1,1,0,3,3,3,3,2,2,2,2,4,4,8,7,1,0,1,0,3,3,3,3,2,2,2,2,4,4,8,7}; // back rank values
   int code;
   int backrank;
   code = 0;
   if(b[5] & MAN) code++;
   if(b[6] & MAN) code+=2;
   if(b[7] & MAN) code+=4; // Golden checker
   if(b[8] & MAN) code+=8;
   if(b[13] == BLK_MAN) code+=16;
   backrank = br[code];
   code = 0;
   if(b[37] & MAN) code+=8;
   if(b[38] & MAN) code+=4; // Golden checker
   if(b[39] & MAN) code+=2;
   if(b[40] & MAN) code++;
   if(b[32] == WHT_MAN) code+=16;
   backrank -= br[code];
   int brv = (phase == OPENING ? 3:1);  // multiplier for back rank -- back rank value
   eval += brv*backrank;

   if ( nbm == nwm ){
            /* balance                */
            /* how equally the pieces are distributed on the left and right sides of the board */

    int balance;
    static int value[17]={0,0,0,0,0,1,256,0,0,16,4096,0,0,0,0,0,0};
    int nbml,nbmr; // number black men left - right
    int nwml,nwmr; // number white men left - right
    // left flank
    code = 0;
    // count file-a men ( on 5,14,23,32 )
    code+=value[b[5]];
    code+=value[b[14]];
    code+=value[b[23]];
    code+=value[b[32]];
    // count file-b men ( on 10,19,28,37 )
    code+=value[b[10]];
    code+=value[b[19]];
    code+=value[b[28]];
    code+=value[b[37]];
    // count file-c men ( on 6,15,24,33 )
    code+=value[b[6]];
    code+=value[b[15]];
    code+=value[b[24]];
    code+=value[b[33]];
    nwml = code & 15;
    nbml = (code>>8) & 15;
    // empty left flank ?
    if ( nwml == 0 ) eval += 10;
    if ( nbml == 0 ) eval -= 10;
    // right flank
    code = 0;
    // count file-f men ( on 12,21,30,39 )
    code+=value[b[12]];
    code+=value[b[21]];
    code+=value[b[30]];
    code+=value[b[39]];
    // count file-g men ( on 8,17,26,35 )
    code+=value[b[8]];
    code+=value[b[17]];
    code+=value[b[26]];
    code+=value[b[35]];
    // count file-h men ( on 13,22,31,40 )
    code+=value[b[13]];
    code+=value[b[22]];
    code+=value[b[31]];
    code+=value[b[40]];
    nwmr = code & 15;
    nbmr = (code>>8) & 15;
    // empty right flank ?
    if ( nwmr == 0 ) eval += 10;
    if ( nbmr == 0 ) eval -= 10;

    balance = abs(nbml-nbmr);
    if ( balance >= 2 )
    eval -= balance<<1;
    balance = abs(nwml-nwmr);
    if ( balance >= 2 )
    eval += balance<<1;

    if ( nbml + nbmr == nbm ) eval -= 4;
    if ( nwml + nwmr == nwm ) eval += 4;

    }

   // developed single corner
   const int devsinglecornerval = 8; // developed single corner value
   if ( !b[5] && !b[10] ){
   eval += devsinglecornerval;
   if ( !b[6] )
   eval -= 5;
                                    }

   if ( !b[40] && !b[35] ){
   eval -= devsinglecornerval;
   if ( !b[39] )
   eval += 5;
                                      }

     /* center control */
     // for black color
                          // f4
     if ( b[21] ){
     if ( (b[21] & BLACK) != 0 ){
     eval++;
     if ( (b[16] & BLACK) != 0 )
     if ( (b[11] & BLACK) != 0 )
     eval++;
     else
     if ( !b[11] )
     eval--;

     if ( (b[17] & BLACK) != 0 )
     if ( (b[13] & BLACK) != 0 )
     eval++;
     else
     if ( !b[13] )
     eval--;
                                               }
                     }

                           // d4
     if ( b[20] ){
     if ( (b[20] & BLACK) != 0 ){
     eval++;
     if ( (b[15] & BLACK) != 0 )
     if ( (b[10] & BLACK) != 0 )
     eval++;
     else
     if ( !b[10] )
     eval--;

     if ( (b[16] & BLACK) != 0 )
     if ( (b[12] & BLACK) != 0 )
     eval++;
     else
     if ( !b[12] )
     eval--;
                                              }
                    }

                            // e3
     if ( b[16] ){
     if ( (b[16] & BLACK) != 0 ){
     eval++;
     if ( (b[11] & BLACK) != 0 )
     if ( (b[6] & BLACK) != 0 )
     eval++;
     else
     if ( !b[6] )
     eval--;

     if ( (b[12] & BLACK) != 0 )
     if ( (b[8] & BLACK) != 0 )
     eval++;
     else
     if ( !b[8] )
     eval--;
                                                }
                   }

                             // c5
     if ( b[24] ){
     if ( (b[24] & BLACK) != 0 && b[28] == FREE && b[29] == FREE && b[39] == FREE ){
     eval+=2;
     if ( (b[19] & BLACK) != 0 )
     if ( (b[14] & BLACK) != 0 )
     eval++;
     else
     if ( !b[14] )
     eval--;

     if ( (b[20] & BLACK) != 0 )
     if ( (b[16] & BLACK) != 0 )
     eval++;
     else
     if ( !b[16] )
     eval--;
                                        }
                           }

                            // d6
     if ( b[29] ){
     if ( (b[29] & BLACK) != 0 && b[38] == FREE ){
     eval+=2;
     if ( (b[24] & BLACK) != 0 )
     if ( (b[19] & BLACK) != 0 )
     eval++;
     else
     if ( !b[19] )
     eval--;

     if ( (b[25] & BLACK) != 0 )
     if ( (b[21] & BLACK) != 0 )
     eval++;
     else
     if ( !b[21] )
     eval--;
                                                    }
                                    }

                        // f6
     if ( b[30] ){
     if ( (b[30] & BLACK) != 0 && b[39] == FREE ){
     eval+=2;
     if ( (b[25] & BLACK) != 0 )
     if ( (b[20] & BLACK) != 0 )
     eval++;
     else
     if ( !b[20] )
     eval--;

     if ( (b[26] & BLACK) != 0 )
     if ( (b[22] & BLACK) != 0 )
     eval++;
     else
     if ( !b[22] )
     eval--;
                                                                 }
                                    }
     // for white color
                    // c5
     if ( b[24] ){
     if ( (b[24] & WHITE) != 0 ){
     eval--;
     if ( (b[28] & WHITE) != 0 )
     if ( (b[32] & WHITE) != 0 )
     eval--;
     else
     if ( !b[32] )
     eval++;

     if ( (b[29] & WHITE) != 0 )
     if ( (b[34] & WHITE) != 0 )
     eval--;
     else
     if ( !b[34] )
     eval++;
                                                  }
                                    }

                      // d6
     if ( b[29] ){
     if ( (b[29] & WHITE) != 0 ){
     eval--;
     if ( (b[33] & WHITE) != 0 )
     if ( (b[37] & WHITE) != 0 )
     eval--;
     else
     if ( !b[37] )
     eval++;

     if ( (b[34] & WHITE) != 0 )
     if ( (b[39] & WHITE) != 0 )
     eval--;
     else
     if ( !b[39] )
     eval++;
                                              }
                                    }
                     // e5
     if ( b[25] ){
     if ( (b[25] & WHITE) != 0 ){
     eval--;
     if ( (b[29] & WHITE) != 0 )
     if ( (b[33] & WHITE) != 0 )
     eval--;
     else
     if ( !b[33]  )
     eval++;

     if ( (b[30] & WHITE) != 0 )
     if ( (b[35] & WHITE) != 0 )
     eval--;
     else
     if ( !b[35] )
     eval++;
                                                  }
                                    }
                   // c3
     if ( b[15] ){
     if ( (b[15] & WHITE) != 0 && b[6] == FREE ){
     eval-=2;
     if ( (b[19] & WHITE) != 0 )
     if ( (b[23] & WHITE) != 0 )
     eval--;
     else
     if ( !b[23] )
     eval++;

     if ( (b[20] & WHITE) != 0 )
     if ( (b[25] & WHITE) != 0 )
     eval--;
     else
     if ( !b[25] )
     eval++;
                                                      }
                                   }
               // e3
     if ( b[16] ){
     if ( (b[16] & WHITE) != 0 && b[7] == FREE ){
     eval-=2;
     if ( (b[20] & WHITE) != 0 )
     if ( (b[24] & WHITE) != 0 )
     eval--;
     else
     if ( !b[24] )
     eval++;

     if ( (b[21] & WHITE) != 0 )
     if ( (b[26] & WHITE) != 0 )
     eval--;
     else
     if ( !b[26] )
     eval++;
                                                  }
                                    }
              // f4
     if ( b[21] ){
     if ( (b[21] & WHITE) != 0 && b[16] == FREE && b[17] == FREE && b[6] == FREE ){
     eval-=2;
     if ( (b[25] & WHITE) != 0 )
     if ( (b[29] & WHITE) != 0 )
     eval--;
     else
     if ( !b[29] )
     eval++;

     if ( (b[26] & WHITE) != 0 )
     if ( (b[31] & WHITE) != 0 )
     eval--;
     else
     if ( !b[31] )
     eval++;
                                                         }
                                   }

                   /*  edge squares         */

                   // h2
     if ( b[13] ){
     if ( (b[13] & BLACK) != 0 )
     eval--;
                      }
                   // h4
     if ( b[22] ){
     if ( (b[22] & BLACK) != 0 )
     eval--;
                     }
                    // a3
     if ( b[14] ){
     if ( (b[14] & BLACK) != 0 )
     eval--;
                     }
                     // a5
     if ( b[23] ){
     if ( (b[23] & BLACK) != 0 )
     eval++;
                     }
                      // h6
     if ( b[31] ){
     if ( (b[31] & BLACK) != 0 && b[39] == FREE )
     eval++;
                        }
                     // a7
     if ( b[32] ){
     if ( (b[32] & BLACK) != 0 )
     eval++;
                     }
     // for white color
                  // a7
     if ( b[32] ){
     if ( (b[32] & WHITE) != 0 )
     eval++;
                     }
                   // a5
     if ( b[23] ){
     if ( (b[23] & WHITE) != 0 )
     eval++;
                      }
                   // h6
     if ( b[31] ){
     if ( (b[31] & WHITE) != 0 )
     eval++;
                       }
                   // a3
     if ( b[14] ){
     if ( (b[14] & WHITE) != 0 && b[6] == FREE )
     eval--;
                     }
                   // h4
     if ( b[22] ){
     if ( (b[22] & WHITE) != 0 )
     eval--;
                       }
                    // h2
     if ( b[13] ){
     if ( (b[13] & WHITE) != 0 )
     eval--;
                     }
                     // e5
     if ( b[25] ){
     if ( (b[25] & BLACK) != 0 )
     if ( phase != OPENING ){
     eval++;
     if (( b[20] & BLACK ) != 0)
     if (( b[15] & BLACK ) != 0)
     eval++;
     else
     if ( !b[15] )
     eval--;

     if (( b[21] & BLACK ) != 0)
     if (( b[17] & BLACK ) != 0)
     eval++;
     else
     if ( !b[17] )
     eval--;
                                         }
     else
     eval-=4;
                                    }

     if ( b[20] ){
     if ( (b[20] & WHITE) != 0 )
     if ( phase != OPENING ){
     eval--;
     if (( b[24] & WHITE ) != 0)
     if (( b[28] & WHITE ) != 0)
     eval--;
     else
     if ( !b[28] )
     eval++;

     if (( b[25] & WHITE ) != 0)
     if (( b[30] & WHITE ) != 0)
     eval--;
     else
     if ( !b[30] )
     eval++;
                                             }
     else
     eval+=4;
                                    }

      // reward checkers that will king on the next move:
      int p_bonus = ( phase == OPENING ? 8:16 ); // promote in one bonus
      for (i = 32; i <= 35; i++){
      if ( (b[i] & BLACK) != 0 )
      if ( (b[i] & MAN) != 0 )
      if ( !b[i + 5] || !b[i + 4] ) eval += color == BLACK ? p_bonus<<1 : p_bonus;
                                            }

      for (i = 10; i <= 13; i++){
      if ( (b[i] & WHITE ) != 0 )
      if ( (b[i] & MAN) != 0 )
      if ( !b[i - 5] || !b[i - 4] ) eval -= color == WHITE ? p_bonus<<1 : p_bonus;
                                            }

    static int LatticeArray[] = {0,0,0,0,0,  // 0 .. 4
                                             0,0,0,0,0,  // 5 .. 9
                                             4,4,4,0,     // 10 .. 13
                                              0,-4,-4,-4,0,   // 14 .. 18
                                              4,4,4,0,           // 19 .. 22
                                              0,-4,-4,-4,0,    // 23 .. 27
                                              4,4,4,0,            //  28 .. 31
                                              0,-4,-4,-4,0  // 32 .. 36
                                            };

    int w_lattice = 0;
    int b_lattice = 0;
    for (i = 10; i <= 35; i++){
    if ( b[i] ){
    if ( b[i] == BLK_MAN )
        b_lattice+=LatticeArray[i];
    else
    if ( b[i] == WHT_MAN )
        w_lattice+=LatticeArray[i];
                 }
                                           }

    w_lattice = abs(w_lattice);
    if (w_lattice) eval += w_lattice - 2;
    b_lattice = abs(b_lattice);
    if (b_lattice) eval -= b_lattice - 2;

    int turn = 2;
    if ( phase == OPENING ){
    if ( b[15] == BLK_MAN )
    eval--;
    if ( b[16] == BLK_MAN )
    eval--;
    if ( b[8] == BLK_MAN )
    eval+=5;
    if ( b[30] == WHT_MAN )
    eval++;
    if ( b[29] == WHT_MAN )
    eval++;
    if ( b[37] == WHT_MAN )
    eval-=5;
                                             }

   // int turn = ( phase == OPENING ? 3:2 );  // color to move gets +turn
   // negamax formulation requires this:
   if(color == BLACK){
    eval+=turn;
    return (eval);
                                   }
   else{
    eval-=turn;
    return (-eval);
         }

      }


static struct  coor numbertocoor(int n)
        {
    /* turns square number n into a coordinate for checkerboard */
  struct coor c;

   switch(n)
      {
      case 5:
         c.x=0;c.y=0;
         break;
      case 6:
         c.x=2;c.y=0;
         break;
      case 7:
         c.x=4;c.y=0;
         break;
      case 8:
         c.x=6;c.y=0;
         break;
      case 10:
         c.x=1;c.y=1;
         break;
      case 11:
         c.x=3;c.y=1;
         break;
      case 12:
         c.x=5;c.y=1;
         break;
      case 13:
         c.x=7;c.y=1;
         break;
      case 14:
         c.x=0;c.y=2;
         break;
      case 15:
         c.x=2;c.y=2;
         break;
      case 16:
         c.x=4;c.y=2;
         break;
      case 17:
         c.x=6;c.y=2;
         break;
      case 19:
         c.x=1;c.y=3;
         break;
      case 20:
         c.x=3;c.y=3;
         break;
      case 21:
         c.x=5;c.y=3;
         break;
      case 22:
         c.x=7;c.y=3;
         break;
      case 23:
         c.x=0;c.y=4;
         break;
      case 24:
         c.x=2;c.y=4;
         break;
      case 25:
         c.x=4;c.y=4;
         break;
      case 26:
         c.x=6;c.y=4;
         break;
      case 28:
         c.x=1;c.y=5;
         break;
      case 29:
         c.x=3;c.y=5;
         break;
      case 30:
         c.x=5;c.y=5;
         break;
      case 31:
         c.x=7;c.y=5;
         break;
      case 32:
         c.x=0;c.y=6;
         break;
      case 33:
         c.x=2;c.y=6;
         break;
      case 34:
         c.x=4;c.y=6;
         break;
      case 35:
         c.x=6;c.y=6;
         break;
      case 37:
         c.x=1;c.y=7;
         break;
      case 38:
         c.x=3;c.y=7;
         break;
      case 39:
         c.x=5;c.y=7;
         break;
      case 40:
         c.x=7;c.y=7;
         break;
      }
     return c;
   }


static void setbestmove( struct move2 move)
{
   int i;
   unsigned int from, to;
   int jumps;
   struct coor c1;

   jumps = move.l -2;

   from = (move.m[0]) & 255;
   to = (move.m[1]) & 255;

   GCBmove.from = numbertocoor(from);
   GCBmove.to = numbertocoor(to);
   GCBmove.jumps = jumps;
   GCBmove.newpiece =  ( (move.m[1]) >> 8 );
   GCBmove.oldpiece =  ( (move.m[0]) >> 8 );


   for ( i = 2; i < move.l; i++ ){
              GCBmove.del[i-2] = numbertocoor( (unsigned)((move.m[i]) & 255 ));
              GCBmove.delpiece[i-2] = ( (unsigned)((move.m[i]) >> 8));
                                             }

  if ( jumps > 1 )
     {
         for ( i = 2; i < move.l; i++ )
                      {
                            c1 = numbertocoor( move.path[i - 1] );
                            GCBmove.path[i - 1] = c1;
                      }
     }
 else
  {
   GCBmove.path[1] = numbertocoor(to);
  }

}


static void movetonotation(struct move2 move,char str[80])
{
   unsigned int j,from,to;
   char c;

   from=(move.m[0]) & 63;
   to=(move.m[1]) & 63;
   from=from-(from/9);
   to=to-(to/9);
   from-=5;
   to-=5;
   j=from%4;from-=j;j=3-j;from+=j;
   j=to%4;to-=j;j=3-j;to+=j;
   from++;
   to++;
   c='-';
   if(move.l>2) c='x'; // capture or normal ?
   sprintf(str,"%2li%c%2li",from,c,to);
  }


             // hash functions
static U64 rand64(void){
   // generates 64 bit random number
   U64 temp = rand();
    temp ^= ((U64)rand() << 15);
    temp ^= ((U64)rand() << 30);
    temp ^= ((U64)rand() << 45);
    temp ^= ((U64)rand() << 60);
    return temp;
}


static void Create_HashFunction(void){
   int p,q;
   srand((unsigned int) time(NULL));
   for ( p=0; p<=46 ; p++ )
      for ( q=0; q <=16 ; q++ )
      ZobristNumbers[p][q] = getUniqueZobrist();
      HashSTM = ZobristNumbers[0][0]; // constant random number - side to move
}


static U64 getUniqueZobrist(void){
  U64 key = rand64();
  while(isZobristUnique(key) != TRUE)
    key = rand64();
    return key;
}


static BOOL isZobristUnique(U64 key)
{
  int p,q;
  for(p=0;p<=46;p++ )
    for(q=0;q<=16;q++ )
        if(ZobristNumbers[p][q] == key)
          return FALSE;
  return TRUE;
}


static U64 Position_to_Hashnumber( int b[46] , int color )
{
  U64 CheckSum = 0;
  int cpos;

  for ( cpos=5; cpos<=40; cpos++ ) {
    if  ( ( b[cpos] != OCCUPIED ) && ( b[cpos] != FREE ) )
     CheckSum ^= ZobristNumbers[cpos][b[cpos]];
                                                         }

     if ( color == BLACK )
          CheckSum ^= HashSTM;

  return (CheckSum);
}


static void update_hash(struct move2 *move){
    // update HASH_KEY incrementally:
    int i;
    int square;
    int contents;

    HASH_KEY ^= HashSTM;
    square = (move->m[0]) & 255;
    contents = ((move->m[0]) >> 8);

    HASH_KEY ^= ZobristNumbers[square][contents];
    square = (move->m[1]) & 255;
    contents = ((move->m[1]) >> 8);
    HASH_KEY ^= ZobristNumbers[square][contents];
    // captured pieces are below:
    for(i=2;i<move->l;i++){
     square = ((move->m[i]) & 255);
     contents = ((move->m[i]) >> 8);
     HASH_KEY ^= ZobristNumbers[square][contents];
                                        }
 }


static void TTableInit( unsigned int hash_mb){
 //
 //
       int j;
       U32 transsize;
       U32 size;
       U32 target;

       if ( ttable != NULL )
       free(ttable);

       target = hash_mb;
       target *= 1024 * 1024;

       j = sizeof(TEntry);
       assert(j==16); // 16 bytes
       for (size = 1; size != 0 && size <= target; size *= 2)
       ;
       size /= 2;
       assert(size>0 && size<=target);
       // allocate table
       size /= 16;
       assert(size!=0&&(size&(size-1))==0); // power of 2
       transsize = size + REHASH - 1; // HACK to avoid testing for end of table
       ttable  = ( TEntry *) malloc(transsize*16+2);
       memset(ttable,0,transsize*16+2);
       MASK = size - 1;
    }


static int rootsearch(int b[46], int alpha, int beta, int depth,int color,int search_type){
//
//
    int old_alpha;
    int bestvalue;
    int value;
    int i;
    int bestindex = 0;
    int fracDepth = 0;
    int bestfrom = -127;
    int bestto = -127;
    static int n;
    static struct move2 movelist[MAXMOVES];
    int SortVals[MAXMOVES];
    int rootnodecount[MAXMOVES];
    int oldnodecount;
    static int rootbest = 0;
    static int capture;
    U64 L_HASH_KEY;
    char string1[255];
    char string2[255];
    char string3[255];
    double elapsed = 0; // time variable

    if (*play) return 0;
    nodes++; // increase node count
    old_alpha = alpha;
    bestvalue = -32767;

    unsigned int Pieces = g_pieces[6] + g_pieces[10] + g_pieces[5] + g_pieces[9];
    if ( !EdNocaptures && Pieces <= EdPieces ){
    EdRoot[color] = EdProbe(b,color);
    if (EdRoot[color] == EdAccess::win) EdRoot[color ^ CHANGECOLOR] = EdAccess::lose;
   else if (EdRoot[color] == EdAccess::lose) EdRoot[color ^ CHANGECOLOR] = EdAccess::win;
   else EdRoot[color ^ CHANGECOLOR] = EdRoot[color];
                                                                          }

    if ( search_type == SearchShort ){
    capture = Test_Capture(b,color);
    if ( capture )
    n = Gen_Captures(b, movelist, color);
    else
    n = Gen_Moves(b,movelist, color);
    if(!n) return -MATE;
                                                        }

    if ( n == 2 ) fracDepth += 16; // extend

    for( i=0;i<MAXMOVES;i++)
    rootnodecount[i] = 0;

   L_HASH_KEY = HASH_KEY; // save HASH_KEY
   MoveToStr(movelist[0],string2); // best move found so far
   // move loop
   for ( i = 0; i<n ; i++ ){
               oldnodecount = nodes;
               domove(b,&movelist[i],color);

   if ( depth >= 13 ){
#ifdef KALLISTO
           MoveToStr(movelist[i],string1); // currently executed move
           if ( bestvalue >-32767 ) rootbest = bestvalue;
           sprintf(string3,"");
           sprintf(string3,"%2i/%2i  ",i+1,n);
           strcat(string3,string1);
           elapsed = (clock()-start)/(double)TICKS;
           if (pfSearchInfo) pfSearchInfo(rootbest, depth, elapsed > 0 ? int(nodes / elapsed / 1000) : 0, string2, string3);
#endif
                            }

               if (search_type == SearchShort || bestvalue == -32767 ){
               value = -PVSearch(b,depth-1,fracDepth,0,-beta,-alpha,color^CHANGECOLOR,NodePV,0,capture);
       if (value <= alpha){ // research
                                     old_alpha = -MATE;
        value = -PVSearch(b,depth-1,fracDepth,0,-beta,MATE,color^CHANGECOLOR,NodePV,0,capture);
       }
       else if (value >= beta){ // research
            value = -PVSearch(b,depth-1,fracDepth,0,-MATE,-alpha,color^CHANGECOLOR,NodePV,0,capture);
       }
                }else{ // other moves
           value = -PVSearch(b,depth-1,fracDepth,0,-alpha-1,-alpha,color^CHANGECOLOR,NodeCut,0,capture);
           if (value > alpha){ // research
           value = -PVSearch(b,depth-1,fracDepth,0,-beta-50,-alpha,color^CHANGECOLOR,NodePV,0,capture);
                                     }
           if (value >= beta){ // research
           value = -PVSearch(b,depth-1,fracDepth,0,-MATE,-alpha,color^CHANGECOLOR,NodePV,0,capture);
                                      }
                        }

         undomove(b,&movelist[i],color);
         // restore HASH_KEY
         HASH_KEY = L_HASH_KEY;
         if (value > bestvalue){
            bestvalue = value;
            bestindex = i;
            if (value > alpha){
            if (search_type == SearchNormal) alpha = value;
            if (value >= beta) break;
                              }
                                          }

         // save nodes taken for this move:
         if ( nodes / 100 > 1000000 )
         rootnodecount[i] = (nodes - oldnodecount)/1000;
         else
         rootnodecount[i] = (nodes - oldnodecount)/100;

            } // end move loop

           bestrootmove = movelist[bestindex];
#ifdef KALLISTO
           if ( i == n ) i--;
           MoveToStr(movelist[i],string1);
           MoveToStr(bestrootmove,string2);
           sprintf(string3,"");
           sprintf(string3,"%2i/%2i  ",i+1,n);
           strcat(string3,string1);
           elapsed = (clock()-start)/(double)TICKS;
           if (pfSearchInfo) pfSearchInfo(bestvalue, depth, elapsed > 0 ? int(nodes / elapsed / 1000) : 0, string2, string3);
#endif
                  bestfrom = ( bestrootmove.m[0] ) & 255;
                  bestto = ( bestrootmove.m[1] ) & 255;
   /* save the position in the hashtable */
                  hashstore( bestvalue,depth,old_alpha,beta, bestfrom, bestto);
                  /* and order the movelist */
                    for (i = 0;i<n;i++)
                    SortVals[i] = rootnodecount[i];
                    QuickSort( SortVals,movelist, 0,(n-1));
                    for ( i = 0;i<n;i++ ){
                    int from = ( movelist[i].m[0] ) & 255;
                    int to =  ( movelist[i].m[1] ) & 255;
                    if ( from == bestfrom && to == bestto ){
                               // swap
                                    struct move2 tmpmove = movelist[0];
                                    movelist[0] = movelist[i];
                                    movelist[i] = tmpmove;
                                    break;
                                                            }
                                          }

        return bestvalue;
}


static int compute( int b[46],int color, int time, char output[256]){
   // compute searches the best move on position b in time time.
   // it uses iterative deepening to drive the PVSearch.
   int depth;
   int i;
   int value;
   int lastvalue = 0;
   int newvalue;
   int dummy,alpha,beta;
   int bestfrom=-127;
   int bestto=-127;
   int bestindex = 0;
   int n;
   double t, elapsed=0; // time variables
   struct move2 movelist[MAXMOVES];
   struct move2 lastbest;
   char str[256];
   char pv[256];
   nodes = 0;
   sprintf(output,"KestoG engine v1.4");
   init(b);
   n = Gen_Captures(b, movelist, color);
   if(!n)
   n = Gen_Moves(b,movelist, color);
   if(!n) return -MATE;

#ifdef KALLISTO
                  if(n==1){ // only one move to do!
                  // return this move instantly
                    value = 0;
   bestrootmove=movelist[0];
   MoveToStr(bestrootmove,str);
   if (pfSearchInfo) pfSearchInfo(value,0, elapsed > 0 ? int(nodes / elapsed / 1000) : 0, str, 0);
                  movetonotation(bestrootmove,str);
                    sprintf(output,"[only move][depth %i][move %s][time %.2fs][eval %i][nodes %i]",0,str,elapsed,value,0);
                    printf("\n%s",output);
                  setbestmove(bestrootmove);
                  domove2(b,&bestrootmove,color);
   return 0;
                              }
#endif

                  //TTableInit(size);
                  searches_performed_in_game++;
                  searches_performed_in_game &= 0x3f;
                  /* clear the history heuristics: */
                  ClearHistory();
                   /* clear the killer moves: */
                  for ( i=0;i<=MAXDEPTH;i++ ){
                  killer_scores1[i] = 0;
                  killer_scores2[i] = 0;
                  killer_scores3[i] = 0;
                  killersf1[i] = 0;
                  killerst1[i] = 0;
                  killersf2[i] = 0;
                  killerst2[i] = 0;
                  killersf3[i] = 0;
                  killerst3[i] = 0;
                                             }
                  maxtime = time;
                  HASH_KEY = Position_to_Hashnumber(b,color);
                  root_depth = 3;
                  realdepth = 0;
                  start = clock();
                  newvalue = rootsearch(b,-MATE,+MATE,3,color,SearchShort);
                  alpha = newvalue - WINDOW;
                  beta = newvalue + WINDOW;
                  value = newvalue;
                  for ( depth = 5;depth<MAXDEPTH;depth+=2){
                                       lastvalue = value;
                                       lastbest = bestrootmove;
                                       repeat_search:
                                       HASH_KEY = Position_to_Hashnumber(b,color);
                                       //init(b);
                                       root_depth = depth;
                                       realdepth = 0;
                                       /* do a search with aspiration window */
                                       newvalue = rootsearch(b,alpha,beta,depth,color,SearchNormal);
                                       elapsed = (clock()-start)/(double)TICKS;
                                       // interrupt by user
                                       if (*play){
                                       value = lastvalue;
                                       bestrootmove = lastbest;
                                       depth-=2;
                                       movetonotation(bestrootmove,str);
                                       if ( nodes < 1048576 )
sprintf(output,"[done][depth %i][move %s][time %.2fs][eval %i][nodes %i]",depth,str,elapsed,value,nodes);
                                       else
sprintf(output,"[done][depth %i][move %s][time %.2fs][eval %i][nodes %.1fM]",depth,str,elapsed,value,(float)nodes/(1024*1024));
                                       printf("\n%s",output);
                                       setbestmove(bestrootmove);
                                       domove2(b,&bestrootmove,color);

#ifdef KALLISTO
   MoveToStr(bestrootmove,str);
   if (pfSearchInfo) pfSearchInfo(value, depth, elapsed > 0 ? int(nodes / elapsed / 1000) : 0, str, 0);
   return value;
#endif
                                                     }

         /* check if aspiration window holds */
         if( newvalue <= alpha ){
         alpha = newvalue - 4 * WINDOW;
         if (alpha < -MATE) alpha = -MATE;
         goto repeat_search;
                                }
         else
         {
         if ( newvalue >= beta ){
         beta = newvalue + 4 * WINDOW;
         if (beta > MATE) beta = MATE;
         goto repeat_search;
                                }
else{       /* now iteration complete, set aspiration window for the next iteration */
           value = newvalue;
           alpha = newvalue - WINDOW;
           if (alpha < -MATE) alpha = -MATE;
           beta = newvalue + WINDOW;
           if (beta > MATE) beta = MATE;
    }
                                 }

      bestfrom = -127;
      bestto = -127;

      // get best move from hashtable:
      hashretrieve(MAXDEPTH, &dummy, &alpha, &beta, &bestfrom, &bestto, &dummy, &dummy);
     // we always store in the last call to PVSearch, so there MUST be a move here!
for ( i = 0; i < n; i++ ){
if ( (( movelist[i].m[0] ) & 255) == bestfrom &&  (( movelist[i].m[1] ) & 255) == bestto ){
             bestindex = i;
             break;                                                                                                                }
             }
      movetonotation(movelist[bestindex],str);
      if ( nodes < 1048576 )
sprintf(output,"[thinking][depth %i][move %s][time %.2fs][eval %i][nodes %i]",depth,str,elapsed,value,nodes);
      else
sprintf(output,"[thinking][depth %i][move %s][time %.2fs][eval %i][nodes %.1fM]",depth,str,elapsed,value,(float)nodes/(1024*1024));
      printf("\n%s",output);
                      // break conditions:
                      // time elapsed
                      t = clock();
/* don't bother going deeper if we've already used 70% of our time since we likely won't finish */
if ( 1000*((t-start) / (double)TICKS) > (0.7*maxtime) && root_depth > 5 ) break;

      // found a win
if ( value >= MATE-root_depth-1 && root_depth > 5 ) break;

     /* reset the killer scores (we can keep the moves for move ordering for
     now, but the scores may not be accurate at higher depths, so we need
     to reset them): */
                                    for ( i=0;i<=MAXDEPTH;i++ ){
                                    killer_scores1[i] = 0;
                                    killer_scores2[i] = 0;
                                    killer_scores3[i] = 0;
                                                               }
                  };  // end for iterative deepening loop

   sprintf(pv,"");
                  retrievepv(b,pv,color);
                  movetonotation(movelist[bestindex],str);
                  if ( nodes < 1048576 )
sprintf(output,"[done][depth %i][move %s][time %.2fs][eval %i][nodes %i][pv %s]",depth,str,elapsed,value,nodes,pv);
   else
sprintf(output,"[done][depth %i][move %s][time %.2fs][eval %i][nodes %.1fM][pv %s]",depth,str,elapsed,value,(float)nodes/(1024*1024),pv);
   printf("\n%s",output);
                  setbestmove(bestrootmove);
                  domove2(b, &bestrootmove,color);
       return value;
   }


                   /* main search function */
                   /* Principal Variation Search */
static int PVSearch(int b[46],int depth,int fracDepth,int trunc,int alpha,int beta,int color,int node_type,int  iid_search,int xcapture){
   int value;
   int i,n;
   int maxvalue;
   int bestindex = 0;
   int played_nb;
   int capture;
   int opt_value = MATE;
   struct move2 movelist[MAXMOVES];

                  /* time check */
                  if ( !( ++nodes & 0x1fff ) )
                  if ( 1000*(clock()-start)/(double)TICKS>(2*maxtime) && root_depth > 5 ){(*play) = 1;}

                    /* return if calculation interrupt */
                  if (*play) return 0;

                  /* stop search if maximal search depth is reached */
   if ( realdepth >= MAXDEPTH )
                  return eval(b,color,alpha,beta,node_type == NodePV);

                  /* check for database use */
                  unsigned int Pieces = g_pieces[6] + g_pieces[10] + g_pieces[5] + g_pieces[9];
                  if ( !EdNocaptures && Pieces <= EdPieces )
   {
      int res = EdProbe(b,color);
      if (res != EdAccess::not_found)
      {
         if (res != EdRoot[color] || !Reversible[realdepth - 1])
         {
      if (res == EdAccess::win)  return 9000 - 100*Pieces; // return ED_WIN - realdepth;
      if (res == EdAccess::lose) return -9000 + 100*Pieces; // return -ED_WIN + realdepth;
      if (res == EdAccess::draw) return 0;
      MessageBox(0, "unknown value from EdAccess", "Error", 0);
         }
      }
   }

                  /* mate distance pruning */

                  value = MATE-realdepth-1;
                  if (value < beta){
                  beta = value;
                  if (value <= alpha) return value;
                                   }

                  /* truncations */
                  trunc = 0; // disabled
                  /* horizon? */
                  if ( depth <= 0 ){
                  capture = Test_Capture(b,color);
                  if ( capture == 0 ){
                  int staticValue = eval(b,color,alpha,beta,node_type == NodePV);
                  if ( staticValue >= beta ) // selective deviation
                  return staticValue;
                  if( staticValue > alpha )
                  alpha = staticValue;
                  n = Gen_Proms(b,movelist,color);
                  if ( n==0) return staticValue;

                  maxvalue = staticValue; // preset maxvalue

                  // move loop
                  for( i=0; i < n; i++ ){
                  domove2(b,&movelist[i],color);
                  if ( g_pieces[(color^3|MAN)] == 0 && g_pieces[(color^3|KING)] == 0 )
                  value = MATE - realdepth;
                  else{
// ******************* recursion*********************************************************************
value = -PVSearch(b,depth-1,fracDepth,trunc,-beta, -alpha,color^CHANGECOLOR,NODE_OPP(node_type),0,0);
// **************************************************************************************************
                      }
                   undomove(b,&movelist[i],color);
                   if ( value > maxvalue ){
                   maxvalue = value;
                   if ( value > alpha ){
                   if ( value >= beta ) break;
                   alpha = value;
                                               }
                                                       }
                                               } // for
                   return maxvalue;
                               }
                  else{
                  n = Gen_Captures(b,movelist,color);
                  int SortVals[MAXMOVES];
                  // sort captures list
                  for ( i=0;i<n;i++ ){
                  if (is_promotion(&movelist[i])) SortVals[i] = 10000;
                  else
                  SortVals[i] = 0;
                  for ( int j=2;j<movelist[i].l;j++ )
                  SortVals[i] += MVALUE[(movelist[i].m[j]) >> 8];
                                             }
                            maxvalue = -32767;    // preset maxvalue
                            // move loop
                            while (pick_next_move( &i,SortVals, n ) != 0){
                            domove2(b,&movelist[i],color);
                            if ( g_pieces[(color^3|MAN)] == 0 && g_pieces[(color^3|KING)] == 0 )
                            value = MATE-realdepth;
                            else
                            {
                            /****************** recursion***************************/
value=-PVSearch(b,depth-1,fracDepth,trunc,-beta, -alpha, color^CHANGECOLOR,NODE_OPP(node_type),0,1);
                            /*******************************************************/
                            }
                            undomove(b,&movelist[i],color);
                            if (value > maxvalue){
                            maxvalue = value;
                            if (value > alpha){
                            if (value >= beta) break;
                            alpha = value;
                                              }
                                                  }
                                                  }
                          return maxvalue;
                                   }
                  } // horizon ?

   int bestfrom = -127; // best move's from square
   int bestto = -127;   // best move's to square
   int try_mcp = 127;   // try multi cut prune,initially enabled
                  // hashlookup
   if ( !iid_search ){
   if ( hashretrieve(depth,&value,&alpha,&beta,&bestfrom,&bestto,&try_mcp,node_type == NodePV))
   return value;
                                          }

                  // check if the side to move has a capture
                  capture = Test_Capture(b,color);   // is there a capture for the side to move?
                  /* check for database use */
                  if ( EdNocaptures && Pieces <= EdPieces && !capture ){
      int res = EdProbe(b,color);
      if (res != EdAccess::not_found)
      {
         if (res != EdRoot[color] || !Reversible[realdepth - 1])
         {
      if (res == EdAccess::win)  return 9000 - 100*Pieces; // return ED_WIN - realdepth;
      if (res == EdAccess::lose) return -9000 + 100*Pieces; // return -ED_WIN + realdepth;
      if (res == EdAccess::draw) return 0;
      MessageBox(0, "unknown value from EdAccess", "Error", 0);
         }
      }
                                                                                                            }

                  if ( capture )
                  n = Gen_Captures(b,movelist,color);
   else
   n = Gen_Moves(b,movelist,color);

                  // if we have no move:
   if ( n == 0 )
   return (realdepth-MATE); // minus sign will be negated in negamax framework

                 U64 L_HASH_KEY = HASH_KEY;  // local variable  for saving position's HASH_KEY

                 if ( n == 1 ){
                     fracDepth += 20; // only 1 move
                     if ( fracDepth > 31 ){
                     fracDepth-=32;
                     depth++;
                                          }
                               }

                  if ( n == 2 ){
                      fracDepth += 8; // only 2 moves
                      if ( fracDepth > 31 ){
                      fracDepth-=32;
                      depth++;
                                           }
                               }

                 /* enhanced transposition cutoffs: do every move and check if the resulting
                 position is in the hashtable. if yes, check if the value there leads to a cutoff
                 if yes, we don't have to search */
                 if ( depth >= ETCDEPTH && beta > -1500 ){
                        int dummy;
                        int ETCalpha;
                        int ETCbeta;
                        int ETCvalue;
                        for( i=min(n-1,10);i>=0;i-- ){
                                    // do move
                                    // domove(b,&movelist[i] );
                                    assert(i>=0);
                                    update_hash(&movelist[i]);
                                        /* do the ETC lookup:
                                        with reduced depth and changed color */
                                    ETCalpha = -beta;
                                    ETCbeta = -alpha;
                                    ETCvalue = -MATE;
                                    if (hashretrieve(depth-1,&ETCvalue,&ETCalpha,&ETCbeta,&dummy,&dummy,&dummy,&dummy)){
                                        /* if one of the values we find is > beta we quit! */
                                        if ( ( -ETCvalue) >= beta ){
                                        /* before we quit: restore all stuff */
                                        // undomove(b,&movelist[i] );
                                        HASH_KEY = L_HASH_KEY;
                                        return (-ETCvalue); // or return beta ?
                                                                    }
                                                                                 } // if
                                       // undomove(b,&movelist[i] );
                                       HASH_KEY = L_HASH_KEY;
                                                     } // for
                                                      } // ETC

            // Internal iterative deepening:

if ( bestfrom == -127 && bestto == -127 ){
if ( n > 1 ){
if ( opt_value == MATE ) opt_value = matval(b,color);
if (( node_type==NodePV && depth>=4 ) || ( node_type!=NodePV && depth>=8 && opt_value>alpha-100)){
                    bestindex = 0;
                    int tempalpha = alpha;
                    int new_depth = node_type == NodePV ? depth -2 : min(depth-2,depth/2);
                    maxvalue = -32767;    // preset maxvalue
                    for( i = 0; i < n; i++){
                    domove(b,&movelist[i],color);
                    if ( g_pieces[(color^3|MAN)] == 0 && g_pieces[(color^3|KING)] == 0 )
                    value = MATE-realdepth;
                    else
value=-PVSearch(b,new_depth-1,fracDepth,trunc,-beta,-tempalpha,color^3,NODE_OPP(node_type),1,capture);
                    undomove(b,&movelist[i],color);
                    // restore HASH_KEY
                    HASH_KEY = L_HASH_KEY;
                    if ( value > maxvalue ){
                    maxvalue = value;
                    bestindex = i;
                    if ( value > tempalpha ) tempalpha = value;
                    if ( value >= beta ) break;
                                                      }
                                                   } // for

                    bestfrom = ( movelist[bestindex].m[0] ) & 255;
                    bestto = ( movelist[bestindex].m[1] ) & 255;
                    bestindex = 0;

                    if ( maxvalue < alpha && maxvalue <= -HASHMATE )
                    return maxvalue;
                    if ( maxvalue >= beta && maxvalue >= HASHMATE )
                    return maxvalue;

                                         }
                                 }
                          }

                         int p,q;
                         int SortVals[MAXMOVES];
                         // sort move list:
                         if ( !capture ){
                         int testhash = (bestfrom == -127 || bestto == -127) ? 0:1;
                         for ( i = 0; i < n; i++ ){ // loop
                         if ( testhash ){
                         p = ( movelist[i].m[0] ) & 255; // from square
                         q =  ( movelist[i].m[1] ) & 255; // to square
if ( ( p == bestfrom  ) && ( q == bestto ) ){ SortVals[i] = 1000;testhash = 0;continue;}
                                             }
                         if ( is_promotion(&movelist[i] )) {SortVals[i] = 800;continue;}
                         p = ( movelist[i].m[0] ) & 255; // from square
                         q =  ( movelist[i].m[1] ) & 255; // to square
if ( ( p == killersf1[realdepth] ) && ( q == killerst1[realdepth] ) ){SortVals[i] = 110;continue;}
if ( ( p == killersf2[realdepth] ) && ( q == killerst2[realdepth] ) ){SortVals[i] = 108;continue;}
if ( ( p == killersf3[realdepth] ) && ( q == killerst3[realdepth] ) ){SortVals[i] = 106;continue;}
                         SortVals[i] = -16384;
                         SortVals[i] += History[p][q];
                                                          } // loop
                                            }
                          else{
                                // assert(capture);
int testhash = (bestfrom == -127 || bestto == -127) ? 0:1; // test for hashmove flag
                                for ( i = 0; i < n; i++ ){ // loop
                                      if ( testhash ){
                                      p = ( movelist[i].m[0] ) & 255; // from square
                                      q =  ( movelist[i].m[1] ) & 255; // to square
if ( ( p == bestfrom  ) && ( q == bestto ) ){ SortVals[i] = 100000;testhash = 0;continue;}
                                                     }
                                      if ( is_promotion(&movelist[i] )) SortVals[i] = 10000;
                                      else
                                      SortVals[i] = 0;
                                      for ( int j=2;j<movelist[i].l;j++ )
                                      SortVals[i]+=MVALUE[(movelist[i].m[j]) >> 8];
                                                                 } // loop
                                }

       // ---------------------------------------------------------------------
       // forward pruning, so-called mc-prune or multicut-prune
       // parameters:
       // m - number of moves to look at when checking for mc-prune
       // c - is the number of cutoffs to cause an mc-prune
       // 1 <= c <= m
       // r - is the search depth reduction for mc-prune searches
       // here m = 8,c = 1, r = adaptive
       // do not works in endgames
       // applied only at expected Cut nodes
       // ----------------------------------------------------------------------

       int pieces = g_pieces[(color|MAN)] + g_pieces[(color|KING)];
       if ( xcapture == 0 && realdepth > 3 && pieces >= 5 && try_mcp && depth >= 4 ){
       if ( n >= 8 && node_type == NodeCut && beta > -1500 ){
       if ( depth == 4 || matval(b,color) >= beta ){
       int mcp_depth = 4;
       QuickSort( SortVals,movelist, 0,(n-1));
       int node_type1 = NodeCut;
       for ( i=0;i<8;i++){
                                    // do move
                                    domove(b,&movelist[i],color);
                                    if ( g_pieces[(color^3|MAN)] == 0 && g_pieces[(color^3|KING)] == 0 )
                                    value = MATE-realdepth;
                                    else
value=-PVSearch(b,depth-mcp_depth,fracDepth,trunc,-beta,-alpha,color^3,NODE_OPP(node_type1),0,capture);
                                     //  undo move
                                    undomove(b,&movelist[i],color);
                                    // restore HASH_KEY
                                    HASH_KEY = L_HASH_KEY;
                                    if ( value >= beta ){
                                    bestfrom = ( movelist[i].m[0] ) & 255;
                                    bestto = ( movelist[i].m[1] ) & 255;
                                    hashstore( value, depth,alpha,beta,bestfrom,bestto);
                                    return value;
                                                         }
                                    node_type1 = NodeAll;
                                                } // for
                                            } // if
                                        } // if
                                      } // if

                  int margin;
                  int reduced;
                  played_nb = 0;
                  int played[MAXMOVES];
                  maxvalue = -32767; // preset maxvalue
                  int oldalpha = alpha;
                  //if ( depth >= 4 ){
                  if (realdepth < MAXDEPTH - 4) {
                  killersf1[realdepth+4]=0;
                  killerst1[realdepth+4]=0;
                  killersf2[realdepth+4]=0;
                  killerst2[realdepth+4]=0;
                  killersf3[realdepth+4]=0;
                  killerst3[realdepth+4]=0;
                                   }

// Ok,let's look now at all moves and pick one with the biggest value

 while (pick_next_move( &i,SortVals, n ) != 0){
   played[played_nb++] = i;
if (realdepth > 2 && node_type != NodePV && played_nb >= 2 + depth && !move_is_dangerous(b,&movelist[i])){
if (ok_to_prune((movelist[i].m[0])&255,(movelist[i].m[1])&255,depth))
continue;
    }
       reduced = 0;
       int new_depth = depth - 1;

   if ( node_type != NodePV && depth >= 3 && played_nb >= 3 ){
   if (ok_to_reduce((movelist[i].m[0])&255, (movelist[i].m[1])&255)){
   if ( !move_is_dangerous(b,&movelist[i])){
   new_depth--;
   reduced = 1;
            }
          }
      }

      // Rebel's style reductions:
   if ( node_type != NodePV && new_depth > 2 && new_depth <= 31 && !capture && !reduced ){
   if ( realdepth > 2 && !ISLOSSSCORE(alpha) && !ISWINSCORE(alpha) && !is_promotion(&movelist[i] ) ){
   if (opt_value == MATE) opt_value = matval(b,color);
   value = opt_value + red_margin[new_depth];
   if ( value < alpha ) { new_depth--;reduced = 1;}
                                                                        }
                                                                      }

       // domove
       domove(b,&movelist[i],color);
       // forward pruning
if (depth>=4 && node_type!=NodePV && n>1 && realdepth>=2 && !capture && !is_promotion(&movelist[i])){
if (!reduced && !move_is_dangerous(b,&movelist[i])){
value=-PVSearch(b,depth-4,fracDepth+4,trunc,-alpha+32-1, -alpha+32, color^3,NODE_OPP(node_type),0,capture);
if (( value < alpha-32 ) && ( value > -HASHMATE ))
goto skip_search;
    }
 }
                             if ( g_pieces[(color^3|MAN)] == 0 && g_pieces[(color^3|KING)] == 0 )
                             value = MATE-realdepth;
                             else
                                    {
 /******************* recursion***********************************************/
                                    if ( node_type != NodePV || maxvalue == -32767 ){ // first move
value = -PVSearch(b, new_depth,fracDepth,trunc,-beta, -alpha, color^CHANGECOLOR,NODE_OPP(node_type),0,capture);
                                    }else{ // other moves
value = -PVSearch(b, new_depth,fracDepth,trunc,-alpha-1, -alpha, color^CHANGECOLOR, NodeCut,0,capture);
      if (  value > alpha )
      value = -PVSearch(b, new_depth,fracDepth,trunc,-beta, -alpha, color^CHANGECOLOR, NodePV,0,capture);
                                            }

     //  history pruning re-search
     if ( value > alpha && reduced ){
     new_depth++;
value=-PVSearch(b, new_depth,fracDepth,trunc,-beta, -alpha, color^CHANGECOLOR,NODE_OPP(node_type),0,capture);
                                     }
 /*******************************************************************************/
                                     }

      skip_search:
      undomove(b,&movelist[i],color);
      // restore HASH_KEY

                                    HASH_KEY = L_HASH_KEY;
      // update best value so far
      // and set alpha and beta bounds
                                    if ( value > maxvalue ){
                                    maxvalue = value;
                                    if ( value > alpha ){
                                    alpha = value;
                                    bestindex = i;
                                    if ( value >= beta ) break;
                                                         }
                                                           }
     // if we were supposed to FH but did not ...
     if ( node_type == NodeCut ) node_type = NodeAll;

     }; // end main recursive loop of forallmoves
         if (  alpha > oldalpha ){
         bestfrom = ( movelist[bestindex].m[0] ) & 255;
         bestto = ( movelist[bestindex].m[1] ) & 255;
         if ( maxvalue >= beta ){
         if ( !capture && !is_promotion(&movelist[bestindex] ) ){
         for (i = 0; i < played_nb - 1; i++){
         int j = played[i];
         if ( !is_promotion(&movelist[j] ) )
         history_bad( ((movelist[j].m[0]) & 255),((movelist[j].m[1]) & 255));
                                            }
         good_move(bestfrom,bestto,realdepth);
         history_good(bestfrom,bestto);
                                 }
                                  }
                                }
        hashstore( maxvalue, depth,oldalpha,beta, bestfrom, bestto);
        return maxvalue;
       }


static void hashstore(int value, int depth, int alpha, int beta, int best_from, int best_to){
     //
     //
                  TEntry *entry;
                  TEntry *replace;
                  U32 lock;
                  register int i;
                  replace = NULL;

                  entry = replace = ttable + (U32)(HASH_KEY&MASK);
                  lock = (U32)(HASH_KEY >> 32);

                  for( i=0;i < REHASH;i++){
                         if ((entry+i)->m_lock == lock ){  // hash hit => update existing entry

                  if ( best_from == -127 && best_to == -127 ){
                  best_from = (entry+i)->m_best_from;
                  best_to = (entry+i)->m_best_to;
                                                             }                                                                                     if(value>=beta){
                                       (entry+i)->m_valuetype = LOWER;
                                       (entry+i)->m_best_from = best_from;
                                       (entry+i)->m_best_to = best_to;
                                       (entry+i)->m_lock = lock;
                                       (entry+i)->m_depth = depth;
                                       (entry+i)->m_value = value;
                                       (entry+i)->m_age = searches_performed_in_game;
                                       return;
                                                           }

                                    if(value<=alpha){
                                       (entry+i)->m_valuetype = UPPER;
                                       (entry+i)->m_best_from = best_from;
                                       (entry+i)->m_best_to = best_to;
                                       (entry+i)->m_lock = lock;
                                       (entry+i)->m_depth = depth;
                                       (entry+i)->m_value = value;
                                       (entry+i)->m_age = searches_performed_in_game;
                                       return;
                                                              }

                                       (entry+i)->m_valuetype = EXACT;
                                       (entry+i)->m_best_from = best_from;
                                       (entry+i)->m_best_to = best_to;
                                       (entry+i)->m_lock = lock;
                                       (entry+i)->m_depth = depth;
                                       (entry+i)->m_value = value;
                                       (entry+i)->m_age = searches_performed_in_game;
                                       return;
                                                               } // hash hit

                                       if(replace->m_age == searches_performed_in_game){
if((entry+i)->m_age != searches_performed_in_game || (entry+i)->m_depth < replace->m_depth)
replace = (entry+i);
                                                                                        }
else if((entry+i)->m_age != searches_performed_in_game && (entry+i)->m_depth < replace->m_depth)
                                      replace = (entry+i);
                                 } // rehash

                                      if(value>=beta){
                                       // save:
                                       replace->m_valuetype = LOWER;
                                       replace->m_best_from = best_from;
                                       replace->m_best_to = best_to;
                                       replace->m_lock = lock;
                                       replace->m_depth = depth;
                                       replace->m_value = value;
                                       replace->m_age = searches_performed_in_game;
                                       return;
      }

                                   if(value<=alpha){
                                       // save:
                                       replace->m_valuetype = UPPER;
                                       replace->m_best_from = best_from;
                                       replace->m_best_to = best_to;
                                       replace->m_lock = lock;
                                       replace->m_depth = depth;
                                       replace->m_value = value;
                                       replace->m_age = searches_performed_in_game;
                                       return;
      }

                                       // save:
                                       replace->m_valuetype = EXACT;
                                       replace->m_best_from = best_from;
                                       replace->m_best_to = best_to;
                                       replace->m_lock = lock;
                                       replace->m_depth = depth;
                                       replace->m_value = value;
                                       replace->m_age = searches_performed_in_game;
   }


static int hashretrieve(int depth,int *value,int *alpha,int *beta,int *best_from,int *best_to,int *try_mcp,bool in_pv){
   //
   //
                  TEntry *entry,*beste;
                  register int i;
                  U32 lock;
                  int found = 0;

                  entry = ttable + (U32)(HASH_KEY & MASK);
                  lock = (U32)(HASH_KEY >> 32);

                  for( i = 0; i < REHASH; i++ ){
                  if ( (entry+i)->m_lock == lock ){
                  found = 1;
                  beste = entry+i;
                                                   }
                                               }

                   if ( !found ) return (0);

                   // set best move

   *best_from = beste->m_best_from;
   *best_to = beste->m_best_to;

                  // forward pruning switch
                  if ( beste->m_valuetype == UPPER )
                  if ( beste->m_depth >= depth - 2 - 1 )
                  if ( beste->m_value < *beta ) *try_mcp = 0;

                  // check if depth this time round is higher
   if( depth > beste->m_depth ){
      // we are searching with a higher remaining depth than what is in the hashtable.
      // all we can do is set the best move for move ordering
                                     return 0;
      }

  // we have sufficient depth in the hashtable to possibly cause a cutoff.
  // if we have an exact value, we don't need to search for a new value.

   if ( beste->m_valuetype == EXACT){
   if ( !in_pv ){
                  *value = beste->m_value;
                   return 1;
                }
                   // if in_pv use only best move from hash
                   // *value = entry->m_value;
                   return 0;
                                    }

   // if we have a lower bound, we might either get a cutoff or raise alpha.
   if ( beste->m_valuetype == LOWER){
   // the value stored in the hashtable is a lower bound, so it's useful
   if ( beste->m_value >= *beta ){
                  // value > beta: we can cutoff!
                  *value = beste->m_value;
                  return 1;
                                  }
                  return 0;
                                     }
   // if we have an upper bound, we can either get a cutoff or lower beta.
   if ( beste->m_valuetype == UPPER){
   // the value stored in the hashtable is an upper bound, so it's useful
                  if( beste->m_value <= *alpha ){
                  // value < alpha: we can cutoff!
                  *value = beste->m_value;
                   return 1;
                                 }
                   return 0;
                                                                   }
        return 0;
                      }


static void retrievepv( int b[46], char *pv, int color){
   // gets the pv from the hashtable
   // get a pv string:
   int n;
   int i;
   int bestfrom;
   int bestto;
   int bestindex = 0;
   struct move2 movelist[MAXMOVES];
   int dummy, alpha, beta;
   char pvmove[256];
   int count = 0;
                  int copy[46];
                  // original board b[46] needs not to be changed
                  for ( i=0;i<46;i++ )
                      copy[i] = b[i];

   bestfrom = -127;
   bestto = -127;
   sprintf(pv,"");
   init(copy);
   HASH_KEY = Position_to_Hashnumber(copy,color);
   hashretrieve(MAXDEPTH, &dummy, &alpha, &beta, &bestfrom, &bestto, 0,false);
   while(bestfrom != -127 && bestto != -127 && count<10){
      n = Gen_Captures( copy, movelist, color);
      if(!n)
      n = Gen_Moves( copy, movelist, color);
                                    if (!n) return;
                                    for ( i = 0; i < n ; i++ ){
if ( (( movelist[i].m[0] ) & 255) == bestfrom &&  (( movelist[i].m[1] ) & 255) == bestto ){
bestindex = i;
break;
                                                }
                                        }

      movetonotation(movelist[bestindex],pvmove);
      domove2( copy, &movelist[bestindex],color );
      strcat(pv," ");
      strcat(pv,pvmove);
      color = color^CHANGECOLOR;
      HASH_KEY = Position_to_Hashnumber(copy,color);
      // look up next move
      bestfrom = -127;
      bestto = -127;
      hashretrieve(MAXDEPTH, &dummy, &alpha, &beta, &bestfrom, &bestto, 0,false);
      count++;
      }
             }


static __inline
int weight( ){
     // depth is important when using history heuristics
     // move searched with bigger depth is more promising
     int score;
     score = max( 0, root_depth - realdepth + 2 );
     score *= score;
     return score;
 }


static void good_move( int from,int to,int ply){
// update killers & history
//
   int p,q;

      /* first, check whether it matches one of the known killers */
      if ( from == killersf1[ply] && to == killerst1[ply])
   {
     killer_scores1[ply]++; // increase popularity counter for 1st killer
   }
      else if ( from == killersf2[ply] && to == killerst2[ply])
   {
     killer_scores2[ply]++; // increase popularity counter for 2nd killer
                    // if 2nd killer more popular than 1st swap them
     if (killer_scores2[ply] > killer_scores1[ply]){
                        SWAPINT (killersf1[ply],killersf2[ply]);
                        SWAPINT (killerst1[ply],killerst2[ply]);
                        SWAPINT (killer_scores1[ply],killer_scores2[ply]);
       }
   }
      else if ( from == killersf3[ply] && to == killerst3[ply])
   {
     killer_scores3[ply]++; // increase popularity counter for 3rd killer
                    // if 3rd killer more popular than 2nd swap them
     if (killer_scores3[ply] > killer_scores2[ply]){
                        SWAPINT (killersf2[ply],killersf3[ply]);
                        SWAPINT (killerst2[ply],killerst3[ply]);
                        SWAPINT (killer_scores2[ply],killer_scores3[ply]);
                       }
   }
   /* if not, replace killer3 */
   else
   {
     killer_scores3[ply] = 1;
     killersf3[ply] = from;
     killerst3[ply] = to;
   }

        History[from][to] += weight();
        if ( History[from][to] >= MAXHIST ){
               for (p=5;p<=40;p++)
                 for (q=5;q<=40;q++)
                     History[p][q] = ( History[p][q] +1 ) / 2;
                                                                   }
}


static __inline void history_bad( int from,int to ){
//
  HistTot[from][to]++;
}


static __inline void history_good( int from,int to ){
//
  HistHit[from][to]++;
}


static __inline bool ok_to_reduce( int from,int to ){
//
  bool value;
     value = (HistHit[from][to] <= 4 * HistTot[from][to]);
     return value;
}


static __inline bool ok_to_prune( int from,int to,int depth ){
//
  bool value;
     value = (depth * HistHit[from][to] < HistTot[from][to]);
     return value;
}


static void ClearHistory(void){
     // clear previous History before each new search
     int p,q;
     for ( p=5;p<=40;p++)
       for ( q=5;q<=40;q++){
         History[p][q] = 0;
         HistHit[p][q] = 1;
         HistTot[p][q] = 1;
           }
 }


static void QuickSort( int SortVals[MAXMOVES],struct move2 movelist[MAXMOVES], int inf, int sup){
    // quick sort algorithm used to sort movelist
        int pivot;
        register int i,j;
        int swap;
        struct move2 temp;
        i = inf;
        j = sup;
        pivot = SortVals[(i+j)/2];
   do {
      while (SortVals[i] > pivot) i++;
      while (SortVals[j] < pivot) j--;
      if (i<j) {
         swap = SortVals[i];
         SortVals[i] = SortVals[j];
         SortVals[j] = swap;
         temp = movelist[i];
         movelist[i] = movelist[j];
         movelist[j] = temp;
      }
      if (i<=j) {
         i++;
         j--;
      }
   } while (i<=j);
    if (inf<j) QuickSort(SortVals,movelist,inf,j); // recurse
    if (i<sup) QuickSort(SortVals,movelist,i,sup); // recurse
}


static int matval(int b[46], int color ){
// returns material balance
    int nbm = g_pieces[6];
    int nbk = g_pieces[10];
    int nwm = g_pieces[5];
    int nwk = g_pieces[9];

    if ( nbm == 0 && nbk == 0 )  return ( color == BLACK ? ( realdepth-MATE):MATE-realdepth);
    if ( nwm == 0 && nwk == 0 ) return( color == WHITE  ? ( realdepth-MATE):MATE-realdepth);

    int phase = get_phase();
    int v1,v2;
    int eval;
    if ( phase == ENDGAME ){

           v1=100*nbm+300*nbk;
           v2=100*nwm+300*nwk;
           eval = v1-v2;
           int White = nwm + nwk;
           int Black = nbm + nbk;

           if ( nbk > 0 && ( White < (2+nbk)) && (eval < 0)) return 0;
           if ( nwk > 0 && ( Black < (2+nwk)) && (eval > 0)) return 0;

                                            int WGL = 0;
                                            int BGL = 0;
                                            for ( int i=5;i<=40;i+=5 ){
                                            if ( b[i] ){
                                            if ( b[i] & MAN ){
                                            WGL=0;
                                            BGL=0;
                                            break;
                                                                        }
                                            if ( b[i] == WHT_KNG ) WGL = 1;
                                            else
                                            if ( b[i] == BLK_KNG ) BGL = 1;
                                                                          }
                                                                                }

                   // surely winning advantage:
                   if ( White == 1 && nwm == 1 && Black >= 4 ) eval = eval + (eval>>1);
                   if ( Black == 1 && nbm == 1 && White>= 4 ) eval = eval + (eval>>1);

                  // scaling down
                  if ( nbk > 0 && eval < 0 ) eval = eval >> 1;
                  if ( nwk > 0 && eval > 0 ) eval = eval >> 1;

                   if ( nbk == 1 && BGL && !WGL && White <= 3 )
                   if ( Black <= 2 || eval < 500 )
                   return (0);

                   if ( nwk == 1 && WGL && !BGL && Black <= 3 )
                   if ( White <= 2 || eval > -500 )
                   return (0);

                // negamax formulation requires this:
                if(color == BLACK){
                eval++;
                return (eval);
                                                }
               else{
               eval--;
               return (-eval);
                     }
                              }  // ENDGAME

           v1 = 100*nbm+250*nbk;
           v2 = 100*nwm+250*nwk;
           eval = v1-v2;

           eval += (200*(v1-v2))/(v1+v2);      /*favor exchanges if in material plus*/

           // king's balance
           if ( nbk != nwk){
           if ( nwk == 0 && nbm >= nwm-2 )
           eval += 200;
           else
           if ( nbk == 0 && nwm >= nbm-2 )
           eval -= 200;
                                     }

            // scaling down
            if ( nbk > 0 && eval < 0 ) eval = ((3*eval) >> 2);
            if ( nwk > 0 && eval > 0 ) eval = ((3*eval) >>2);

           // negamax formulation requires this:
           if(color == BLACK){
           eval+=2; // turn
           return (eval);
                                           }
           else{
           eval-=2; // turn
           return (-eval);
                }

    }


static int __inline is_promotion(struct move2 *move){
//
//
        if ( ((move->m[0]) >> 8) == ((move->m[1]) >> 8) )
           return 0;
           return 1;
}


static __inline int get_phase(){
//
//
    const int mat = g_pieces[6]+g_pieces[5]+g_pieces[9]+g_pieces[10];
    if ( mat > 19 ){ // 20-24
            return OPENING;
                          }
    else if ( mat > 8 ){ // 9-19
           return MIDGAME;
                              }
    else{ // 0-8
           return ENDGAME;
          }
}


static int pick_next_move( int *marker,int SortVals[MAXMOVES],int n ){
     // a function to give pick the top move order, one at a time on each call.
     // Will return 1 while there are still moves left, 0 after all moves
     // have been used

          int i,best = -32767;

          *marker = -32767;

          for ( i=0; i < n; i++ ){
                if ( SortVals[i] > best ){
                     *marker = i;
                      best = SortVals[i];
                                                   }
                                          }

          if ( *marker > -32767 ){
                 SortVals[*marker] = -32767;
                 return 1;
                                              }
          else
                 return 0;
     }


static void Sort( int start,int num, int SortVals[],struct move2 movelist[] ){
        // do a linear search through the current ply's movelist starting at start
        // and swap the best one with start

        int best = SortVals[start];
        int besti = start; // best index
        int i,j;

        for( i = start+1; i < num; i++){
             if( SortVals[i]  > best){
                  best = SortVals[i];
                  besti = i;
                                     }
                                       }

        if ( besti > start ){ // swap
        struct move2 m = movelist[start];
        movelist[start] = movelist[besti];
        movelist[besti] = m;
        j = SortVals[start];
        SortVals[start] = SortVals[besti];
        SortVals[besti] = j;
                                }
          }

          
static bool move_is_dangerous( int b[46],struct move2 *move )
{
//
int to,piece;

to = (move->m[1]) & 255;
piece = b[(move->m[0]) & 255];
  
if (( piece == BLK_MAN ) && (to > 31) && (to < 37)){
if ( b[to+5] == FREE || b[to+5] == OCCUPIED )
if ( b[to+4] == FREE || b[to+4] == OCCUPIED )
return true;
                                                   }

if (( piece == WHT_MAN ) && (to < 14) && (to > 9)){
if ( b[to-4] == FREE || b[to-4] == OCCUPIED )
if ( b[to-5] == FREE || b[to-5] == OCCUPIED )
return true;
                                                  }
return false;
}


static void  init( int b[46] ){
int i;
int color;
// int index;

num_wpieces = 0;
num_bpieces = 0;

g_pieces[5] = 0; //
g_pieces[6] = 0; //
g_pieces[9] = 0; //
g_pieces[10] = 0; //

/*
     for ( i=0;i<=31;i++ )
        plist[i] = 0;

     for ( i=5;i<=40;i++ ){
        if  ( ( b[i] != OCCUPIED ) && ( b[i] != FREE ) ){
               index=pos[i];
               plist[index] = i;
                                        }
                                     }
*/

  for ( i=0;i<=40;i++ ){
       indices[i] = 0;
                       }

  for ( i=0;i<=12;i++ ){
       p_list[1][0] = 0;
       p_list[2][0] = 0;
                       }

  for ( i=5;i<=40;i++ ){
       if  ( ( b[i] != OCCUPIED ) && ( b[i] != FREE ) ){
       g_pieces[b[i]]++;
       color = ( b[i] & WHITE ) ? WHITE:BLACK;
       if ( color == WHITE ){
       num_wpieces += 1;
       p_list[1][num_wpieces] = i;
       indices[i] = num_wpieces;
                                          }
       else{
       num_bpieces += 1;
       p_list[2][num_bpieces] = i;
       indices[i] = num_bpieces;
             }
                                                                                    }
       else
       indices[i] = 0;
                                 }
}


static void Perft(int b[46],int color,int depth){
    int capture;
    int i;
    int n;
    struct move2 movelist[50];
                  capture = Test_Capture(b,color);
                  if ( capture )
   n = Gen_Captures(b,movelist,color);
   else
   n = Gen_Moves(b,movelist,color);

    --depth;
    for ( i = 0;i < n;i++){
          domove2(b,&movelist[i],color );
          if ( depth ) Perft(b,color^CHANGECOLOR,depth);
          else  ++PerftNodes;
          undomove(b,&movelist[i],color );
                                   }
}


int EdProbe(int c[46],int color)
{
   if (!ED) return EdAccess::not_found;

   unsigned i;
   EdAccess::EdBoard1 b;
                  static const int Map_32_to_45[32] = {
       8,  7,  6,  5,
      13, 12, 11, 10,
      17, 16, 15, 14,
      22, 21, 20, 19,
      26, 25, 24, 23,
      31, 30, 29, 28,
      35, 34, 33, 32,
      40, 39, 38, 37
   };

   if (color == WHITE)
   {
      for (i = 0; i < 32; i++)
      {
         switch (c[Map_32_to_45[i]])
         {
            case FREE           : b.board[i] = EdAccess::empty; break;
            case WHITE | MAN : b.board[i] = EdAccess::white; break;
            case BLACK | MAN : b.board[i] = EdAccess::black; break;
            case WHITE | KING: b.board[i] = EdAccess::white | EdAccess::king; break;
            case BLACK | KING: b.board[i] = EdAccess::black | EdAccess::king; break;
         }
      }
   }
   else
   {
      // при ходе черных "переворачиваем" доску
      for (i = 0; i < 32; i++)
      {
         switch (c[Map_32_to_45[31 - i]])
         {
            case FREE           : b.board[i] = EdAccess::empty; break;
            case WHITE | MAN : b.board[i] = EdAccess::black; break;
            case BLACK | MAN : b.board[i] = EdAccess::white; break;
            case WHITE | KING: b.board[i] = EdAccess::black | EdAccess::king; break;
            case BLACK | KING: b.board[i] = EdAccess::white | EdAccess::king; break;
         }
      }
   }

   return ED->GetResult(&b, 0);
}


//*************************************************
//*                                               *
//*               Kallisto support                *
//*                                               *
//*************************************************

int AnalyseMode = 0;
int PlayNow = 0;

int Board[46];
int Color;
int TimeRemaining;
int IncTime;

void Wait(int &v)
{
   while(v) Sleep(10);
}

void SquareToStr(short sq, char *s)
{
   static const int Square64[] = {
       0,  0,  0,  0,  0,
       7,  5,  3,  1,  0,
      14, 12, 10,  8,
      23, 21, 19, 17,  0,
      30, 28, 26, 24,
      39, 37, 35, 33,  0,
      46, 44, 42, 40,
      55, 53, 51, 49,  0,
      62, 60, 58, 56
   };

   sq = Square64[sq];
   s[0] = sq % 8 + 'a';
   s[1] = 8 - sq / 8 + '0';
   s[2] = 0;
}

void MoveToStr(move2 m, char *s)
{
   SquareToStr(m.m[0] & 255, s);
   for (int i = 2; i < m.l; i++) {
      strcat(s, ":");
      SquareToStr(m.m[i] & 255, s + strlen(s));
   }
   if (m.l > 2) strcat(s, ":");
   SquareToStr(m.m[1] & 255, s + strlen(s));
}

// Сделать ход move
// Формат ходов: "a3b4" и "a3:b4:d6:e7". Такой формат позволяет устранить все неоднозначности при взятиях
void KALLISTOAPI EI_MakeMove(char *move)
{
   if (AnalyseMode){
      PlayNow = 1;
      Wait(AnalyseMode);
   EnterCriticalSection(&AnalyseSection);
   LeaveCriticalSection(&AnalyseSection);
                                              }
   move2 ml[MAXMOVES];
   init(Board);
   int n = Gen_Captures(Board, ml, Color);
   if (!n) n = Gen_Moves(Board, ml, Color);
   for (int i = 0; i < n; i++) {
      char s[128];
      MoveToStr(ml[i], s);
      if (!strcmp(s, move)) {
         domove2(Board, &ml[i],Color);
         if (Color == WHITE) Color = BLACK;
         else Color = WHITE;
         return;
      }
   }
   MessageBox(0, "KestoG: move not found", move, MB_OK);
}

// Начать вычисления. После вернуть лучший ход
// Эта функция может выполняться как угодно долго
// Но надо иметь в виду количество оставшегося времени (см. EI_SetTimeControl и EI_SetTime)
char *KALLISTOAPI EI_Think()
{
   //Create_HashFunction();
   char dummy[256];
   PlayNow = 0;
   play = &PlayNow;
   int time_limit = TimeRemaining / 20 + IncTime;
   if (time_limit * 2 > TimeRemaining) time_limit = TimeRemaining / 2;
                  // init(Board);
   compute(Board, Color, time_limit, dummy);
                  static char s[128];
   MoveToStr(bestrootmove, s);
   if (Color == WHITE) Color = BLACK;
   else Color = WHITE;
   return s;
}

// Здесь можно делать что угодно и как угодно долго
// Эта функция вызывается в момент когда противник думает над своим ходом
void KALLISTOAPI EI_Ponder()
{
   // здесь можно ничего и не делать :)
   return;
}

// Противник делает ход move
// Перед этим вызывалась функция Ponder
// Здесь сразу можно вернуть ход на основе вычиcлений сделанных в Ponder
// Можно подумать еще и только после этого вернуть ход
char *KALLISTOAPI EI_PonderHit(char *move)
{
   EI_MakeMove(move);
   return EI_Think();
}

// Инициализация движка
// si - см. выше описание PF_SearchInfo
// mem_lim - лимит памяти, которую может использовать движок
// здесь в основном имеется ввиду размер хэш-таблицы
void  KALLISTOAPI EI_Initialization(PF_SearchInfo si, int mem_lim)
{
   pfSearchInfo = si;
                  InitializeCriticalSection(&AnalyseSection);
                  size = (unsigned int)mem_lim;
//                MessageBox(0, "Initialization", "", MB_OK);
//   Create_HashFunction();
}

// Закончить вычисления и выйти из функций EI_Think, EI_Ponder, EI_PonderHit или EI_Analyse
void KALLISTOAPI EI_Stop()
{
   PlayNow = 1;
}

// Установить позицию pos на доске
// пример: начальная позиция bbbbbbbbbbbb........wwwwwwwwwwwww
// b - простая черная
// B - черная дамка
// w - простая белая
// W - белая дамка
// . - пустое поле
// поля перечисляются так: b8, d8, f8, h8, a7, c7, ..., a1, c1, e1, g1
// последний символ определяет очередность хода
// w - белые, b - черные
void KALLISTOAPI EI_SetupBoard(char *p)
{
   static const int Map[32] = {
       8,  7,  6,  5,
      13, 12, 11, 10,
      17, 16, 15, 14,
      22, 21, 20, 19,
      26, 25, 24, 23,
      31, 30, 29, 28,
      35, 34, 33, 32,
      40, 39, 38, 37
   };

   int i;
   for(i = 0; i < 46; i++) Board[i] = OCCUPIED;
   for (i = 0; i < 32; i++) {
      switch (p[i]) {
      case 'w': Board[Map[i]] = WHITE | MAN; break;
      case 'W': Board[Map[i]] = WHITE | KING; break;
      case 'b': Board[Map[i]] = BLACK | MAN; break;
      case 'B': Board[Map[i]] = BLACK | KING; break;
      case '.': Board[Map[i]] = FREE; break;
      }
   }
   if (p[32] == 'w') Color = WHITE;
   else Color = BLACK;
                  TTableInit(size);
                  searches_performed_in_game = 0;
                  Create_HashFunction();
}

void KALLISTOAPI EI_NewGame()
{
          EI_SetupBoard("bbbbbbbbbbbb........wwwwwwwwwwwww");
//                  TTableInit(size);
//                  searches_performed_in_game = 0;
//                  Create_HashFunction();
}

// Установить контроль времени
// time минут на партию
// inc секунд - бонус за каждый сделанный ход (часы Фишера)
void KALLISTOAPI EI_SetTimeControl(int time, int inc)
{
   TimeRemaining = time * 60 * 1000;
   IncTime = inc;
}

// Установить время в миллисекундах оставшееся на часах
// time - свое время
// otime - время противника
void KALLISTOAPI EI_SetTime(int time, int otime)
{
   TimeRemaining = time;
}

// Вернуть название движка
char *KALLISTOAPI EI_GetName()
{
   return "KestoG 1.4 Moscow";
}

// Вызывается перед выгрузкой движка
void KALLISTOAPI EI_OnExit()
{
}

// Анализировать текущую позицию
// Выход из режима анализа осуществляется при получении команд Stop или МакеMove
void KALLISTOAPI EI_Analyse()
{
   AnalyseMode = true;
   char dummy[256];
   PlayNow = 0;
   play = &PlayNow;
   EnterCriticalSection(&AnalyseSection);
                  // init(Board);
                  TTableInit(size);
                  searches_performed_in_game = 0;
                  Create_HashFunction();
   compute(Board, Color, 2000000000, dummy); // infinity thinking
   undomove(Board, &bestrootmove,Color);
   AnalyseMode = false;
                  LeaveCriticalSection(&AnalyseSection);
}


// функция интерфейса экспортируемая из dll
void KALLISTOAPI EI_EGDB(EdAccess *eda)
{

   ED = eda;
   if (ED)
   {
                    EdPieces = ED->Load("russian");
                    if (strstr(ED->GetBaseType(), "nocaptures"))
                    EdNocaptures = true;
   }
//                  assert(EdPieces != 0);
                  // if ( EdPieces == 0 ) exit(0);

}
