#ifndef ED_ACCESS_H
#define ED_ACCESS_H

// ������������ ����� ��� ������� � ����������� �����
struct EdAccess
{
	// ��������� ��� ������������� �������� ������� 

	struct EdBoard1
	{
		// ��� ���� ���� �� ������� b8, d8, � �.�. �� g1
		unsigned char board[32];
	};

	struct EdBoard2
	{
		unsigned char *wman;
		unsigned wman_cnt;

		unsigned char *wkings;
		unsigned wkings_cnt;

		unsigned char *bman;
		unsigned bman_cnt;

		unsigned char *bkings;
		unsigned bkings_cnt;
	};


	// ������������ ��������

    enum
	{
		draw      =      0,
		win       =  10000,
		lose      = -10000,
		not_found =  32000
	};


	// ���� ��� ����������� �����

	enum
	{
		white = 1,
		black = 2,
		empty = 4,
		king  = 8
	};


	// ����� ��� ������� � ����

	enum
	{
		in_mem = 1
	};


	// ��������� ����
	// ���� ����� ���� ���:
	// 		russian
	// 		russianlosers
	// 		brazil
	// 		brazillosers
	// 		pool
	// 		poollosers
	// 		checkers
	// 		checkerslosers

	virtual unsigned __stdcall Load(char *game_type) = 0;


	// �������� ��� ����

	virtual char * __stdcall GetBaseType() = 0;


	// ������ ������� (������ ��� �����)

	virtual int __stdcall GetResult(EdBoard1 *board, unsigned flags) = 0;
	virtual int __stdcall GetResult(EdBoard2 *board, unsigned flags) = 0;


	// �������� ��������� �� ������� �� ���������

	virtual unsigned __stdcall GetTable(unsigned wm, unsigned wk, unsigned bm, unsigned bk) = 0;


	// �������� ��������� �� ������� �� ��������� � �� �������� ����������� �����

	virtual unsigned __stdcall GetTable(unsigned wm, unsigned wk, unsigned bm, unsigned bk, unsigned rank) = 0;


	// �������� ������������� ������� ������� � ������

	virtual unsigned __stdcall IsTableInMemory(unsigned table) = 0;

	
	// �������� ������ � �������

	virtual unsigned __int64 __stdcall GetIndex(EdBoard1 *board) = 0;
	virtual unsigned __int64 __stdcall GetIndex(EdBoard2 *board) = 0;

	
	// �������� ������ �� ��������� �� ������� � �������

	virtual int __stdcall GetResult(unsigned table, unsigned __int64 index, unsigned flags) = 0;
};

#endif