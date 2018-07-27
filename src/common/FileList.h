/* --- --- ---
 * Copyright (C) 2008--2010 Idiap Research Institute (.....@idiap.ch)
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
// FileList.h: interface for the CFileList class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FILE_LIST_H_)
#define _FILE_LIST_H_

#include <ctime>						// clock
#include <cstdlib>						// C standard library
#include <cstdio>						// C I/O (for sscanf)
#include <cstring>						// string manipulation
#include <fstream>						// file I/O
#include <cmath>						// math includes
#include <iostream>						// I/O streams
#include <string>

using namespace std;


class CFileList
{
public:
	/* import the file name list */
	bool Read(char *file_list_name);

	/* get list length */
	int GetListLength();

	/* get the line number of the given file name with the index file_loc */
	int GetFileNameLocation(int file_loc);

	/* get the file name given the index file_loc */
	char* GetFileName(int file_loc);

	/* get all the file names */
	char** GetFileNames();

	/* get the file list name */
	char* GetFileListName();

	/* export file names into a file name list */
	bool Write(char *file_list_name, char** file_names, int length);

	/* export file names into a file name list */
	bool Write(char *file_list_name=NULL);

	/* initialize data */
	void Init(int list_length);

	/* set the length of all file names in the list */
	void SetListLength(int list_length);

	/* set file name on the particular location */
	void SetFileName(int file_loc, char* file_name);

	/* set all file names */
	void SetFileNames(char** file_names, int length);

	CFileList(char *file_list_name);

	CFileList();
	virtual ~CFileList();

private:
	/* remove the left and right blank charaters in one line */
	void FileLineFilter(char* line);

	/* clean data */
	void CleanData();

	char* m_pFileListName;
	char** m_ppFileNames;
	int* m_pFileNameLineNOs;
	int m_nListLength;
	int m_nMaxListLength;
};

#endif // !defined(_FILE_LIST_H_)
