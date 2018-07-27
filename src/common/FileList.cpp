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
// FileList.cpp: implementation of the CFileList class.
//
//////////////////////////////////////////////////////////////////////

#include "FileList.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CFileList::CFileList()
{
	m_pFileListName = NULL;
	m_ppFileNames = NULL;
	m_pFileNameLineNOs = NULL;
	m_nListLength = 0;
	m_nMaxListLength = 0;
}

CFileList::CFileList(char *file_list_name)
{
	m_pFileListName = NULL;
	m_ppFileNames = NULL;
	m_pFileNameLineNOs = NULL;
	m_nListLength = 0;
	m_nMaxListLength = 0;

	Read(file_list_name);
}

CFileList::~CFileList()
{
	CleanData();
}

void CFileList::Init(int list_length)
{
	if ( list_length > m_nMaxListLength ) {
		char** tmpFileNames = new char*[list_length];

		if ( m_pFileNameLineNOs )
			delete [] m_pFileNameLineNOs;
		m_pFileNameLineNOs = new int[list_length];

		for ( int i = 0 ; i < list_length ; i++ ) {
			m_pFileNameLineNOs[i] = i+1;
			tmpFileNames[i] = new char[1024];
		}

		for ( int i = 0 ; i < list_length ; i++ ) {
			tmpFileNames[i] = new char[1024];
			if ( i < m_nMaxListLength )
				strcpy(tmpFileNames[i], m_ppFileNames[i]);
		}

		for ( int i = 0 ; i < m_nMaxListLength ; i++ )
			delete [] m_ppFileNames[i];
		delete [] m_ppFileNames;

		m_ppFileNames = tmpFileNames;

		m_nMaxListLength = list_length;
	}
}

void CFileList::SetListLength(int list_length)
{
	if ( list_length > m_nMaxListLength ) {
		Init(list_length);
	}
	m_nListLength = list_length;
}

void CFileList::SetFileName(int file_loc, char* file_name)
{
	if ( file_loc+1 > m_nListLength )
		m_nListLength = file_loc+1;

	strcpy(m_ppFileNames[file_loc], file_name);
}

void CFileList::SetFileNames(char** file_names, int length)
{
	CleanData();
	m_nListLength = length;
	Init(length);
	for ( int i = 0 ; i < length ; i++ )
		strcpy(m_ppFileNames[i], file_names[i]);
}

int CFileList::GetListLength()
{
	return m_nListLength;
}

char* CFileList::GetFileName(int file_loc)
{
	if ( file_loc < 0 || file_loc >= m_nMaxListLength ) {
		//printf("Please provide correct file location in range [0, %d] \n", m_nListLength);
		return NULL;
	}

	return m_ppFileNames[file_loc];
}

int CFileList::GetFileNameLocation(int file_loc)
{
	if ( file_loc < 0 || file_loc >= m_nMaxListLength ) {
		//printf("Please provide correct file location in range [0, %d] \n", m_nListLength);
		return 0;
	}

	return m_pFileNameLineNOs[file_loc];
}

char** CFileList::GetFileNames()
{
	return m_ppFileNames;
}

char* CFileList::GetFileListName()
{
    return m_pFileListName;
}

void CFileList::CleanData()
{
	if ( m_pFileListName )
		delete [] m_pFileListName;
	if ( m_ppFileNames ) {
		for ( int i = 0 ; i < m_nMaxListLength ; i++ )
            delete [] m_ppFileNames[i];
		delete [] m_ppFileNames;
	}
	if ( m_pFileNameLineNOs )
		delete [] m_pFileNameLineNOs;

	m_pFileListName = NULL;
	m_ppFileNames = NULL;
	m_pFileNameLineNOs = NULL;
	m_nListLength = 0;
	m_nMaxListLength = 0;
}


bool CFileList::Write(char *file_list_name, char** file_names, int length)
{
    if ( !file_list_name || strlen(file_list_name) == 0 )
        return false;

	ofstream fout;
	fout.open(file_list_name, ios::out);

	if (fout.fail()) {
		printf("Error opening file name list %s.\n", file_list_name);
		fout.close();
		exit(0);
	}
	else {
		for ( int i = 0 ; i < length ; i++ ) {
			fout << file_names[i] << endl;
		}

		fout.close();
	}

	return true;
}

bool CFileList::Write(char *file_list_name)
{
    if ( file_list_name && strlen(file_list_name) > 0 )
        return Write(file_list_name, m_ppFileNames, m_nListLength);
    else
        return Write(m_pFileListName, m_ppFileNames, m_nListLength);
    return false;
}

bool CFileList::Read(char *file_list_name)
{
	/* clean the data */
	CleanData();

	if ( !file_list_name || strlen(file_list_name) == 0 )
		return false;

	ifstream fin;
	fin.open(file_list_name, ios::in);

    /* get the file list name */
    m_pFileListName = new char[strlen(file_list_name)+1];
    sprintf(m_pFileListName, "%s", file_list_name);

	if (fin.fail()) {
		//printf("Error opening file name list %s.\n", file_list_name);
		fin.close();
		return false;
		//exit(0);
	}
	else {
		/* get the valid file name lenght */
		char temp[1024];
		m_nListLength = 0;
		/* alopez: replaced old code with bug, last line was not read*/
/*		while (1)
		{
			fin.getline(temp, 1024);
			FileLineFilter(temp);
			if ( fin.eof() ) //alopez: eof checked just after last line is read, a line is typically missed
				break;
			else if ( strlen(temp) > 0 && temp[0] != '#' )
				m_nListLength++;
		}*/
//		std::string msg;
		while(!fin.eof())
		{
			fin.getline(temp, 1024);
			FileLineFilter(temp);
			if ( strlen(temp) > 0 && temp[0] != '#' )
				m_nListLength++;
		}
		//if lines end with \n last line  is empty

		std::cout << "File list length " << m_nListLength << std::endl;

		/* go back to file start point */
		fin.clear();
		fin.seekg(0, ios::beg);

		/* allocate memories */
		m_ppFileNames = new char*[m_nListLength];
		m_pFileNameLineNOs = new int[m_nListLength];

		/* get the valid file names */
		m_nListLength = 0;
		int line_no = 0;
		while (!fin.eof())
		{
			fin.getline(temp, 1024);
			line_no++;
			FileLineFilter(temp);
/*			if ( fin.eof() )
				break;*/
			if ( strlen(temp) > 0 && temp[0] != '#' )
			{
				m_ppFileNames[m_nListLength] = new char[strlen(temp)+1];
				sprintf(m_ppFileNames[m_nListLength], "%s", temp);
				m_pFileNameLineNOs[m_nListLength] = line_no;
				m_nListLength++;
			}
		}

		fin.close();
	}

	m_nMaxListLength = m_nListLength;

	return true;
}

void CFileList::FileLineFilter(char* line)
{
	// Remove leading and trailing whitespace
	static const char whitespace[] = " \n\t\v\r\f";
	string s(line,1024);
	s.erase( 0, s.find_first_not_of(whitespace) );
	s.erase( s.find_last_not_of(whitespace) + 1U );

	sprintf(line, "%s", s.c_str());
	return;

	/////////////////////////////////////

	int i;

	/* remove the left blank charaters */
	int left_blank_num = 0;
	int line_len = strlen(line);

	for ( i = 0 ; i < line_len ; i++ ) {
		if ( line[i] != ' ' && line[i] != '\t' )
			break;
		left_blank_num++;
	}

	if ( left_blank_num ) {
		printf("remove left blanks\n");
		for ( i = 0 ; i < line_len - left_blank_num ; i++ )
			line[i] = line[i+left_blank_num];
		line[i] = '\0';
	}

	/* remove the right blank charaters */
	int right_blank_num = 0;
	line_len = strlen(line);
	for ( i = line_len-1 ; i >= 0 ; i++ ) {
		if ( line[i] != '\t' && line[i] != ' ' )
			break;
		right_blank_num++;
	}
	if ( right_blank_num )
		line[strlen(line)+1-right_blank_num] = '\0';
}
