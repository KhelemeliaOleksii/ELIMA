#ifndef READNLINE_H
#define READNLINE_H

#include<conio.h>
#include<io.h>
#include <stdlib.h> // string
#include<stdio.h>

#define MSGLENGTH 500 
// ������� readNline ��������� 
// n-�� ������ �� ��������� ����� "in"
// ���������� ������ ���� double ������� count, 
// ������������ � �������� ������
// �. � ��������� ��������� ����
//		������������ �������� 1
//		����������� ������ ������ msg
//	�����:
//		������������ �������� 0
int readNline (FILE **in, unsigned int number, float *par, int count, char *msg[]);
#endif