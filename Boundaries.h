// Boundaries.h определяет гарницы интегирования 

#include"InputData.h"

#ifndef BOUNDARIES_H
#define BOUNDARIES_H
#define A 1.E5
#define B 1.E5
// for ELIA
#define LOWX -A
#define UPX A
#define LOWY -A
#define UPY A
#define LOWZ -B
#define UPZ B

//#define LOWZ 0.
//#define UPZ 1.14711678
//#define UPZ 0.1*TAU/DELTA/DELTA
// по Y
//#define LOWX 0.
//#define UPX 1.
//#define LOWY 0.
//#define UPY 1.
//#define LOWZ 0.
//#define UPZ 1.
// по Z
// если интегрирование по некоторой переменной отсутствует, то задаем границы типа
//#define LOW 0
//#define UP 1
//, тогда проведенное интегрирование по этой переменной дает 1, что не влияет на результат
#endif