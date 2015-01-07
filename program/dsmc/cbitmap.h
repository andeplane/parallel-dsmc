// pixmap.h
#pragma once

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <string>

using namespace std;


typedef struct {
   unsigned short int type;                 /* Magic identifier            */
   unsigned int filesize;                       /* File size in chars          */
   unsigned short int reserved2;
   unsigned int offset;                     /* Offset to image data, chars */

   unsigned int size;               /* Header size in chars      */
   int width,height;                /* Width and height of image */
   unsigned short int planes;       /* Number of colour planes   */
   unsigned short int bits;         /* Bits per pixel            */
   unsigned int compression;        /* Compression type          */
   unsigned int imagesize;          /* Image size in chars       */
   int xresolution,yresolution;     /* Pixels per meter          */
   unsigned int ncolours;           /* Number of colours         */
   unsigned int importantcolours;   /* Important colours         */

} CBMPHeader;


//typedef char char;
//typedef unsigned char unsigned char;

class CBitMap
{
 public:
  unsigned width, height;
  unsigned char *data;
  
  enum {RED = 0, GREEN, BLUE};
  

  CBitMap();
  
  CBitMap(char *fname);
         
  ~CBitMap();
  
  void LoadBMP(const char *fname);
  void SaveBMP(string filename);

  void Create(int w, int h);
  
  unsigned char pixel_elem(int x, int y, int elem);
  unsigned char *pixel_pos(int x, int y);	
};



