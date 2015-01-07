// CBitMap.cpp
#include <CBitMap.h>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <CUtil.h>
#include <CVector.h>
#include <memory.h>



void CBitMap::Create(int w, int h) {
    if (data!=0)
      delete[] data;
    data = new unsigned char[3*w*h];
    for (int i=0;i<3*w*h;i++)
      data[i] = 0;
    width = w;
    height = h;
  }


void CBitMap::SaveBMP(string filename) {
  
  CBMPHeader b;
  
  b.type = 19778;
  b.offset = 54;
  
  b.filesize = width*height*3 + 54; // file size
  b.reserved2=0;
  b.size = 40;
  b.width = width; 
  b.height = height;
  b.planes = 1;
  b.bits=24;
  b.compression = 0; // 0-compressed
  b.imagesize=0; // 0 for compressed
  
  int winWidth = width;
  int winHeight = height;
  
  
  ofstream f(filename.c_str(), ios::out | ios::binary);
  
  int n = 0;
  
  f.write ((char *)&b.type, 2 );
  f.write ((char *)&b.filesize, 4 );
  f.write ((char *)&b.reserved2, 2 );
  f.write ((char *)&b.reserved2, 2 );
  f.write ((char *)&b.offset, 4 );
  f.write ((char *)&b.size, 4 ); // 14
  f.write ((char *)&b.width, 4 );  // 18
  f.write ((char *)&b.height, 4 );  // 22
  f.write ((char *)&b.planes, 2 );  // 24
  f.write ((char *)&b.bits, 2 );  // 26
  f.write ((char *)&b.compression, 4 ); //30
  f.write ((char *)&b.imagesize, 4 ); // 34
  f.write ((char *)&n, 4 );  //38
  f.write ((char *)&n, 4 );  // 42
  f.write ((char *)&n, 4 );  //46
  f.write ((char *)&n, 4 ); // 50
  
  
  f.write ((char *)data , winWidth*winHeight*3);
  f.write ((char *)&n, 2); // 50
  
  f.close ();
  
  return;  
  
  
}



CBitMap::CBitMap(char *fname)
	: width(0), height(0), data(0)
{
	this->LoadBMP(fname);
}

CBitMap::CBitMap()
	: width(0), height(0), data(0) {}

CBitMap::~CBitMap()
{
	if( data ) {
    	delete[] data;
    }
 }

void CBitMap::LoadBMP(const char *fname)
{
	  using namespace std;

	  unsigned short planes;	// number of planes in image (must be 1) 
	  unsigned short bpp;			// number of bits per pixel (must be 24)
	  
     CUtil::verify_file(fname);
	  ifstream fin(fname, ios::in | ios::binary);
	  
	  fin.seekg(18, ios::cur);
	  
	  fin.read((char *)&width, sizeof(unsigned));
	  fin.read((char *)&height, sizeof(unsigned));
	  //cout << "width: " << width << " height: " << height << '\n';
	  
	  fin.read((char *)&planes, sizeof(unsigned short));
	  if( planes != 1 )
	    {
	      //cout << "Planes from " << fname << " is not 1: " << planes << "\n";
	      //exit(1);
          throw string("Planes in texture is incorrect: " +string(fname));
	    }
	  
	  fin.read((char *)&bpp, sizeof(unsigned short));
	  if( bpp != 24 )
	    {
          throw string("bpp is not 24 in texture: " +string(fname));
	      //cout << "Bpp from " << fname << " is not 24: " << bpp << "\n";
	      //exit(1);
	    }
	  
	  fin.seekg(24, ios::cur);
	  
	  unsigned size(width * height * 3);				// size of the image in chars (3 is to RGB component).
	  data = new unsigned char[size];
	  fin.read((char *)data, size);
	  
	  unsigned char tmp;					// temporary color storage for bgr-rgb conversion.
	  for( unsigned int i(0); i < size; i += 3 )
	    {
	      tmp = data[i];
	      data[i] = data[i+2];
	      data[i+2] = tmp;
	    }
}

unsigned char CBitMap::pixel_elem(int x, int y, int elem)
{
	int pos = (y*width+x) * 3 + elem;
	return data[pos];
}

unsigned char *CBitMap::pixel_pos(int x, int y)
{
	int pos = (y * width + x) * 3;
	return &data[pos];
}