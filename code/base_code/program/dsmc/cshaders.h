#pragma once

#include <copengl.h>
#include <oglshader.h>

class CShaderParent {
 protected:
  string frag;
  string vert;
  CShaderObject* Shader;
  int tangent_loc;
 public:
  GLuint ScreenTexture;
  void Initialize(string s);
  CShaderParent();
  virtual void Start() = 0;
  virtual void End() = 0;

  static COpenGL* COpenGLPointer;
  void ScreenToTexture(GLuint texture);

};

class CBlurAndToonShader: public CShaderParent{
 public:
  const static int BLUR = 0;
  const static int TOON1 = 1;
  const static int TOON2 = 2;
  const static int TOON3 = 3;
  const static int TOON4 = 4; // only pixelizing
  float Offsets[25 * 2];
  float Kernel[25*4];

  int Type;
  float Value, Width;
  CVector Haze;
  string frag_start, frag_end;
  string frag_type1, frag_type2, frag_type3;

  CBlurAndToonShader();
  void Initialize(string name, int type, double value, CVector haze, double width);
  void Start();
  void End();
  void Render();
  void ScreenToTexture();
  void NewKernel(double);
};

class CMathShader : public CShaderParent {
 public:
  CMathShader();
  const static int ADD = 0;
  const static int SUB = 1;
  const static int MUL = 2;
  const static int DIV = 3;
  const static int MYBODY = 4;

  int Type;
  double Weight;
  double Distort;

  string mybody;


  void Initialize(string name, int type, double weight);
  void Start();
  void End();
  void Render();
  GLuint ScreenTexture2;

};

class CBumpShader : public CShaderParent {
 public:
  GLuint ColorMap, NormalMap;
  float InvertRadius;
  float NormalMapScale;

  CBumpShader();
  void Tangent(CVector t);
  void Start();
  void End();
};

class CBillboardLeafShader: public CShaderParent {
 public:
  CBillboardLeafShader();
  GLuint leaftexture;
  void Start();
  void End();
};

class CThinLensShader: public CShaderParent {
 public:
  CThinLensShader();
  GLuint texture;
  CVector center;
  float lensStrength;
  float Dd;
  float type;
  void Start();
  void End();
  void Render();
};

class C3DLensShader: public CShaderParent {
 public:
  C3DLensShader();
  GLuint textureX;
  GLuint textureY;
  GLuint textureZ;
  //CVector center;
  float lensStrength, scale;
  CVector delta, min, camera;
  //float Dd;
  float type;
  float width;
  void Start();
  void End();
  void Render();
};

