#pragma once
#include <cshaders.h>

class CVector;

class TestShader : public CShaderParent {
 public:
  TestShader();

  CVector lightpos;
  CVector targetdir;

  void Initialize(string);
  void Start();
  void End();

};
