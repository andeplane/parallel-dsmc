#include <TestShader.h>

void TestShader::Initialize(string s) {

  ((CShaderParent*)this)->Initialize(s);

}

void TestShader::Start() {
  Shader->begin();

  Shader->sendUniform3f((char*)"lightpos",lightpos.x, lightpos.y, lightpos.z);

  Shader->sendUniform3f((char*)"targetdir",targetdir.x, targetdir.y, targetdir.z);

}

void TestShader::End() {
  Shader->disable_multitextures();
  Shader->end();
}

TestShader::TestShader() {
  vert = string(
    "uniform vec3 lightpos; \n"
    "uniform vec3 targetdir; \n"
    "varying vec3 normal; \n"
    "varying vec3 myPos; \n"
    "void main(void) \n"
    "{ \n"
    "	  normal = gl_Normal; \n"
    "	  myPos = gl_Vertex.xyz; \n"
    "   gl_TexCoord[0] = gl_MultiTexCoord0;\n"
    "   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
    "}\n");

  frag = string(
    "uniform vec3 lightpos; \n"
    "varying vec3 normal; \n"
    "uniform vec3 targetdir; \n"
    "varying vec3 myPos; \n"
    "void main(void)\n"
    "{\n "
    // "  vec4 val = vec4(0.5+myPos.x*0.1,0.5+myPos.y*0.1,0.5,1.0);"
    "  vec4 val = vec4(89.0/255.0, 193.0/255.0, 235.0/255.0, 1.0);"
    // "  if(myPos.z > 1.1) {"
    // "      val = vec4(10.0*(myPos.z-1.0),10.0*(myPos.z-1.0),10.0*(myPos.z-1.0),1.0);"
    // "  }                  "
    "  float light = clamp(pow(dot(normalize(lightpos), normal),0.3), 0.3, 1.0);"
    "  float shininess = 10.0;"
    "  float specular = pow(clamp(dot(reflect(-normalize(lightpos), normal), targetdir), 0.0, 1.0), shininess);"
    "  gl_FragColor = val*light + vec4(0.2,0.2,0.2,1)*specular; \n"
    "  gl_FragColor.w = 1.0;"
    "}\n");




}
