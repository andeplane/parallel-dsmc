#include <TestShader.h>

void TestShader::Initialize(string s) {

  ((CShaderParent*)this)->Initialize(s);

}

void TestShader::Start() {
  Shader->begin();

  Shader->sendUniform3f((char*)"lightpos",lightpos.x, lightpos.y, lightpos.z);

  // Shader->sendUniform3f((char*)"targetdir",targetdir.x, targetdir.y, targetdir.z);

}

void TestShader::End() {
  Shader->disable_multitextures();
  Shader->end();
}

TestShader::TestShader() {
  vert = string(
    "uniform vec3 lightpos; \n"
    // "uniform vec3 targetdir; \n"
    "varying vec3 normal; \n"
    "varying vec4 myPos; \n"
    "void main(void) \n"
    "{ \n"
    "	normal = gl_Normal; \n"
    "   myPos = gl_Vertex;"
    "   gl_TexCoord[0] = gl_MultiTexCoord0;\n"
    "   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
    "}\n");

  frag = string(
    "uniform vec3 lightpos; \n"
    "varying vec3 normal; \n"
    //"uniform vec3 targetdir; \n"
    "varying vec4 myPos; \n"
    "void main(void)\n"
    "{\n "
    "  float r = 151.0/255.0;"
    "  float g = 171.0/255.0;"
    "  float b = 210.0/255.0;"
    "  vec3 ambient = 0.6*vec3(r,g,b);"
    "  vec3 diffuse = 1.5*ambient;"
    "  float distance = length(lightpos - myPos.xyz);"
    "  float one_over_color_cutoff = 0.1;"
    "  float color_factor = max( (1.0-distance*one_over_color_cutoff),0.0);"
    "  vec3 toLight = normalize(lightpos - myPos.xyz);"
    "  float angle  = max(dot(normal, toLight), 0.0);"
    "  gl_FragColor = vec4(color_factor*(ambient + angle*diffuse), 1.0)\n;"
    // "  gl_FragColor = vec4(angle,0,0,1);"
    "  gl_FragColor.w = 1.0;"
    "}\n");

//  frag = string(
//      "uniform vec3 lightpos; \n"
//      "varying vec3 normal; \n"
//      "uniform vec3 targetdir; \n"
//      "varying vec3 myPos; \n"
//      "void main(void)\n"
//      "{\n "
//      "  vec4 val = vec4(0.7,0.5,0.3,1.0);"
////      "  if(myPos.z > 1.1) {"
////      "      val = vec4(10.0*(myPos.z-1.0),10.0*(myPos.z-1.0),10.0*(myPos.z-1.0),1.0);"
////      "  }                  "
//      "  float light = clamp(dot(normalize(lightpos), normal), 0.4, 1.0);"
//      "  float shininess = 10.0;"
//      // "  float specular = pow(clamp(dot(reflect(-normalize(lightpos), normal), targetdir), 0.0, 1.0), shininess);"
//      "  float specular = 0.0;\n"
//      "  gl_FragColor = val*light + vec4(1,1,1,1)*specular; \n"
//      "  gl_FragColor.w = 1.0;"
//      "}\n");




}
