#include <cshaders.h>
#include <camera.h>

COpenGL* CShaderParent::COpenGLPointer;

CShaderParent::CShaderParent() {
  Shader = 0;
}
void CShaderParent::Initialize(string s) {
  if (COpenGLPointer==0)
    throw string("Error in CShaderParent::Initialize - COpenGLPointer not set");
  if (!Shader) {
    COpenGLPointer->shadercontainer.addshaderfrommemory(s,vert, frag);
    Shader = COpenGLPointer->shadercontainer.getshader(s);

    glBindAttribLocation(Shader->ShaderObject, 6, "vTangent");
    glLinkProgram(Shader->ShaderObject);
    tangent_loc = glGetAttribLocation(Shader->ShaderObject, "vTangent");

  }

  glEnable(GL_TEXTURE);
  glGenTextures(1, &ScreenTexture);
  glBindTexture(GL_TEXTURE_2D, ScreenTexture);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);


}

void CShaderParent::ScreenToTexture(GLuint texture) {
  COpenGLPointer->buffer2texture(texture, COpenGLPointer->window_width, COpenGLPointer->window_height, COpenGL::NOMIPMAP);
}


CBlurAndToonShader::CBlurAndToonShader() {
  Width = -1;
  vert = string("uniform vec4 haze_color;\n"
        "uniform vec4 camera_input;\n "
        "uniform float width; \n"
        "varying vec4 camera;\n"
        "void main(void) {\n"
        "   gl_TexCoord[0] = gl_MultiTexCoord0;\n"
        "   gl_Position = gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
        "   camera = camera_input;\n"
        "}\n"
        );
  frag_start = string(
        "uniform sampler2D tex1;\n"
        "varying vec4 camera;\n"
        "uniform float width; \n"
        "uniform vec4 haze_color;\n"
        "uniform float blurscale;\n"
        "uniform float type;\n"
        "uniform vec2 o[25];\n"
        "uniform vec4 k[25];\n"
        "vec4 colorclamp(vec4 vv) \n"
        "{ \n"
        "     vec4 v1 = vv;                  \n   "
        "     for (int i=0;i<3;i++) { \n"
        "      if (v1[i]<0.1) v1[i]=0.1; \n"
        "      if (v1[i]>=0.2 && v1[i]<0.4) v1[i]=0.3; \n"
        "      if (v1[i]>=0.4 && v1[i]<0.6) v1[i]=0.5; \n"
        "      if (v1[i]>=0.6 && v1[i]<0.8) v1[i]=0.7; \n"
        "      if (v1[i]>=0.8) v1[i]=1.0;\n "
        "      }  	 \n"
        " //     v1[1] = v1[0];\n "
        " //     v1[2] = v1[0];\n "
        "     return v1;   \n "
        "    }\n"
        "vec4 kernelblur() {\n"
        "  vec4 total = vec4(0.0);\n"
        "  vec2 tx = gl_TexCoord[0].st;\n"
        "  if (width!=-1.0) { \n"
        "    tx[0]=float(int(tx[0]*width))/width; \n"
        "    tx[1]=float(int(tx[1]*width))/width; \n"
        //"  return texture2D(tex1, tx);\n"
        "  }  \n"
        "  if (type==4.0) return texture2D(tex1, tx);\n"
        "  for (int i=0; i<25; i++) {\n"
        "    vec4 tmp = texture2D(tex1, tx + o[i]);\n"
        "    total += tmp * k[i];\n"
        "  }\n ");
  frag_type1 = string("if (type!=0.0) {total.y = total.x; total.z=total.x;}\n"
                 "  vec4 v = texture2D(tex1, tx);\n"
                 "    total = (blurscale * (1.0-total)*v + v*(1.0-blurscale));\n");
  frag_type2 = string("if (type!=0.0) {total.y = total.x; total.z=total.x;}\n"
                 "  vec4 v = texture2D(tex1, tx);\n"
                 "  if (length(total)>blurscale)\n"
                 "      total=vec4(0,0,0,v.w);//total - vec4(length(total));\n"
                 "    else total = v;//colorclamp(v);\n");
  frag_type3 = string("if (type!=0.0) {total.y = total.x; total.z=total.x;}\n"
              "  vec4 v = texture2D(tex1, tx);\n"
                 "    total = vec4(1,1,1,1)-total;                \n "
                 "    total = (blurscale * (1.0-total)*v + v*(1.0-blurscale));\n");
  frag_end = string("  return total*(haze_color);\n"
               "}\n"
               "void main(void)\n"
               "{\n "
               "  gl_FragColor = kernelblur();\n"
               "}\n");


}

void CBlurAndToonShader::Render() {
  if (!Shader)
    return;
  COpenGLPointer->SetOrthographicProjection();

  COpenGLPointer->push();
  glLoadIdentity();

  glDisable(GL_DEPTH_TEST);
  glEnable(GL_TEXTURE_2D);
  glDisable(GL_FOG);
  glDisable(GL_CULL_FACE);

  Start();

  CVector c(1,1,1);
  glBegin(GL_QUADS);
  glTexCoord2d(0,1);
  glVertex3f(0,0,0);
  glTexCoord2d(1,1);
  glVertex3f(COpenGLPointer->window_width,0,0);
  glTexCoord2d(1,0);
  glVertex3f(COpenGLPointer->window_width,COpenGLPointer->window_height,0);
  glTexCoord2d(0,0);
  glVertex3f(0,COpenGLPointer->window_height,0);
  glEnd();

  End();

  COpenGLPointer->pop();

  COpenGLPointer->ResetPerspectiveProjection();

  glDisable(GL_TEXTURE_2D);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

}

void CBlurAndToonShader::Start() {

  Shader->begin();
  Shader->setup_shader_texture(ScreenTexture, 0, "tex1");
  Shader->sendUniform4f((char*)"haze_color", Haze.x,Haze.y,Haze.z, 1.0);
  Shader->sendUniform4f((char*)"camera_input", COpenGLPointer->camera->position.x,COpenGLPointer->camera->position.y,COpenGLPointer->camera->position.z, 1.0);
  if (Type!=BLUR) {
    if (Type!=TOON4) Shader->sendUniform1f((char*)"blurscale",Value);
    Shader->sendUniform1f((char*)"type",Type);
  }
  Shader->sendUniform1f((char*)"width",Width);
  Shader->sendUniform4ivArray((char*)"k", 25, Kernel);
  Shader->sendUniform2fv((char*)"o", 25, Offsets);



}

void CBlurAndToonShader::End() {

  Shader->disable_multitextures();
  Shader->end();

}



void CBlurAndToonShader::NewKernel(double width) {
  for(int x=0; x<5; x++) {
    for(int y=0; y<5; y++) {
      Offsets[2*(x + y*5)    ] = ((float)x - 2.0) / width;
      Offsets[2*(x + y*5) + 1] = ((float)y - 2.0) /width;
    }
  }
  float k[25];
  if (Type==BLUR) {
    k[ 0]=01; k[ 1]=04; k[ 2]=07; k[ 3]=04; k[ 4]=01;
    k[ 5]=04; k[ 6]=20; k[ 7]=33; k[ 8]=20; k[ 9]=04;
    k[10]=07; k[11]=33; k[12]=55; k[13]=33; k[14]=07;
    k[15]=04; k[16]=20; k[17]=33; k[18]=20; k[19]=04;
    k[20]=01; k[21]=04; k[22]=07; k[23]=04; k[24]=01;
  }
  else
  {
    k[ 0]=00; k[ 1]=00; k[ 2]=-1; k[ 3]=00; k[ 4]=00;
    k[ 5]=00; k[ 6]=-1; k[ 7]=-2; k[ 8]=-1; k[ 9]=00;
    k[10]=-1; k[11]=-2; k[12]=16; k[13]=-2; k[14]=-1;
    k[15]=00; k[16]=-1; k[17]=-2; k[18]=-1; k[19]=00;
    k[20]=00; k[21]=00; k[22]=-1; k[23]=00; k[24]=00;
  }
  double sum =1.0;

  if (Type==BLUR)
    for (int i=0;i<25;i++) sum+=k[i];
  for (int i=0;i<25;i++) {
    Kernel[4*i + 0] = k[i]/sum;
    Kernel[4*i + 1] = k[i]/sum;
    Kernel[4*i + 2] = k[i]/sum;
    Kernel[4*i + 3] = k[i]/sum;
  }


}


void CBlurAndToonShader::Initialize(string name, int type, double value, CVector haze, double width) {
  if (type==BLUR) frag = frag_start + frag_end;
  if (type==TOON1) frag = frag_start +frag_type1 + frag_end;
  if (type==TOON2) frag = frag_start +frag_type2 + frag_end;
  if (type==TOON3) frag = frag_start +frag_type3 + frag_end;
  if (type==TOON4) frag = frag_start  + frag_end;

  ((CShaderParent*)this)->Initialize(name);


  Haze = haze;
  Value = value;
  Type = type;

  // setup kernels & offsets
  NewKernel(width);
}

void CMathShader::Initialize(string s, int type, double weight) {


  if (type==MYBODY)
    frag = mybody;

  ((CShaderParent*)this)->Initialize(s);
  Type = type;
  Weight = weight;
  Distort = 0;

  glEnable(GL_TEXTURE);
  glGenTextures(1, &ScreenTexture2);
  glBindTexture(GL_TEXTURE_2D, ScreenTexture2);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);


}



void CBlurAndToonShader::ScreenToTexture() {
  COpenGLPointer->buffer2texture(ScreenTexture, COpenGLPointer->window_width, COpenGLPointer->window_height, COpenGL::NOMIPMAP);
}

CBumpShader::CBumpShader() {
  vert = string("varying vec3 lightVec; \n "
        "varying vec3 eyeVec; \n"
        "varying float fogFactor; \n"
        "varying vec2 texCoord; \n"
        "varying vec2 texCoord2; \n"
        "varying mat3 rotmat;"
        "uniform float NormalScale; \n"
        "attribute vec3 vTangent; \n"
        "void main(void) \n"
        "{ \n"
        " gl_Position = ftransform();\n"
        " texCoord = gl_MultiTexCoord0.xy;\n"
        " texCoord2 = gl_MultiTexCoord0.xy*NormalScale;\n"
        " \n"
        " vec3 n = normalize(gl_NormalMatrix * gl_Normal); \n"
        " vec3 t = normalize(gl_NormalMatrix * vTangent); \n"
        " vec3 b = normalize(cross(n, t)); \n"
        "  \n "
        " vec3 vVertex = vec3(gl_ModelViewMatrix * gl_Vertex); \n"
        " vec3 tmpVec = gl_LightSource[0].position.xyz - vVertex; \n"
        " \n "
        " lightVec.x = dot(tmpVec, t); \n"
        " lightVec.y = dot(tmpVec, b); \n"
        " lightVec.z = dot(tmpVec, n); \n"
        " \n"
        " tmpVec = -vVertex;\n"
        " eyeVec.x = dot(tmpVec, t);\n"
        " eyeVec.y = dot(tmpVec, b);\n"
        " eyeVec.z = dot(tmpVec, n);\n"
        " rotmat = mat3(t,b,n);"

        "const float LOG2 = 1.442695*0.07; \n"
        "gl_FogFragCoord = length(vVertex);\n"
        "fogFactor = exp2( -gl_Fog.density * gl_Fog.density * gl_FogFragCoord * gl_FogFragCoord * LOG2 ); \n "
        "fogFactor = clamp(fogFactor, 0.0, 1.0);\n"
        "}\n");

  frag = string("varying vec3 lightVec; \n"
        "varying vec3 eyeVec; \n"
        "varying vec2 texCoord; \n"
        "varying vec2 texCoord2; \n"
        "varying float fogFactor; \n"
        "varying mat3 rotmat; \n"
        "uniform sampler2D colorMap; \n"
        "uniform sampler2D normalMap; \n"
        "uniform float invRadius; \n"
        " \n"
        "void main (void)\n"
        "{ \n"
        "  float distSqr = dot(lightVec, lightVec);\n"
        "  float att = clamp(1.0 - invRadius * sqrt(distSqr), 0.0, 1.0); \n"
        "  vec3 lVec =lightVec* inversesqrt(distSqr);\n"
        " \n"
        "  vec3 vVec = normalize(eyeVec);\n"
        " \n"
        "  vec4 base = texture2D(colorMap, texCoord);\n"
        //"  if (base.w<0.2) discard;"
        "  vec3 bump = normalize( texture2D(normalMap, texCoord2).xyz * 2.0 - 1.0);\n"

        "  vec4 vAmbient = gl_LightSource[0].ambient * gl_FrontMaterial.ambient;\n"

        "  float diffuse = max( dot(lVec, bump), 0.0 );\n"
        "  vec3 rbump = normalize(rotmat*bump);"
        // all the next ones are really cool
        //" diffuse=pow(abs(diffuse), 0.25);\n"
        //" diffuse=pow(abs(diffuse), 10.0);\n"
        "  vec4 vDiffuse = gl_LightSource[0].diffuse * gl_FrontMaterial.diffuse * diffuse;  \n"

        "  float specular = pow(clamp(dot(reflect(-lVec, bump), vVec), 0.0, 1.0), gl_FrontMaterial.shininess ); \n"

        "  vec4 vSpecular = gl_LightSource[0].specular * gl_FrontMaterial.specular * specular;	\n"
        "  vec4 fcol = ( vAmbient*base + vDiffuse*base + vSpecular) * att; \n"
        "  fcol = mix(gl_Fog.color, fcol, fogFactor ); \n"
        "  gl_FragColor = fcol;"
        //"  if (rbump.z>0.0) gl_FragColor = fcol;"
        //"  else discard;"

        "}\n");
  NormalMapScale = 1.0;
}

void CBumpShader::Start() {
  //attribute vec3 vTangent;
  //"uniform sampler2D colorMap; \n"
  //"uniform sampler2D normalMap; \n"
  //"uniform float invRadius; \n"

  Shader->begin();
  Shader->setup_shader_texture(NormalMap, 0, "normalMap");
  Shader->setup_shader_texture(ColorMap, 1, "colorMap");
  InvertRadius = 0.0;
  Shader->sendUniform1f((char*)"invRadius", InvertRadius);
  Shader->sendUniform1f((char*)"NormalScale", NormalMapScale);
}

void CBumpShader::End() {
  Shader->disable_multitextures();
  Shader->end();
}

void CBumpShader::Tangent(CVector t) {
  t = t.normalize();
  Shader->sendAttribute3f((char*)"vTangent", t.x, t.y, t.z,tangent_loc);
}


void CBillboardLeafShader::Start() {
  Shader->begin();
  Shader->setup_shader_texture(leaftexture, 0, "colorMap");

}

void CBillboardLeafShader::End() {
  Shader->disable_multitextures();
  Shader->end();
}

CBillboardLeafShader::CBillboardLeafShader() {
  vert = string("varying vec2 texCoord;\n"
        "varying float fogFactor;\n"
        "varying vec3 vVertex; \n"
        "varying vec3 light;"
        "varying vec3 normal;"
        "varying vec4 color;"
                "void main(void) \n"
        "{ \n"
        " gl_Position = ftransform();\n"
        " texCoord = gl_MultiTexCoord0.xy;\n"
        "const float LOG2 = 1.442695*0.07; \n"
        "vVertex = vec3(gl_ModelViewMatrix * gl_Vertex); \n"
        "gl_FogFragCoord = length(vVertex);\n"
        "fogFactor = exp2( -gl_Fog.density * gl_Fog.density * gl_FogFragCoord * gl_FogFragCoord * LOG2 ); \n "
        "fogFactor = clamp(fogFactor, 0.0, 1.0);\n"
        "light = normalize(gl_LightSource[0].position).xyz;"
        "normal = normalize(gl_NormalMatrix  * gl_Normal); \n"
        "color = gl_Color;\n"
        "}\n");

  frag = string("varying vec2 texCoord;"
        "varying vec3 light;"
        "varying float fogFactor;"
        "varying vec3 vVertex; \n"
        "varying vec3 normal; \n"
        "varying vec4 color; \n"
        "uniform sampler2D colorMap; \n"
                "void main(void) \n"
        "{ \n"
        "  vec4 base = texture2D(colorMap, texCoord);\n"
        "  float alpha = base.w; "

        "  vec2 ln = (normalize(light)).xy;\n"
        "  float l = 2.5*(texCoord.y-0.3); //0.5*length(ln);\n"
        //"  float n = 0.8*(0.8 - dot(normal, normalize(light)));\n"
        "  float n = 0.8*(1.0 - dot(normal, normalize(light)));\n"
        "  n = clamp(n,0.8,1.5); \n"

        "  l = - 0.2 + 1.3*n*(length(light.xy - vec2(texCoord.x-0.5,-(texCoord.y-0.5))  )) + n*0.5;"
        "  l = clamp(l,0.4,1.5); \n"
        "  if (base.w<0.15) discard;"
        //		"  vec4 vSpecular = gl_LightSource[0].specular * gl_FrontMaterial.specular * specular;	\n"
        //		"  vec4 fcol = ( vAmbient*base + vDiffuse*base + vSpecular) * att; \n"
        "  vec4 vAmbient = gl_LightSource[0].diffuse * gl_FrontMaterial.diffuse;\n"
        "  vec4 fcol  = base*l*vAmbient*color; \n"
        "  fcol.w = alpha; \n"
        "  fcol = mix(gl_Fog.color, fcol, fogFactor ); \n"
        "  gl_FragColor = fcol;"
        "}\n");

}

void CMathShader::Start() {
  Shader->begin();

  if (Type!=MYBODY) {
    Shader->sendUniform1f((char*)"distort",Distort);
    Shader->sendUniform1f((char*)"type",Type);
    Shader->sendUniform1f((char*)"weight",Weight);
    Shader->setup_shader_texture(ScreenTexture2, 1, "screen2");
  }

  Shader->setup_shader_texture(ScreenTexture, 0, "screen");

}

void CMathShader::End() {
  Shader->disable_multitextures();
  Shader->end();
}

CMathShader::CMathShader() {
  vert = string("uniform sampler2D screen; \n"
        "uniform sampler2D screen2; \n"
        "uniform float weight; \n"
        "uniform float distort; \n"
        "uniform float type; \n"
                 "void main(void) \n"
        "{ \n"
        "   gl_TexCoord[0] = gl_MultiTexCoord0;\n"
        "   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
        //		"  color = gl_Color;\n"
        "}\n");

  mybody = string("uniform sampler2D screen; \n"
                "void main(void) \n"
        "{ \n"
        "  vec4 base = texture2D(screen, gl_TexCoord[0].st);\n"
        //"  vec4 base2 = texture2D(screen2, gl_TexCoord[0].st);\n"
        //"  base[3]=1.0;   \n"
        "  float d = 0.99*length(vec3(base)); \n"
        "  float l = 1.0/log(255.0); \n"
        "  d = 2.0*log(255.0*d)*l - 1.0; \n"
        "  float w = 10.0; \n"
        //"  if (d>1.0) d= 1.0; \n"
        "  base = vec4(0.0,0.0,0.0,1.0); \n"
        "  base[2] = exp( -w*(d-0.33)*(d-0.33)); \n "
        "  base[1] = exp( -w*(d-0.66)*(d-0.66)); \n "
        "  base[0] = exp( -w*(d-0.99)*(d-0.99)); \n "
        //		"  for (int i=0;i<3;i++) base[i] = log(255.0*base[i])*l; \n"
        //"   base[0] = 1.0; \n"
        "  gl_FragColor = base;"
        "}\n");

  frag = string("uniform sampler2D screen; \n"
        "uniform sampler2D screen2; \n"
        "uniform float weight; \n"
        "uniform float distort; \n"
        "uniform float type; \n"
        "float rand(vec2 co){"
        "return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);"
        "}"
        "void main(void) \n"
        "{ \n"

        " vec2 v = gl_TexCoord[0].st; \n"
        " v.x=0.0; \n"
        //" v.y=0.0; \n"
        " float d = rand(v)*distort - distort/2.0;"
        "  vec4 base = texture2D(screen, gl_TexCoord[0].st + vec2(d,0));\n"
        "  vec4 base2 = texture2D(screen2, gl_TexCoord[0].st+ vec2(d,0));\n"

        "  if (type==0.0) gl_FragColor = 2.0*(base*weight + base2*(1.0-weight)); "
        "  if (type==1.0) gl_FragColor = 4.0*(base*weight - base2*(1.0-weight)); "
        "  if (type==2.0) gl_FragColor = (base*weight * base2*(1.0-weight)); "
        "  if (type==3.0) gl_FragColor = 5.0*(base*weight / base2*(1.0-weight)); "
        "}\n");



}

void CMathShader::Render() {
  if (!Shader)
    return;
  COpenGLPointer->SetOrthographicProjection();

  COpenGLPointer->push();
  glLoadIdentity();

  glDisable(GL_DEPTH_TEST);
  glEnable(GL_TEXTURE_2D);
  glDisable(GL_FOG);
  glDisable(GL_CULL_FACE);

  Start();

  CVector c(1,1,1);
  glBegin(GL_QUADS);
  glTexCoord2d(0,1);
  glVertex3f(0,0,0);
  glTexCoord2d(1,1);
  glVertex3f(COpenGLPointer->window_width,0,0);
  glTexCoord2d(1,0);
  glVertex3f(COpenGLPointer->window_width,COpenGLPointer->window_height,0);
  glTexCoord2d(0,0);
  glVertex3f(0,COpenGLPointer->window_height,0);
  glEnd();

  End();

  COpenGLPointer->pop();

  COpenGLPointer->ResetPerspectiveProjection();

  glDisable(GL_TEXTURE_2D);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

}








void CThinLensShader::Start() {
  Shader->begin();

  Shader->setup_shader_texture(texture, 0, "screen");
  //cout << center.x << "   " << center.y << endl;
  Shader->sendUniform3f("center", center.x, center.y, center.z); //!< send vec3 to program
  Shader->sendUniform1f("lensStrength", lensStrength); //!< send vec3 to program
  Shader->sendUniform1f("Dd", Dd); //!< send vec3 to program
  Shader->sendUniform1f("type", type); //!< send vec3 to program


}

void CThinLensShader::End() {
  Shader->disable_multitextures();
  Shader->end();
}

CThinLensShader::CThinLensShader() {
  lensStrength = 0.3;
  Dd = 2.5;
  type = 0;
  vert = string(//"uniform sampler2D screen; \n"
        "uniform sampler2D texture; \n"
        " uniform vec3 center; \n"
                 "void main(void) \n"
        "{ \n"
        "   gl_TexCoord[0] = gl_MultiTexCoord0;\n"
        "   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
        //		"  color = gl_Color;\n"
        "}\n");

  frag = string(//"uniform sampler2D screen; \n"
        " uniform sampler2D texture; \n"
        " uniform vec3 center; \n"
        " uniform float type; \n"
        " uniform float lensStrength; \n"
        " float Ds = 1.0; \n"
        " uniform float Dd; \n"
        " float Dds = Dd + Ds; \n"
        " float getMass(vec2 p) {  \n"
        "   float m = 0.0;  \n"
        "   float l = length(p)*2.0;"
        "   if (abs(l)>1.0) return 0.0; \n "
        "   m = lensStrength*acos(l); \n"
                "   return m; } \n"

        " float getImageMass(vec2 p) {  \n"
        "   float m = 0.0;  \n"
        "   float l = clamp(length(texture2D(texture, p).xyz) - 0.05, 0.0, 1.0); \n "
        "   m = lensStrength*l; \n"
                "   return m; } \n"

        " float getGeneralMass(vec2 p) { \n"
        "   if (type==0.0) return getMass(p); \n"
                "   if (type==1.0) return getImageMass(p); \n "
                "  } \n"



        " float getPointSource(vec2 p) {  \n"
        "   float m = 0.0;  \n"
        "   float l = length(p)*2.0;"
        "   if (abs(l)<0.1) return 100.0; \n "
        "   m = 0.0; \n"
                "   return m; } \n"

          " vec2 getAlpha(vec2 org) {  \n "
        "  vec2 cen = center.xy*Dds + vec2(0.0, 0.0); \n "
          //"  vec2 org = vec2(rad.x*cos(rad.y), rad.x*sin(rad.y)); \n"
        "  float n = 1.0; \n"
        "  vec2 alpha = vec2(0,0); \n"
        "  int w = 20; float x = 1.0/float(w)/Dds; \n"
        "  if (length(org-cen)<1.0 + type*1.5)     \n"
        "  for (int i=0;i<w;i++) \n"
        "    for (int j=0;j<w;j++)   {\n"
        "       vec2 p = org + vec2(float(i-w/2)*x,float(j-w/2)*x); \n "
        "       //if (length(p-cen)<0.1) \n"
        "        {    \n "
        "         alpha = alpha - (p-cen)*getGeneralMass(cen - p);  \n"
        "         n+=1.0; \n "
        "       } \n"
        "  } \n"
        "  alpha = alpha/n; \n"
          //" if (length(alpha)<0.001) return vec2(0,0); \n"
          //"  return vec2(length(alpha), atan(alpha.y, alpha.x));; \n "
          "return alpha; \n"
        "  } \n"
        " \n"
                "void main(void) \n"
        "{ \n"
          /*		"  vec2 org = gl_TexCoord[0].st; \n "
        "  vec2 theta = vec2(length(org), atan(org.y, org.x)); \n"
        "  vec2 beta = theta - getAlpha(theta); \n"
        "  vec2 pos = vec2(beta.x*cos(beta.y), beta.x*sin(beta.y)); \n"
          */
        "  vec2 org = gl_TexCoord[0].st; \n "
        "  vec2 pos = org + getAlpha(Dds*org); \n"

          "  vec4 base = texture2D(texture, pos);\n"
        "  gl_FragColor = vec4(base.x,base.y,base.z,1.0); \n"
        "}\n");



}

void CThinLensShader::Render() {
  if (!Shader)
    return;
  COpenGLPointer->SetOrthographicProjection();

  COpenGLPointer->push();
  glLoadIdentity();

  glDisable(GL_DEPTH_TEST);
  glEnable(GL_TEXTURE_2D);
  glDisable(GL_FOG);
  glDisable(GL_CULL_FACE);

  Start();

  CVector c(1,1,1);
  glBegin(GL_QUADS);
  glTexCoord2d(0,1);
  glVertex3f(0,0,0);
  glTexCoord2d(1,1);
  glVertex3f(COpenGLPointer->window_width,0,0);
  glTexCoord2d(1,0);
  glVertex3f(COpenGLPointer->window_width,COpenGLPointer->window_height,0);
  glTexCoord2d(0,0);
  glVertex3f(0,COpenGLPointer->window_height,0);
  glEnd();

  End();

  COpenGLPointer->pop();

  COpenGLPointer->ResetPerspectiveProjection();

  glDisable(GL_TEXTURE_2D);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

}

/*
 3D lens shader
 */

void C3DLensShader::Start() {
  Shader->begin();

  Shader->setup_shader_texture3D(textureX, 0, "textureX");
  Shader->setup_shader_texture3D(textureY, 1, "textureY");
  Shader->setup_shader_texture3D(textureZ, 2, "textureZ");
  //cout << center.x << "   " << center.y << endl;
  Shader->sendUniform1f("lensStrength", lensStrength); //!< send vec3 to program
  Shader->sendUniform3f("Delta", delta.x, delta.y, delta.z); //!< send vec3 to program
  Shader->sendUniform3f("Camera", camera.x, camera.y, camera.z); //!< send vec3 to program
  Shader->sendUniform3f("Min", min.x, min.y, min.z); //!< send vec3 to program
  Shader->sendUniform1f("type", type); //!< send vec3 to program
  Shader->sendUniform1f("Width", width); //!< send vec3 to program
  Shader->sendUniform1f("Scale", scale); //!< send vec3 to program


}

void C3DLensShader::End() {
  Shader->disable_multitextures();
  Shader->end();
}

C3DLensShader::C3DLensShader() {
  lensStrength = 0.3;
  type = 0;
  vert = string(//"uniform sampler2D screen; \n"
        " uniform sampler3D textureX; \n"
        " uniform sampler3D textureY; \n"
        " uniform sampler3D textureZ; \n"
        " uniform float type; \n"
        " uniform vec3 Camera; \n"
        " uniform float Scale; \n"
        " uniform vec3 Min; \n"
        " uniform vec3 Delta; \n"
        " uniform float lensStrength; \n"
        " uniform float Width; \n"
        " varying float val;\n"
        " vec3 direction; \n"
        " vec3 rayPosition; \n"

        " float distanceTooPlane(in vec3 N, in vec3 PlaneP, in vec3 RayP) { \n"
        "   return dot(N, RayP - PlaneP); \n "
        "} \n"


        " float RayPlaneIntersection(in vec3 P0, in vec3 V, in vec3 N, in float d) { \n "
        "   return -(dot(P0, N) + d)/(dot(V,N)); \n "
        " } \n"

        " vec3 getAcceleration(in vec3 P) { \n"
        "   float f = 0.5; \n"
        "   float v = 1.0; \n"
        "   vec3 box_position = (P/Scale + Min*0.0)/(Width*Delta) + vec3(f, f, f); \n"
        "   if (box_position.x>v || box_position.x<0.0) return vec3(0,0,0); \n"
        "   if (box_position.y>v || box_position.y<0.0) return vec3(0,0,0); \n"
        "   if (box_position.z>v || box_position.z<0.0) return vec3(0,0,0); \n"
        //"   box_position = vec3(1,1,1) - box_position; \n"
        "   vec3 A = vec3(0,0,0); \n"
        "   A.x = texture3D(textureX, box_position.xzy).x; \n"
        "   A.y = texture3D(textureY, box_position.xzy).x; \n"
        "   A.z = texture3D(textureZ, box_position.xzy).x; \n"
        "   return A; \n "
        " } \n"

        " vec3 propagateRay() { \n"
        "   float n = 100.0; \n "
        "   vec3 V = direction/n; \n "
        "     float l = length(V); \n "
        "   for (int i=0;i<int(1.0*n);i++) { \n"
        "     rayPosition = rayPosition - V; \n "
        //"     V=normalize(V-getAcceleration(rayPosition)*0.07)*l; \n "
        "     V+=getAcceleration(rayPosition)*0.1; \n "
        "     float d = distanceTooPlane(normalize(direction), gl_Vertex.xyz, rayPosition); \n"
        //"     if (d<0.0) return rayPosition; \n "

        "   } \n"
        "   return rayPosition; \n"
        " } \n"

                " void main(void) \n"
        "{ \n"
        "   gl_TexCoord[0] = gl_MultiTexCoord0;\n"
        "   direction = ( Camera-gl_Vertex.xyz ); \n"
                "   rayPosition = Camera; \n"
        "   float distance = length(direction); \n"
                "   vec3 N = normalize(direction); \n"
        /*		"   float f = -0.5; \n"
        "   vec3 box_position = (gl_Vertex.xyz + Min)/(Width*Delta*Scale) - vec3(f, f, f); \n"
        "   vec4 v = texture3D(textureX, box_position); "
        "   val = v.x*1.0;  "
        */
        "  // float t = RayPlaneIntersection(rayPosition, direction, N, distance); \n"
        "  // vec3 P = rayPosition + t*direction; \n "
        "   vec3 P = propagateRay(); \n "
        "   gl_Position = gl_ModelViewProjectionMatrix * vec4(P, gl_Vertex.w);\n"
        //		"  color = gl_Color;\n"
        "}\n");

  frag = string(//"uniform sampler2D screen; \n"
        "varying float val; \n"
        "void main(void) \n"
        "{ \n"
          "  vec4 base = vec4(1.0, 1.0, 1.0, 1.0);\n"
        "  gl_FragColor = vec4(base.x,base.y,base.z,0.2); \n"
        "}\n");



}

void  C3DLensShader::Render() {
  if (!Shader)
    return;
  COpenGLPointer->SetOrthographicProjection();

  COpenGLPointer->push();
  glLoadIdentity();

  glDisable(GL_DEPTH_TEST);
  glEnable(GL_TEXTURE_2D);
  glDisable(GL_FOG);
  glDisable(GL_CULL_FACE);

  Start();

  CVector c(1,1,1);
  glBegin(GL_QUADS);
  glTexCoord2d(0,1);
  glVertex3f(0,0,0);
  glTexCoord2d(1,1);
  glVertex3f(COpenGLPointer->window_width,0,0);
  glTexCoord2d(1,0);
  glVertex3f(COpenGLPointer->window_width,COpenGLPointer->window_height,0);
  glTexCoord2d(0,0);
  glVertex3f(0,COpenGLPointer->window_height,0);
  glEnd();

  End();

  COpenGLPointer->pop();

  COpenGLPointer->ResetPerspectiveProjection();

  glDisable(GL_TEXTURE_2D);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

}
