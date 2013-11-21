#include <oglshader.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string.h>

//#include "aGL_Extensions.cpp"

using namespace std;

bool CuseGLSL = false;


//-----------------------------------------------------------------------------
// Error Strings
char* CaGLSLErrorString[] = {
        "(e0000) GLSL not enabled",
        "(e0001) not a valid program object",
        "(e0002) not a valid object",
        "(e0003) out of memory",
        "(e0004) unknown compiler error"};
//-----------------------------------------------------------------------------

// GL ERROR CHECK
int CCheckGLError(char *file, int line)
{
    GLenum glErr;
    int    retCode = 0;

    glErr = glGetError();
    while (glErr != GL_NO_ERROR) {
        cout << "GL Error #" << glErr << "(" << gluErrorString(glErr) << ") " << " in File " << file << " at line: " << line << endl;
        retCode = 1;
        glErr = glGetError();
    }
    return retCode;
}
#define CHECK_GL_ERROR() CCheckGLError(__FILE__, __LINE__)


//-----------------------------------------------------------------------------

bool CinitGLSL(void)
{
 int error = 0;

  if (CuseGLSL) return true;  // already initialized

/*   if (!init_ARB_fragment_shader()) error = -1;
   if (!init_ARB_vertex_shader()) error = -1;;
   if (!init_ARB_shader_objects()) error = -1;;
 */

  if (error)
  {
    CuseGLSL = false;
  }
  else
  {
    CuseGLSL = true;
  }

return CuseGLSL;
}

//-----------------------------------------------------------------------------

// Test if a certain Extension exists, returns true if it exists
// example: mcTextExtension("ARB_vertex_shader");
bool CmcTestExtension(const char* extension_name)
{
    const GLubyte *str;
    str = glGetString(GL_EXTENSIONS);

    if (strstr((const char *)str, extension_name) != 0)  return true;

#ifdef GLSL_WINDOWS
#ifdef USE_WGLEXT
    PFNWGLGETEXTENSIONSSTRINGARBPROC   wglGetExtensionsStringARB;
    wglGetExtensionsStringARB = (PFNWGLGETEXTENSIONSSTRINGARBPROC)wglGetProcAddress("wglGetExtensionsStringARB");

    if(wglGetExtensionsStringARB == NULL)
    {
       return false;
    }

    str = (const GLubyte*) wglGetExtensionsStringARB(wglGetCurrentDC());
    if (strstr((const char *)str, extension_name) != 0) return true;
#endif
#else
// ignore additional platform extensions
#endif

    return false;
}

// ************************************************************************
// CGLExtensions : Create a STL Vector with all available CGLExtensions names
// ************************************************************************

CGLExtensions::CGLExtensions()
{
const GLubyte *str;

    str = glGetString(GL_EXTENSIONS);
    GLString_Convert((char*)str);

#ifdef GLSL_WINDOWS
#ifdef USE_WGLEXT
    PFNWGLGETEXTENSIONSSTRINGARBPROC   wglGetExtensionsStringARB;
    wglGetExtensionsStringARB = (PFNWGLGETEXTENSIONSSTRINGARBPROC)wglGetProcAddress("wglGetExtensionsStringARB");

    if(wglGetExtensionsStringARB == NULL)
    {
       return;
    }

    str = (const GLubyte*) wglGetExtensionsStringARB(wglGetCurrentDC());
    GLString_Convert((char*)str);
#endif
#else
// ignore additional extensions....
#endif

}

//-----------------------------------------------------------------------------

CGLExtensions::~CGLExtensions()
{
 for (unsigned int i=0;i<ExtensionList.size();i++)
 {
    delete[] ExtensionList[i];  // free string
 }
}

//-----------------------------------------------------------------------------

void CGLExtensions::GLString_Convert(char* str)
{
#define ELaAX_BUF_SIZE 80
char buf[ELaAX_BUF_SIZE];
int i=0;

   char* d = (char*) str;
        while (*d!=0)
        {
            if (*d==' ')
            {
                buf[i] = 0;
                i=0;
                addElement(buf);
            }
            else
            {
                buf[i++]= *d;
                if (i>ELaAX_BUF_SIZE-1) i=ELaAX_BUF_SIZE-1;
            }
            d++;
        }
}

//-----------------------------------------------------------------------------

void CGLExtensions::addElement(char* str)
{
  char* s = new char[strlen(str)+1];
#ifdef _WIN32
  strcpy_s(s, strlen(str)+1,str);
#else
  strcpy(s,str);
#endif
  ExtensionList.push_back(s);

}

//-----------------------------------------------------------------------------

void CGLExtensions::print(ostream& out)
{
   out << "OpenGL Vendor: " << (char*) glGetString(GL_VENDOR) << "\n";
   out << "OpenGL Renderer: " << (char*) glGetString(GL_RENDERER) << "\n";
   out << "OpenGL Version: " << (char*) glGetString(GL_VERSION) << "\n\n";

   for (unsigned int i=0;i<ExtensionList.size();i++)
   {
        out << ExtensionList[i] << endl;
   }

}

//-----------------------------------------------------------------------------

bool CGLExtensions::checkextensions(char* extension_name){
   for (unsigned int i=0;i<ExtensionList.size();i++)
   {
        if (strstr((const char *)ExtensionList[i], extension_name) != 0)
        {
            return true;
        }
   }

   return false;
}

//-----------------------------------------------------------------------------

bool CGLExtensions::init(void)
{
    //init_ARB_extensions();
    return true;
}

//-----------------------------------------------------------------------------

// Bubble Sort to sort extension names.
void CGLExtensions::sort(void)
{
char* tmp;

     for (unsigned int i = 0;  i < ExtensionList.size(); i++)
     {
        for (unsigned int j = (unsigned int) ExtensionList.size()-1; j>i; j--)
        {
            if (strcmp(ExtensionList[j-1],ExtensionList[j]) > 0)
            {
                tmp = ExtensionList[j-1];
                ExtensionList[j-1] = ExtensionList[j];
                ExtensionList[j] = tmp;
            }
        }
    }
}


// ************************************************************************
// Implementation der Klasse CShaderObject
// ************************************************************************


CShaderObject::CShaderObject()
{
  ShaderObject = 0;
  linker_log = 0;
  _mM = false;


  if (!CinitGLSL())
  {
    cout << "Error initializing OpenGL Shading Language function pointers" << endl;
  }

  if (!CmcTestExtension("GL_ARB_shader_objects"))
        cout << "**warning** GL_ARB_shader_objects not defined!!\n";
  if (!CmcTestExtension("GL_ARB_vertex_shader"))
        cout << "**warning** GL_ARB_vertex_shader not defined!!\n";
  if (!CmcTestExtension("GL_ARB_fragment_shader"))
        cout << "**warning** GL_ARB_fragment_shader not defined!!\n";

  if (!CmcTestExtension("GL_ARB_shading_language_100"))
        cout << "**warning** GL_ARB_shading_language_100 not defined!!\n";

  if (CinitGLSL())
  {
      ShaderObject = glCreateProgramObjectARB();
  }
  is_linked = false;
}

//-----------------------------------------------------------------------------

CShaderObject::~CShaderObject()
{
    if (linker_log!=0) free(linker_log);
    if (CuseGLSL)
    {
       for (unsigned int i=0;i<ShaderList.size();i++)
       {
            glDetachObjectARB(ShaderObject, ShaderList[i]->ProgramObject);
            CHECK_GL_ERROR(); // if you get an error here, you deleted the Program object first and then
                           // the ShaderObject! Always delete ShaderObjects last!
            if (_mM) delete ShaderList[i];
       }

       glDeleteObjectARB(ShaderObject);
       CHECK_GL_ERROR();
    }

}

//-----------------------------------------------------------------------------

bool CShaderObject::oglslEnabled(void)
{
   return CuseGLSL;
}

//-----------------------------------------------------------------------------

void CShaderObject::addShader(CShaderProgram* ShaderProgram)
{
if (!CuseGLSL) return;

   if (ShaderProgram==0) return;


   if (!ShaderProgram->is_compiled)
   {
        cout << "**warning** please compile program before adding object! trying to compile now...\n";
        if (!ShaderProgram->compile())
        {
            cout << "...compile ERROR!\n";
            return;
        }
        else
        {
            cout << "...ok!\n";
        }
   }

   ShaderList.push_back(ShaderProgram);

}
// quick set-up multi texturing

void CShaderObject::setup_shader_texture(GLuint texture, int texture_number, string variable_name) {
        GLuint txt = GL_TEXTURE0_ARB;
        switch (texture_number) {
               case 0: txt = GL_TEXTURE0_ARB; break;
               case 1: txt = GL_TEXTURE1_ARB; break;
               case 2: txt = GL_TEXTURE2_ARB; break;
               case 3: txt = GL_TEXTURE3_ARB; break;
               case 4: txt = GL_TEXTURE4_ARB; break;
               case 5: txt = GL_TEXTURE5_ARB; break;
               case 6: txt = GL_TEXTURE6_ARB; break;
               case 7: txt = GL_TEXTURE7_ARB; break;
        }

        GLuint ref = GetUniLoc(variable_name.c_str());
        glUniform1iARB(ref, texture_number);
        glActiveTextureARB(txt);
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, texture);
}

void CShaderObject::setup_shader_texture3D(GLuint texture, int texture_number, string variable_name) {
        GLuint txt = GL_TEXTURE0_ARB;
        switch (texture_number) {
               case 0: txt = GL_TEXTURE0_ARB; break;
               case 1: txt = GL_TEXTURE1_ARB; break;
               case 2: txt = GL_TEXTURE2_ARB; break;
               case 3: txt = GL_TEXTURE3_ARB; break;
               case 4: txt = GL_TEXTURE4_ARB; break;
               case 5: txt = GL_TEXTURE5_ARB; break;
               case 6: txt = GL_TEXTURE6_ARB; break;
               case 7: txt = GL_TEXTURE7_ARB; break;
        }

        GLuint ref = GetUniLoc(variable_name.c_str());
        glUniform1iARB(ref, texture_number);
        glActiveTextureARB(txt);
        glEnable(GL_TEXTURE_3D);
        glBindTexture(GL_TEXTURE_3D, texture);
}



void CShaderObject::disable_multitextures() {
       glActiveTexture(GL_TEXTURE1);
       glDisable(GL_TEXTURE_2D);
       glActiveTexture(GL_TEXTURE2);
       glDisable(GL_TEXTURE_2D);
       glActiveTexture(GL_TEXTURE3);
       glDisable(GL_TEXTURE_2D);
       glActiveTexture(GL_TEXTURE4);
       glDisable(GL_TEXTURE_2D);
       glActiveTexture(GL_TEXTURE5);
       glDisable(GL_TEXTURE_2D);
       glActiveTexture(GL_TEXTURE6);
       glDisable(GL_TEXTURE_2D);
       glActiveTexture(GL_TEXTURE7);
       glDisable(GL_TEXTURE_2D);
       glActiveTexture(GL_TEXTURE8);
       glDisable(GL_TEXTURE_2D);
       glActiveTexture(GL_TEXTURE9);
       glDisable(GL_TEXTURE_2D);
       glActiveTexture(GL_TEXTURE10);
       glDisable(GL_TEXTURE_2D);


       glActiveTexture(GL_TEXTURE0);
       glDisable(GL_TEXTURE_2D);

}

//-----------------------------------------------------------------------------

bool CShaderObject::link(void)
{
if (!CuseGLSL) return false;

unsigned int i;

    if (is_linked)  // already linked, detach everything first
    {
       cout << "**warning** Object is already linked, trying to link again" << endl;
       for (i=0;i<ShaderList.size();i++)
       {
            glDetachObjectARB(ShaderObject, ShaderList[i]->ProgramObject);
            CHECK_GL_ERROR();
       }
    }

    for (i=0;i<ShaderList.size();i++)
    {
        glAttachObjectARB(ShaderObject, ShaderList[i]->ProgramObject);
        CHECK_GL_ERROR();
        //cout << "attaching ProgramObj [" << i << "] @ 0x" << hex << ShaderList[i]->ProgramObject << " in ShaderObj @ 0x"  << ShaderObject << endl;
    }

    int linked;
    glLinkProgramARB(ShaderObject);
    CHECK_GL_ERROR();
    glGetObjectParameterivARB(ShaderObject, GL_OBJECT_LINK_STATUS_ARB, &linked);
    CHECK_GL_ERROR();

    if (linked)
    {
        is_linked = true;
        return true;
    }
    else
    {
        cout << "**linker error**\n";
    }

return false;
}

//-----------------------------------------------------------------------------
// Compiler Log: Ausgabe der Compiler Meldungen in String

char* CShaderObject::getLinkerLog(void)
{
if (!CuseGLSL) return CaGLSLErrorString[0];

 int blen = 0;
 int slen = 0;


 if (ShaderObject==0) return CaGLSLErrorString[2];

 glGetObjectParameterivARB(ShaderObject, GL_OBJECT_INFO_LOG_LENGTH_ARB , &blen);
 CHECK_GL_ERROR();

 if (blen > 1)
 {
    if (linker_log!=0)
    {
        free(linker_log);
        linker_log =0;
    }
    if ((linker_log = (GLcharARB*)malloc(blen)) == NULL)
     {
        printf("ERROR: Could not allocate compiler_log buffer\n");
        return CaGLSLErrorString[3];
    }

     glGetInfoLogARB(ShaderObject, blen, &slen, linker_log);
     CHECK_GL_ERROR();

 }
 if (linker_log!=0)
    return (char*) linker_log;

    return CaGLSLErrorString[4];
}

void CShaderObject::begin(void)
{
try {
if (!CuseGLSL) return;
if (ShaderObject == 0) return;
if (!_noshader) return;

    if (is_linked)
    {
        glUseProgramObjectARB(ShaderObject);
        //CHECK_GL_ERROR();
    }
} catch(...) { throw string("Error using shader!"); }

}

//-----------------------------------------------------------------------------

void CShaderObject::end(void)
{
if (!CuseGLSL) return;
if (!_noshader) return;


    glUseProgramObjectARB(0);
    //CHECK_GL_ERROR();
}

//-----------------------------------------------------------------------------

bool CShaderObject::sendUniform1f(char* varname, GLfloat v0)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform1fARB(loc, v0);

    return true;
}

//-----------------------------------------------------------------------------

bool CShaderObject::sendUniform2f(char* varname, GLfloat v0, GLfloat v1)
{
   if (!CuseGLSL) return false; // GLSL not available
   if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform2fARB(loc, v0, v1);

    return true;
}

//-----------------------------------------------------------------------------

bool CShaderObject::sendUniform3f(char* varname, GLfloat v0, GLfloat v1, GLfloat v2)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform3fARB(loc, v0, v1, v2);

    return true;
}

bool CShaderObject::sendAttribute3f(char* varname, GLfloat v0, GLfloat v1, GLfloat v2, int slot)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    //int loc = glGetAttribLocation(ShaderObject,varname);
    //cout << loc << endl;
    glVertexAttrib3fARB(slot, v0, v1, v2);

    return true;
}



//-----------------------------------------------------------------------------

bool CShaderObject::sendUniform4f(char* varname, GLfloat v0, GLfloat v1, GLfloat v2, GLfloat v3)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform4fARB(loc, v0, v1, v2, v3);

    return true;
}


bool CShaderObject::sendUniform2ivArray(char* varname, GLsizei count, GLfloat *value)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform2fvARB(loc,count, value);
    //glUniform2fvARB(uniform_location, 25, offsets);

    return true;
}

bool CShaderObject::sendUniform4ivArray(char* varname, GLsizei count, GLfloat *value)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform4fvARB(loc,count, value);
    //glUniform2fvARB(uniform_location, 25, offsets);

    return true;
}


//-----------------------------------------------------------------------------

bool CShaderObject::sendUniform1i(char* varname, GLint v0)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform1iARB(loc, v0);

    return true;
}
bool CShaderObject::sendUniform2i(char* varname, GLint v0, GLint v1)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform2iARB(loc, v0, v1);


    return true;
}

//-----------------------------------------------------------------------------

bool CShaderObject::sendUniform3i(char* varname, GLint v0, GLint v1, GLint v2)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform3iARB(loc, v0, v1, v2);

    return true;
}
bool CShaderObject::sendUniform4i(char* varname, GLint v0, GLint v1, GLint v2, GLint v3)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform4iARB(loc, v0, v1, v2, v3);

    return true;
}

//-----------------------------------------------------------------------------

bool CShaderObject::sendUniform1fv(char* varname, GLsizei count, GLfloat *value)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform1fvARB(loc, count, value);

    return true;
}
bool CShaderObject::sendUniform2fv(char* varname, GLsizei count, GLfloat *value)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform2fvARB(loc, count, value);

    return true;
}

//-----------------------------------------------------------------------------

bool CShaderObject::sendUniform3fv(char* varname, GLsizei count, GLfloat *value)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform3fvARB(loc, count, value);

    return true;
}

//-----------------------------------------------------------------------------

bool CShaderObject::sendUniform4fv(char* varname, GLsizei count, GLfloat *value)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform4fvARB(loc, count, value);

    return true;
}

//-----------------------------------------------------------------------------

bool CShaderObject::sendUniform1iv(char* varname, GLsizei count, GLint *value)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform1ivARB(loc, count, value);

    return true;
}

//-----------------------------------------------------------------------------

bool CShaderObject::sendUniform2iv(char* varname, GLsizei count, GLint *value)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform2ivARB(loc, count, value);

    return true;
}

//-----------------------------------------------------------------------------

bool CShaderObject::sendUniform3iv(char* varname, GLsizei count, GLint *value)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform3ivARB(loc, count, value);

    return true;
}

//-----------------------------------------------------------------------------

bool CShaderObject::sendUniform4iv(char* varname, GLsizei count, GLint *value)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniform4ivARB(loc, count, value);

    return true;
}

//-----------------------------------------------------------------------------

bool CShaderObject::sendUniformMatrix2fv(char* varname, GLsizei count, GLboolean transpose, GLfloat *value)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniformMatrix2fvARB(loc, count, transpose, value);

    return true;
}

//-----------------------------------------------------------------------------

bool CShaderObject::sendUniformMatrix3fv(char* varname, GLsizei count, GLboolean transpose, GLfloat *value)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniformMatrix3fvARB(loc, count, transpose, value);

    return true;
}

//-----------------------------------------------------------------------------

bool CShaderObject::sendUniformMatrix4fv(char* varname, GLsizei count, GLboolean transpose, GLfloat *value)
{
    if (!CuseGLSL) return false; // GLSL not available
    if (!_noshader) return true;

    GLint loc = GetUniLoc(varname);
    if (loc==-1) return false;  // can't find variable

    glUniformMatrix4fvARB(loc, count, transpose, value);

    return true;
}


//-----------------------------------------------------------------------------

GLint CShaderObject::GetUniLoc(const GLcharARB *name)
{
    GLint loc;

    loc = glGetUniformLocationARB(ShaderObject, name);
    if (loc == -1)
    {
        cout << "Error: can't find uniform variable \"" << name << "\"\n";
    }
    //CHECK_GL_ERROR();
    return loc;
}


GLint CShaderObject::GetAttribLoc(const GLcharARB *name)
{
    GLint loc = 0;

    //loc = glBindAttribLocationARB(ShaderObject, name);
    //loc = glBindAttribLocationARB(ShaderObject, name);
    if (loc == -1)
    {
        cout << "Error: can't find attrib variable \"" << name << "\"\n";
    }
    //CHECK_GL_ERROR();
    return loc;
}

//-----------------------------------------------------------------------------

void CShaderObject::GetUniformfv(char* name, GLfloat* values)
{
if (!CuseGLSL) return;
    GLint loc;

    loc = glGetUniformLocationARB(ShaderObject, name);
    if (loc == -1)
    {
        cout << "Error: can't find uniform variable \"" << name << "\"\n";
    }
    glGetUniformfvARB(ShaderObject, loc, values);

}

//-----------------------------------------------------------------------------

void CShaderObject::GetUniformiv(char* name, GLint* values)
{
if (!CuseGLSL) return;


    GLint loc;

    loc = glGetUniformLocationARB(ShaderObject, name);
    if (loc == -1)
    {
        cout << "Error: can't find uniform variable \"" << name << "\"\n";
    }

    glGetUniformivARB(ShaderObject, loc, values);

}

bool  CShaderObject::_noshader = true;

// ************************************************************************
// Shader Program : Manage Shader Programs (Vertex/Fragment)
// ************************************************************************


CShaderProgram::CShaderProgram()
{
    CinitGLSL();
    compiler_log = 0;
    is_compiled = false;
    program_type = 0;
    ProgramObject = 0;
    ShaderSource = 0;
    _memalloc = false;

}

//-----------------------------------------------------------------------------

CShaderProgram::~CShaderProgram()
{
   if (compiler_log!=0) free(compiler_log);
   if (ShaderSource!=0)
   {
        if (_memalloc)
            delete[] ShaderSource;  // free ASCII Source
   }

   if (is_compiled)
   {
        glDeleteObjectARB(ProgramObject);
        CHECK_GL_ERROR();
   }
}

//-----------------------------------------------------------------------------
unsigned long CgetFileLength(ifstream& file)
{
    if(!file.good()) return 0;

    //unsigned long pos=file.tellg();
    file.seekg(0,ios::end);
    unsigned long len = (unsigned long)file.tellg();
    file.seekg(ios::beg);

    return len;
}


//-----------------------------------------------------------------------------

int CShaderProgram::load(char* filename)
{
    ifstream file;
    file.open(filename, ios::in);
    if(!file) return -1;

    unsigned long len = CgetFileLength(file);

    if (len==0) return -2;   // "Empty File"

    if (ShaderSource!=0)    // there is already a source loaded, free it!
    {
        if (_memalloc)
        delete[] ShaderSource;
    }

    ShaderSource = (GLubyte*) new char[len+1];
    if (ShaderSource == 0) return -3;   // can't reserve memory
    _memalloc = true;


    ShaderSource[len] = 0;  // len isn't always strlen cause some characters are stripped in ascii read...
                            // it is important to 0-terminate the real length later, len is just max possible value...

    unsigned int i=0;
    while (file.good())
    {
        ShaderSource[i++] = file.get();       // get character from file
        if (i>len) i=len;   // coding guidelines...
    }

    ShaderSource[i] = 0;  // 0 terminate it.

    file.close();

return 0;
}

//-----------------------------------------------------------------------------

void CShaderProgram::loadFromMemory(const char* program)
{
    if (ShaderSource!=0)    // there is already a source loaded, free it!
    {
        if (_memalloc)
        delete[] ShaderSource;
    }
   _memalloc = false;
   ShaderSource = (GLubyte*) program;

}


// ----------------------------------------------------
// Compiler Log: Ausgabe der Compiler Meldungen in String

char* CShaderProgram::getCompilerLog(void)
{
if (!CuseGLSL) return CaGLSLErrorString[0];

 int blen = 0;
 int slen = 0;


 if (ProgramObject==0) return CaGLSLErrorString[1];

 glGetObjectParameterivARB(ProgramObject, GL_OBJECT_INFO_LOG_LENGTH_ARB , &blen);
 CHECK_GL_ERROR();

 if (blen > 1)
 {
    if (compiler_log!=0)
    {
        free(compiler_log);
        compiler_log =0;
    }
    if ((compiler_log = (GLcharARB*)malloc(blen)) == NULL)
     {
        printf("ERROR: Could not allocate compiler_log buffer\n");
        return CaGLSLErrorString[3];
    }

     glGetInfoLogARB(ProgramObject, blen, &slen, compiler_log);
     CHECK_GL_ERROR();
     //cout << "compiler_log: \n", compiler_log);
 }
 if (compiler_log!=0)
    return (char*) compiler_log;

    return CaGLSLErrorString[4];
}

// ----------------------------------------------------

bool CShaderProgram::compile(void)
{
if (!CuseGLSL) return false;

is_compiled = false;

int compiled = 0;

  if (ShaderSource==0) return false;

  GLint	length = (GLint) strlen((const char*)ShaderSource);
  glShaderSourceARB(ProgramObject, 1, (const GLcharARB **)&ShaderSource, &length);
  CHECK_GL_ERROR();

  glCompileShaderARB(ProgramObject);
  CHECK_GL_ERROR();
  glGetObjectParameterivARB(ProgramObject, GL_OBJECT_COMPILE_STATUS_ARB, &compiled);
  CHECK_GL_ERROR();

  if (compiled) is_compiled=true;

return is_compiled;
}

// ----------------------------------------------------


CVertexShader::CVertexShader()
{
  program_type = 1;
   if (CuseGLSL)
   {
       ProgramObject = glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);
       CHECK_GL_ERROR();
   }
}

// ----------------------------------------------------

CVertexShader::~CVertexShader()
{
}

// ----------------------------------------------------

CFragmentShader::CFragmentShader()
{
    program_type = 2;
    if (CuseGLSL)
    {
        ProgramObject = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);
        CHECK_GL_ERROR();
    }
}

// ----------------------------------------------------

CFragmentShader::~CFragmentShader()
{
}

// ----------------------------------------------------------------------------
// ShaderManager: Easy use of (multiple) Shaders

CShaderManager::CShaderManager()
{

}

CShaderManager::~CShaderManager()
{
   // free objects
   vector<CShaderObject*>::iterator  i=_shaderObjectList.begin();
   while (i!=_shaderObjectList.end())
   {
        CShaderObject* o = *i;
        i=_shaderObjectList.erase(i);
        delete o;
   }
}

// ----------------------------------------------------------------------------

CShaderObject* CShaderManager::loadfromFile(char* vertexFile, char* fragmentFile)
{
   CShaderObject* o = new CShaderObject();

   CVertexShader* tVertexShader = new CVertexShader;
   CFragmentShader* tFragmentShader = new CFragmentShader;

    // load vertex program
   if (vertexFile!=0)
   if (tVertexShader->load(vertexFile) != 0)
   {
     cout << "error: can't load vertex shader!\n";
     delete o;
     delete tVertexShader;
     delete tFragmentShader;
     return 0;
   }

  // Load fragment program
  if (fragmentFile!=0)
  if (tFragmentShader->load(fragmentFile) != 0)
  {
     cout << "error: can't load fragment shader!\n";
     delete o;
     delete tVertexShader;
     delete tFragmentShader;
     return 0;
  }

  // Compile vertex program
  if (vertexFile!=0)
  if (!tVertexShader->compile())
  {
      cout << "\n\n***SEVERE COMPILER ERROR (Vertex Shader):\n";
      cout << tVertexShader->getCompilerLog() << endl;
      delete o;
      delete tVertexShader;
      delete tFragmentShader;
      throw string("Vertex shader compile error in '" + string(vertexFile)+"'");
      return 0;
  }

  // Compile fragment program
  if (fragmentFile!=0)
  if (!tFragmentShader->compile())
  {
     cout << "\n\n***SEVERE COMPILER ERROR (Fragment Shader):\n";
     cout << tFragmentShader->getCompilerLog() << endl;

     delete o;
     delete tVertexShader;
     delete tFragmentShader;
     throw string("Fragment shader compile error in '" + string(fragmentFile) + "'");
     return 0;

  }

  // Add to object
  if (vertexFile!=0) o->addShader(tVertexShader);
  if (fragmentFile!=0) o->addShader(tFragmentShader);

  // link
  if (!o->link())
  {
     cout << "**LINKER ERROR\n";
     cout << o->getLinkerLog() << endl;
     delete o;
     delete tVertexShader;
     delete tFragmentShader;
     return 0;
  }

  _shaderObjectList.push_back(o);
  o->manageMemory();

   return o;
}

// ----------------------------------------------------------------------------

CShaderObject* CShaderManager::loadfromMemory(const char* vertexMem, const char* fragmentMem)
{
  CShaderObject* o = new CShaderObject();

  CVertexShader* tVertexShader = new CVertexShader;
  CFragmentShader* tFragmentShader = new CFragmentShader;

  // get vertex program
  if (vertexMem!=0)
     tVertexShader->loadFromMemory(vertexMem);

  // get fragment program
  if (fragmentMem!=0)
     tFragmentShader->loadFromMemory(fragmentMem);

  // Compile vertex program
  if (vertexMem!=0)
  if (!tVertexShader->compile())
  {
      cout << "***COMPILER ERROR (Vertex Shader):\n";
      cout << tVertexShader->getCompilerLog() << endl;
      delete o;
      delete tVertexShader;
      delete tFragmentShader;
      return 0;
  }

  // Compile fragment program
  if (fragmentMem!=0)
  if (!tFragmentShader->compile())
  {
     cout << "***COMPILER ERROR (Fragment Shader):\n";
     cout << tFragmentShader->getCompilerLog() << endl;

     delete o;
     delete tVertexShader;
     delete tFragmentShader;
     return 0;

  }

  // Add to object
  if (vertexMem!=0) o->addShader(tVertexShader);
  if (fragmentMem!=0) o->addShader(tFragmentShader);

  // link
  if (!o->link())
  {
     cout << "**LINKER ERROR\n";
     cout << o->getLinkerLog() << endl;
     delete o;
     delete tVertexShader;
     delete tFragmentShader;
     return 0;
  }

  _shaderObjectList.push_back(o);
  o->manageMemory();

   return o;
}

// ----------------------------------------------------------------------------

 bool  CShaderManager::free(CShaderObject* o)
 {
   vector<CShaderObject*>::iterator  i=_shaderObjectList.begin();
   while (i!=_shaderObjectList.end())
   {
        if ((*i)==o)
        {
            _shaderObjectList.erase(i);
            delete o;
            return true;
        }
        i++;
   }
   return false;
 }

// ----------------------------------------------------------------------------


