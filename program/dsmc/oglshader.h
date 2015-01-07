#pragma once


//#define USE_WGLEXT          // if you dont want WGL Support, just remove this define
#ifdef _WIN32
#include <GL/glew.h>
#include <GL/glext.h>
#include <GL/gl.h>
#include <GL/glu.h>
#else
#include <GL/glew.h>

#ifdef LINUX
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#endif

#endif

#include <cutil.h>
#include <cvector.h>
#include <vector>
#include <iostream>
#include <string>

using namespace std;


// useful Macros:

#ifdef GLSL_WINDOWS
#define CLoadExtension(functype,funcname) ((funcname = (functype) wglGetProcAddress( #funcname )) == NULL)
#endif

#ifdef GLSL_LINUX
#define CLoadExtension(functype,funcname) ((funcname = (functype) glXGetProcAddressARB( #funcname )) == NULL)
#endif

// useful helper functions:
bool mcTestExtension(const char* extension_name);



// ***************************
// CGLExtensions - Helper Class
// ***************************

class CGLExtensions
{
public:
    CGLExtensions();
    ~CGLExtensions();

    void print(std::ostream& out=std::cout); //!< output list to ostream, standard: console
    bool checkextensions(char* extension_name);        //!< returns true if extension exists...

    void sort(void);                        //!< sort extensions (alphabetical)

    bool init(void);						//!< init all extensions

    std::vector<char*>   ExtensionList;     //!< List of all available OpenGL Extensions

private:
    void GLString_Convert(char* str);
    void addElement(char* str);

};

//-----------------------------------------------------------------------------

class CShaderProgram
{
    friend class CShaderObject;

public:
    CShaderProgram();
    ~CShaderProgram();

    int load(char* filename);   //!< read file, if result is 0 everything is ok. -1: File not found, -2: Empty File, -3: no memory
    void loadFromMemory(const char* program); //!< load program from char array, make sure program is 0 terminated!


    bool compile(void);         //!< compile program

    char* getCompilerLog(void);  //!< get compiler messages

    GLint GetUniLoc(const GLcharARB *name);      // get location of a variable

protected:

    int                 program_type;          //!< 1=Vertex Program, 2=Fragment Program, 0=none

    GLhandleARB         ProgramObject;         //!< Program Object
    GLubyte*            ShaderSource;          //!< ASCII Source-Code

    GLcharARB*          compiler_log;

    bool                is_compiled;            //!< true if compiled
    bool                _memalloc;               //!< true if shader allocated memory


};

// --------------------------------------------------------------

class CVertexShader : public CShaderProgram
{
  public:
       CVertexShader();
       ~CVertexShader();
};

// --------------------------------------------------------------

class CFragmentShader : public CShaderProgram
{
 public:
    CFragmentShader();
    ~CFragmentShader();

};

//-----------------------------------------------------------------------------

class CShaderObject
{
public:
    CShaderObject();            // Standard Constructor
    ~CShaderObject();           // Destructor

    void addShader(CShaderProgram* ShaderProgram); //!< add a Vertex or Fragment Program

    bool link(void);            //!< Link all Shaders
    char* getLinkerLog(void);   //!< get Linker messages

    void begin();	//!< use Shader. OpenGL calls will go through shader.
    void end();		//!< Stop using this shader. OpenGL calls will go through regular pipeline.

    bool oglslEnabled(void);    //!< returns true if OGLSL is enabled. It is possible user hardware doesn't support OGLSL!

    // Send Variables to Program

    void setup_shader_texture(GLuint texture, int texture_number, std::string variable_name);
    void setup_shader_texture3D(GLuint texture, int texture_number, std::string variable_name);
    void disable_multitextures();


    bool sendUniform1f(char* varname, GLfloat v0); //!< send float to program
    bool sendUniform2f(char* varname, GLfloat v0, GLfloat v1); //!< send vec2 to program
    bool sendUniform3f(char* varname, GLfloat v0, GLfloat v1, GLfloat v2); //!< send vec3 to program
    bool sendAttribute3f(char* varname, GLfloat v0, GLfloat v1, GLfloat v2, int slot); //!< send vec3 to program
    bool sendUniform4f(char* varname, GLfloat v0, GLfloat v1, GLfloat v2, GLfloat v3); //!< send vec4 to program

    bool sendUniform1i(char* varname, GLint v0);
    bool sendUniform2i(char* varname, GLint v0, GLint v1);
    bool sendUniform3i(char* varname, GLint v0, GLint v1, GLint v2);
    bool sendUniform4i(char* varname, GLint v0, GLint v1, GLint v2, GLint v3);

    bool sendUniform1fv(char* varname, GLsizei count, GLfloat *value);
    bool sendUniform2fv(char* varname, GLsizei count, GLfloat *value);
    bool sendUniform3fv(char* varname, GLsizei count, GLfloat *value);
    bool sendUniform4fv(char* varname, GLsizei count, GLfloat *value);

    bool sendUniform1iv(char* varname, GLsizei count, GLint *value);
    bool sendUniform2iv(char* varname, GLsizei count, GLint *value);
    bool sendUniform3iv(char* varname, GLsizei count, GLint *value);
    bool sendUniform4iv(char* varname, GLsizei count, GLint *value);

    bool sendUniformMatrix2fv(char* varname, GLsizei count, GLboolean transpose, GLfloat *value);
    bool sendUniformMatrix3fv(char* varname, GLsizei count, GLboolean transpose, GLfloat *value);
    bool sendUniformMatrix4fv(char* varname, GLsizei count, GLboolean transpose, GLfloat *value);

    bool sendUniform2ivArray(char* varname, GLsizei count, GLfloat *value);
    bool sendUniform4ivArray(char* varname, GLsizei count, GLfloat *value);

    // Receive Uniform variables:
    void GetUniformfv(char* name, GLfloat* values);
    void GetUniformiv(char* name, GLint* values);


    void manageMemory(void){_mM = true;}

    // Turn off all Shaders:
    static void useShader(bool b)		//!< Deactivate Shader
    {
        _noshader = b;
    }

    GLint GetUniLoc(const GLcharARB *name);      // get location of a variable
    GLint GetAttribLoc(const GLcharARB *name);      // get location of a variable

    GLhandleARB         ShaderObject;            // Shader Object

private:



    GLcharARB*          linker_log;
    bool                is_linked;
    std::vector<CShaderProgram*> ShaderList;     // List of all Shader Programs

    bool                _mM;
    static bool         _noshader;

};

//-----------------------------------------------------------------------------
// To simplify the process loading/compiling/linking shaders I created this
// high level interface to setup a vertex/fragment shader.

class CShaderManager
{
public:
    CShaderManager();
    ~CShaderManager();

    CShaderObject* loadfromFile(char* vertexFile, char* fragmentFile);    // load vertex/fragment shader from file
    CShaderObject* loadfromMemory(const char* vertexMem, const char* fragmentMem);

    //CShaderObject* arb_loadfromFile(char* vertexFile, char* fragmentFile);
    //CShaderObject* arb_loadfromMemory(const char* vertexMem, const char* fragmentMem);


    bool           free(CShaderObject* o);

private:
    std::vector<CShaderObject*>  _shaderObjectList;
};

struct CShader {
  CShaderObject* s;
  string name;
};


class CShaderContainer {
public:
   string shader_dir;
   CShaderManager shader_manager;
   // the shaders
   vector<CShader> shaders;

   CShaderObject* addshader(string name, string vert, string frag) {
     if (!glewIsSupported("GL_VERSION_2_0"))
       throw string("Error! Opengl 2.0 not supported!");

      vert = shader_dir + vert;
      frag = shader_dir + frag;
      CUtil::verify_file(vert);
      CUtil::verify_file(frag);

      CShader shader;
      shader.s = shader_manager.loadfromFile((char*)vert.c_str(), (char*)frag.c_str());
      shader.name = name;
      shaders.push_back(shader);
      return shader.s;
   }

   CShaderObject* addshaderfrommemory(string name, string vert, string frag) {
     if (!glewIsSupported("GL_VERSION_2_0"))
       throw string("Error! Opengl 2.0 not supported!");

     /*vert = shader_dir + vert;
      frag = shader_dir + frag;
      CUtil::verify_file(vert);
      CUtil::verify_file(frag);*/

      CShader shader;
      shader.s = shader_manager.loadfromMemory((char*)vert.c_str(), (char*)frag.c_str());
      shader.name = name;
      shaders.push_back(shader);
      return shader.s;
   }

   CShaderObject* getshader(string name) {
      for (unsigned int i=0;i<shaders.size();i++) {
         if (shaders[i].name == name)
            return shaders[i].s;
          }
      return 0;
   }

   CShaderContainer() {
      shader_dir = "shaders/";

   }
   void initialize() {

   }


};
