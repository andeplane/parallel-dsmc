#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

class CUtil {
  
 public:
  static void Tokenize(const string& str,
		       vector<string>& tokens,
		       const string& delimiters = " ");
  
  
  static void string2char(string s, char* to);
  static string toString(double d, string param); 
  static string toString(double d);
  static string toString(int d);
  static const char* read_textfile(string filename);
  static void verify_file(string filename);
  static bool verify_file_bool(string filename);
  static string trim(string s);

};

