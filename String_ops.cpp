#include"String_ops.h"
bool compare(string original, string input)
{
 unsigned int i;
 char c;
 bool res=false;
 string str1,str2;
 str1=original;
 str2=input;
 for(i=0;i<str1.length();i++)
 {
  c=str1[i];
  str1[i]=tolower(c);
 }
 for(i=0;i<str2.length();i++)
 {
  c=str2[i];
  str2[i]=tolower(c);
 }
 if(str1==str2){res=true;}
 return res;
}
void lowercase(string &original)
{
 unsigned int i;
 char c;
 for(i=0;i<original.length();i++)
 {
  c=original[i];
  original[i]=tolower(c);
 }

}
