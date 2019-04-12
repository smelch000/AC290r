#ifndef READARGS_H
#define READARGS_H
 
/* 
  18/10/2016   P.Miocchi
  
WARNING: templates methods must not go into a .cpp file
*/


#include <string>
#include <cstring>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>



/*

Read options values from the command line.

GetParameter gives false if parameter has not been given
or of wrong type.

Example of usage:

#include <iostream>
#include <string>
#include "readargs.h"

int main(int argc, char *argv[]){
  int numargs;
  string options[]={"I", "r", "N", "i"};
  string descriptions[]={"Input file(s)=","r_core=", "Num. particles=","simulation Id="};

  ReadArgs inputPars(options,descriptions,4,numargs,argc, argv,&cout);
  cout << numargs << endl;
  int npar;
  string files[10];
  double rcore;
  string descr, id;

  if (!inputPars.isset("N")) {
    cerr << "Error: option -N not given!" << endl;
  }
  nfiles=inputPars.GetParameter<string>("I",descr,files,10);
  if (nfiles>0) {
    cout << descr;
    for (int i=0; i<nfiles;i++) cout << files[i]<<" ";
    cout << endl;
  }
  if (inputPars.GetParameter<double>("r",descr,rcore)) cout << descr << rcore << endl;
  if (inputPars.GetParameter<string>("i",descr,id)) cout << descr << id << endl;
}
*/
class ReadArgs
{
  private:
    std::map<std::string,std::string> val;
    std::map<std::string,std::string> descr;
  public:

/*
  * options[] (IN)= contains a list of all possible 1-char options.
  * descriptions[] (IN) = the corresponding list of descriptions.
  *                 WARNING: the descriptions of the options requiring a value
  *                (which is read from the subsequent argument) MUST end with '='
  *                otherwise is treated as a "flag" option, without a value to set to.
  * numoptions (IN)= must contain the no. of admitted options
  * numGivenOptions (OUT)= return the no. of option that have been actually set in the command line. 
  * argc, argv (IN) must be passed from the main program.
  * logfile = if specified will be written into with the pairs "<description of option><?>" each per line,
  *           with <?> that is ": enabled" or "= value" for a non-valued or valued option, respectively.
  *           If the description is empty ("" or " " or "=") then it is not written into logfile.         
  */   
  ReadArgs(std::string options[], std::string descriptions[], int numoptions,
                       int & numGivenOptions, int argc, char *argv[], std::ostream * logfile=0);
   
  
  ~ReadArgs() {}
  

  template <class T>
  inline bool GetParameter(std::string const & option, std::string & description, T & value) {
    if (this->descr.find(option) != descr.end()) description=this->descr[option];
    if (this->val.find(option) != val.end()) { 
      std::istringstream iss(this->val[option]);
      iss >> value;
//      value=ReadArgsS::convert<T>(this->val[option]);
      return true;
    }
    return false;

  }

  template <class T>
  inline bool GetParameter(std::string const & option, T & value) {
    if (this->val.find(option) != val.end()) { 
      std::istringstream iss(this->val[option]);
      iss >> value;
      return true;
    }
    return false;
  }
  
  /*
   * return in value[] the set of values (separated by spaces) given for the <option>
   * and the number of these values. Return -1 if the option wasn't given at all.
   * <maxitems> is the maximum no. of values accepted 
   */ 
  template <class T>
  inline int GetParameter(std::string const & option, std::string & description, T value[], int maxitems) {
    if (this->descr.find(option) != descr.end()) description=this->descr[option];
    if (this->val.find(option) != val.end()) {
      int i;
      std::istringstream iss(this->val[option]);
      for (i=0; i<maxitems && !iss.eof(); i++)  iss >> value[i];
      return (i);
    }
    return -1;

  }

  inline bool isset(std::string const & option) const {
     return (this->val.count(option) != 0);
  }
  
  // return the whole commandline as a string
  inline std::string GetCommandLine(int argc, char *argv[]) const {
    std::string commandline="";
    std::string ar;
    for (int i=0; i<argc; i++) {ar=argv[i]; commandline += ar + " ";}
    return commandline;
  }

//  static int string2array(string const & s, string v[], int maxitems);
  
};

#endif