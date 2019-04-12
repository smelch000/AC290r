//#include "/home/miocchi/Sources/include/readargs.h"
#include "./readargs.h"
/*
  * options[] (IN)= contains a list of all possible 1-char options.
  * descriptions[] (IN) = the corresponding list of descriptions.
  *                 WARNING: the descriptions of the options requiring a value
  *                (which is read from the subsequent argument) MUST end with '='
  *                otherwise is treated as a "flag" option, without a value to set to.
  * 
  * numoptions (IN)= must contain the no. of admitted options
  * numGivenOptions (OUT)= return the no. of option that have been actually set in the command line. 
  * argc, argv (IN) must be passed from the main program.
  * logfile = if specified will be written into with the pairs "<description of option><?>" each per line,
  *           with <?> that is ": enabled" or "= value" for a non-valued or valued option, respectively.
  *           If the description is empty ("" or " " or "=") then it is not written into logfile.         
  */ 

ReadArgs::ReadArgs(std::string options[], std::string descriptions[], int numoptions,
                      int & numGivenOptions, int argc, char *argv[], std::ostream * logfile) {
  numGivenOptions=0;
  std::string option;
  int i,j;
  std::map<std::string, bool> valued;
  //store all options' descriptions into the associative array descr
  for (i=0; i<numoptions; i++) {
    if (descriptions[i]=="" || descriptions[i] == "=" || descriptions[i]==" ") {
      descr[options[i]]=" "; //no descr. will be written in logfile
      valued[options[i]]=(descriptions[i] == "=");
    } else {
      descr[options[i]]= descriptions[i];
      valued[options[i]]= (descriptions[i][descriptions[i].length()-1] == '=');
    }
  }
  i=1;
  while (i<argc) {
    if (argv[i][0]=='-') {
      // the argument is an option or a list of options (i.e., starting with '-')
      // consider every options in the list, but the last.
      for (j=1; j< strlen(argv[i])-1; j++) {
        option=argv[i][j];
        if (descr.find(option) != descr.end()) {
          //the option has been described
          if (!valued[option]) {
          //the option does not require a value
            numGivenOptions++;
            //the option within the list (but the last) has no value  
            val[option]="";
            if (logfile && descr[option] != " ")
              (*logfile) << descr[option] << ": enabled\n";
          } else {
            std::cerr << "WARNING: option -"<<option<<" was not set.\n";
          }
        } else {
          std::cerr << "WARNING: unknown option -" << option << std::endl; 
        }
      }
      //consider the last option in the list
      option=argv[i][j];
      if (descr.find(option) != descr.end()) {
        //the option has been described
        if (valued[option]) {
          //the option requires a value
          //take the subsequent argument as the value of the last option in the list 
          if (i+1<argc) {
            numGivenOptions++;
            val[option]=argv[i+1];
            //read a possible multi-valued option value
            while (i+2 <argc && argv[i+2][0] != '-') {
              i++;
              val[option]= val[option]+" "+argv[i+1];
            }
            if (logfile && descr[option] != " ")
              (*logfile) << descr[option] << " " << val[option] << std::endl;
            i++; //thus the value is jumped
          } else {
            std::cerr << "WARNING: option -"<<option<<" was not set.\n";
          }
        } else {
          numGivenOptions++;
          val[option]="";
          if (logfile && descr[option] != " ")
            (*logfile) << descr[option] << ": enabled\n";
        }
      } else {
        std::cerr << "Warning: unknown option -" << option << std::endl; 
      }
    }
    i++;
  }
  
}


  
