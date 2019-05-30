// Based on example code Copyright Vladimir Prus 2002-2004.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

/* Shows how to use both command line and config file. */

#include "ambitious/ostreamer/ostreamer.hpp"

#include <string.h>
#include <iostream>
#include "ambitious/po_tool/po_tool.hpp"

namespace po = boost::program_options;


class Commandline : public CommandlineBase {
protected:
  int opt_ = 10;
  uint64_t seed_ = 0;
  std::vector<std::string> include_path_;
  std::vector<std::string> input_file_;
  std::vector<std::string> output_file_;
  std::string F_;
  bool B_;

public:
  Commandline(int ac, char* av[], const char * config_file)
    : CommandlineBase(ac, av, config_file) {

    // Declare a group of options that will be allowed only on command
    // line
    add_generic();
    
    // Declare a group of options that will be allowed both on command
    // line and in config file
    add_config()
      ("optimization,O",
       po::value<int>(&opt_)->default_value(opt_), 
       "Optimization level")
      ("seed,S",
       po::value<uint64_t>(&seed_)->default_value(seed_), 
       "Random seed")
      ("include-path,I", 
       po::value< std::vector<std::string> >(&include_path_)->composing(), 
       "Include path")
      ("F,F", 
       po::value< std::string >(&F_)->default_value(""), 
       "string test")
      ("B,B", 
       po::value< bool >(&B_)->default_value(false), 
       "bool test")
      ;
    
    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    add_hidden()
      ("input-file", po::value< std::vector<std::string> >(&input_file_), "Input file")
      ("output-file", po::value< std::vector<std::string> >(&output_file_), "Output file")
      ;

    add_positional("input-file", 2);
    add_positional("output-file", 2);
  }

  int opt() { return opt_; }
  uint64_t seed() { return seed_; }
  const std::vector<std::string> & include_path() { return include_path_; }
  const std::vector<std::string> & input_file() { return input_file_; }
  const std::vector<std::string> & output_file() { return output_file_; }

};












int main(int ac, char* av[])
{
  try {
    Commandline cmdline(ac, av, "config.cfg");
    if (cmdline.handle_base_options()) {
      return 0;
    }
        
    if (cmdline.vm().count("include-path")) {
      std::cout << "Include paths are: " << cmdline.include_path() << std::endl;
    }
    
    if (cmdline.vm().count("input-file")) {
      std::cout << "Input files are: " << cmdline.input_file() << std::endl;
    }
    
    if (cmdline.vm().count("output-file")) {
      std::cout << "Output files are: " << cmdline.output_file() << std::endl;
    }
    
    std::cout << "Optimization level is " << cmdline.opt() << std::endl;
    std::cout << "Seed is " << cmdline.seed() << std::endl;

  }
  catch(std::exception& e)
    {
      std::cout << e.what() << std::endl;
      return 1;
    }    
  return 0;
}
