#ifndef PO_TOOL_HPP
#define PO_TOOL_HPP

#include "ambitious/versions/versions.hpp"
VERSIONS_ADD(po_tool_hpp, "$Revision: 1258 $")

// Based on example code Copyright Vladimir Prus 2002-2004.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

/* Shows how to use both command line and config file. */

// Might stop compiler errors with boost:
#include <utility>

#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <boost/filesystem.hpp>
#include <boost/any.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;



class CommandlineBase {
protected:
  bool is_processed_ = false;
  // Config file name
  std::string config_file_;
  // Config file output name
  std::string config_file_out_;
  // Base option containters used for setup
  boost::program_options::options_description generic_;
  boost::program_options::options_description config_;
  boost::program_options::options_description hidden_;
  boost::program_options::positional_options_description positional_;
  // Aggregated options used for internal processing 
  boost::program_options::options_description cmdline_options_;
  boost::program_options::options_description config_file_options_;
  boost::program_options::options_description visible_;
  // Storage
  boost::program_options::variables_map vm_;
  int ac_;
  char** av_;
  // Store positional argument information separately, to simplify usage printing
  typedef std::list< std::pair<std::string, int> > positional_list_class;
  positional_list_class positional_list_;

  std::string program_name_;
  std::string program_version_;

  int recommended_return_value_;
  
  void PrintVariableMap(std::ostream &output_stream);
  void process_config_file(std::istream &input_stream);

  void add_base_generic();
  
public:
  CommandlineBase(int ac, char* av[], const std::string& config_file);
  CommandlineBase(int ac, char* av[], const std::string& config_file,
                  const std::string& program_name,
                  const std::string& program_version);
  int ret() const {
    return recommended_return_value_;
  }
  
  void process();
  
  boost::program_options::variables_map& vm() {
    if (!is_processed_) process();
    return vm_;
  }

  std::string config_file() {
    if (!is_processed_) process();
    return config_file_;
  }

  int handle_help();
  void set_name_version(const std::string& program_name,
                        const std::string& program_version);
  int handle_version();
  int handle_printconfig();

  /** Handle help&config printing, etc
   * Returns non-zero if recommending termination of the program.
   *
   * Use ret to access the recommended the exit code.
   */
  int handle_base_options();

  /** Handle help&config printing, catching exceptions
   *
   * Returns non-zero if recommending termination of the program.
   *
   * If returning -1, the termination is intentional.
   * If returning a positive number, an error occurred.
   *
   * Use ret to access the recommended the exit code.
   */
  int handle_options();



  /** Add a generic option.
   *
   * Only allowed before process() is called.
   */ 
  boost::program_options::options_description_easy_init add_generic();
  /** Add a config option.
   *
   * Only allowed before process() is called.
   */ 
  boost::program_options::options_description_easy_init add_config();
  /** Add a hidden option.
   *
   * Only allowed before process() is called.
   */ 
  boost::program_options::options_description_easy_init add_hidden();

  /** Add a positional option.
   *
   * Specifies that up to 'max_count' next positional options should
   * be given the 'name'. The value of '-1' means 'unlimited'. No
   * calls to 'add' can be made after call with 'max_value' equal to
   * '-1'.
   *
   * Only allowed before process() is called.
   */ 
  void add_positional(const char * name, int max_count);
  
  /** Returns the maximum number of positional options that can be
   *  present. Can return numeric_limits<unsigned>::max() to indicate
   *  an unlimited number.
   */
  unsigned positional_max_total_count() const;
  /** Returns the name that should be associated with positional
   *  options at 'position'. Precondition: position <
   *  max_total_count()
   */
  const std::string & positional_name_for_position(unsigned position) const;


    unsigned max_total_count() const;


    const std::string & name_for_position(unsigned position) const;


	

  
};


#endif
