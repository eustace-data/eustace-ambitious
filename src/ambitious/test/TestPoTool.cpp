#include <cppunit/extensions/HelperMacros.h>

#include <ambitious/ostreamer/ostreamer.hpp>
#include <ambitious/po_tool/po_tool.hpp>


#include <string>

namespace po = boost::program_options;


class MockupArgs {
protected:
  int argc_;
  std::vector< char* > argv_;

public:
  MockupArgs() : argc_(0), argv_(0) {
    argv_.push_back(NULL);
  }
  MockupArgs(std::string input) : argc_(0), argv_(0) {
    assign(input);
  }

  int argc() { return argc_; }
  char** argv() { return argv_.data(); }

public:
  MockupArgs& operator=(std::string input) {
    return assign(input);
  }
  MockupArgs& assign(std::string input) {
    clear();

    argv_.resize(0); // Remove the argv_[0] NULL pointer
    argc_ = 0;
    
    size_t pos_end = 0;
    while (pos_end < input.size()) {
      size_t pos_start = input.find_first_not_of(" ", pos_end);
      if (pos_start >= input.size()) {
	break;
      }

      // Parse token or quoted string (does not handle escaped " or ' characters inside strings)
      bool is_string_dbl = (input.data()[pos_start] == '"');
      bool is_string_sgl = (input.data()[pos_start] == '\'');
      if (is_string_dbl) {
	pos_start += 1;
	pos_end = input.find_first_of('"', pos_start);
      } else if (is_string_sgl) {
	pos_start += 1;
	pos_end = input.find_first_of("'", pos_start);
      } else {
	pos_end = input.find_first_of(" =", pos_start);
      }

      // Extract token
      std::string token;
      if (pos_start < pos_end) {
	token = input.substr(pos_start, pos_end-pos_start);
      } else {
	token = "";
      }
      argv_.push_back(new char[token.length() + 1]);
      std::strcpy(argv_[argc_], token.c_str());
      ++argc_;

      // Move to first position after a string or =
      if ((is_string_sgl || is_string_dbl) && (pos_end < input.size())) {
	++pos_end;
      }
      if (pos_end < input.size()) {
	if (input.data()[pos_end] == '=') {
	  ++pos_end;
	}
      }
    }

    // Terminate argv_ with a NULL pointer, so that argv_[argc_] == NULL
    argv_.push_back(NULL);

    return *this;
  }


protected:
  void clear() {
    for (uint64_t i=0; i < argv_.size(); ++i) {
      if (argv_[i] != NULL) {
	delete[] argv_[i];
	argv_[i] = NULL;
      }
    }
    argv_.resize(1);
    argv_[0] = NULL;
    argc_ = 0;
  }
public:
  ~MockupArgs() {
    clear();
  }
};


class Commandline : public CommandlineBase {
protected:
  int param_int_;
  std::string param_string_;
  std::vector<std::string> composing_string_;
  std::vector<std::string> positional_string_;

public:
  Commandline(int argc, char* argv[], const char * config_file)
    : CommandlineBase(argc, argv, config_file),
      param_int_(10),
      param_string_(""),
      composing_string_(0),
      positional_string_(0)
  {

    // Declare a group of options that will be allowed only on command
    // line
    add_generic();
    
    // Declare a group of options that will be allowed both on command
    // line and in config file
    add_config()
      ("int,I",
       po::value<int>(&param_int_)->default_value(param_int_), 
       "Integer")
      ("string,s",
       po::value< std::string >(&param_string_)->default_value(param_string_), 
       "String")
      ("composing-string", 
       po::value< std::vector<std::string> >(&composing_string_)->composing(), 
       "Composing string")
      ;
    
    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    add_hidden()
      ("positional-string",
       po::value< std::vector<std::string> >(&positional_string_),
       "Positional string")
      ;

    add_positional("positional-string", 2);
  }

  const int & param_int() { return param_int_; }
  const std::string & param_string() { return param_string_; }
  const std::vector<std::string> & composing_string() { return composing_string_; }
  const std::vector<std::string> & positional_string() { return positional_string_; }

};


class TestPoTool  : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE( TestPoTool );
  CPPUNIT_TEST( testPositionalCount );
  CPPUNIT_TEST_SUITE_END();

public:
  
  void testPositionalCount()
  {
    MockupArgs args(std::string("program string1  --composing-string \"string3\" string2"));
    Commandline cmdline(args.argc(), args.argv(), "");

    CPPUNIT_ASSERT_EQUAL_MESSAGE("handle_base_options", 0, cmdline.handle_base_options());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("positional_string.size",
				 size_t(2), cmdline.positional_string().size());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("composing_string.size",
				 size_t(1), cmdline.composing_string().size());
    
  }


};

CPPUNIT_TEST_SUITE_REGISTRATION( TestPoTool );
