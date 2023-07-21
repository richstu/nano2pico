#ifndef H_GENERATE_TREE_CLASSES
#define H_GENERATE_TREE_CLASSES

#include <vector>
#include <string>
#include <regex>
#include <iomanip>
#include <sstream>

class Variable{
public:
  Variable():
    type_(""),
    name_(""),
    version_(0.){
  }

  Variable(const std::string &type,
           const std::string &base_type,
           const std::string &name,
           const float &version):
    type_(type),
    base_type_(base_type),
    name_(name),
    version_(version){
  }

  std::string get_out_name() const{
    std::string out_name = name_;
    if (version_>0.01) {
      std::stringstream version_stream;
      version_stream << std::fixed << std::setprecision(1) << version_;
      out_name += "_"+version_stream.str();
    }
    out_name = std::regex_replace(out_name, std::regex("\\."), "p");
    return out_name;
  }

  bool operator<(const Variable& var) const{
    //return type_<var.type_ || (type_==var.type_ && name_<var.name_);
    return type_<var.type_ || (type_==var.type_ && get_out_name()<var.get_out_name());
  }

  std::string type_, base_type_, name_;
  float version_;
};

bool Contains(const std::string &text, const std::string &pattern);

void WriteHeader(const std::vector<Variable> &vars, 
                 const std::string name);
void WriteSource(const std::vector<Variable> &vars, 
                 const std::string name);

void WriteNanoHeader(const std::vector<Variable> &vars);
void WriteNanoSource(const std::vector<Variable> &vars);

int GetArrayLength(const std::string var_name);

#endif
