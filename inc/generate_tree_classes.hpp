#ifndef H_GENERATE_TREE_CLASSES
#define H_GENERATE_TREE_CLASSES

#include <vector>
#include <string>

class Variable{
public:
  Variable():
    type_(""),
    name_(""),
    out_name_(""){
  }

  Variable(const std::string &type,
           const std::string &base_type,
           const std::string &name,
           const std::string &out_name):
    type_(type),
    base_type_(base_type),
    name_(name),
    out_name_(out_name){
  }

  bool operator<(const Variable& var) const{
    //return type_<var.type_ || (type_==var.type_ && name_<var.name_);
    return type_<var.type_ || (type_==var.type_ && out_name_<var.out_name_);
  }

  std::string type_, base_type_, name_, out_name_;
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
