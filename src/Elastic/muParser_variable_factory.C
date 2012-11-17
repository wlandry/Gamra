#include <list>

double* muParser_variable_factory(const char *, void *)
{
  static std::list<double> variables;
  variables.push_back(0);
  return &variables.back();
}
