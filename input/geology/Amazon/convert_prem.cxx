#include <fstream>
#include <string>

int main(int argc, char *argv[])
{
  if(argc<3)
    abort();

  std::ifstream infile(argv[1]);
  std::string s;
  getline(infile,s);
  getline(infile,s);
  getline(infile,s);

  infile >> depth
         >> lambda
         >> mu;

  std::vector<tuple<double,double,double>> v;
  
  while(infile)
    {
      v.push_back(make_tuple(depth,lambda,mu));
      infile >> depth
             >> lambda
             >> mu;
    }

  const int N(100);
  double min(3.000001), max(2.8909999e+03), step((max-min)/N);
  std::ofstream outfile(argv[2]);
  int i=0;
  for(double d=min;d<max;d+=step)
    {
      while(v[i+1]<d)
        {
          ++i;
        }
      double delta=get<0>(v[i+1]) - get<0>(v[i]);
      

    }
}
