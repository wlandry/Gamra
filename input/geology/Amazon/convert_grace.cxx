#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

namespace fs = boost::filesystem;

int main()
{
  fs::directory_iterator end_iter;
  for( fs::directory_iterator dir_iter("GRACEf") ; dir_iter != end_iter ; ++dir_iter)
    {
      fs::ifstream infile(*dir_iter);
      const int nx(56), ny(77);
      double data[nx][ny];
      double x,y;

      std::string uid(dir_iter->path().leaf().stem().string().substr(0,8));
      fs::path outfile_path="Amazon/Amazon_" + uid + ".input";
      fs::copy_file("Amazon_begin",outfile_path,
                    fs::copy_option::overwrite_if_exists);
      fs::ofstream outfile(outfile_path,std::ios_base::app);
      outfile << uid;
      outfile << fs::ifstream("Amazon_middle").rdbuf();

      for(int i=nx-1;i>=0;--i)
        for(int j=ny-1;j>=0;--j)
          {
            infile >> x >> y >> data[i][j];
          }
      for(int i=0;i<nx;++i)
        for(int j=0;j<ny;++j)
          {
            outfile << data[i][j];
            if(i!=nx-1 || j!=ny-1)
              outfile << ",\n";
          }

      outfile << fs::ifstream("Amazon_end").rdbuf();
    }
}
