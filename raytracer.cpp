#include "picture.h"
#include "raytracing/raytracing.h"

int main(int argc, const char** argv) 
{
  std::unordered_map<std::string, std::string> cmdLineParams;

  for(int i=0; i<argc; i++)
  {
    std::string key(argv[i]);

    if(key.size() > 0 && key[0]=='-')
    {
      if(i != argc-1) // not last argument
      {
        cmdLineParams[key] = argv[i+1];
        i++;
      }
      else
        cmdLineParams[key] = "";
    }
  }

std::string outFilePath = "out.ppm";
  if(cmdLineParams.find("-out") != cmdLineParams.end())
    outFilePath = cmdLineParams["-out"];
    
  try {
      picture::Picture my_pic;
      raytracing::render(outFilePath.c_str(), my_pic.spheres, my_pic.lights, my_pic.cubes);
      std::cout << "result is saved in build directory in file " << outFilePath.c_str() << std::endl;
  }
  catch (const std::exception& er){
      std::cout << er.what();
  }
  return 0;
}

