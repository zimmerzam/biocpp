#ifndef DEFINITION_H
#define DEFINITION_H

#include <libconfig.h++>

namespace BioCpp{

class definition{
  public:
    virtual void importSetting(libconfig::Setting& setting) =0;
};

}
#endif
