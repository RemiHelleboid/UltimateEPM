#include <string>
#include <vector>

#include "tclap/CmdLine.h"

using namespace TCLAP;
using namespace std;

// https://sourceforge.net/p/tclap/bugs/30/
int main() {
    CmdLine             cmd("test empty argv");
    std::vector<string> args;
    cmd.parse(args);
}
