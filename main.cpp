
#include <clues/CETestGraphGenerator.h>
#include "Makros.h"


int main( int argc, char **argv  ) {
    std::ios_base::sync_with_stdio(0);
    std::cin.tie(NULL);


    CETestGraphGenerator gen;
    gen.generateAllInstances4(); // final test cases


    return 0;
}

