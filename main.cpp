#include <cstdio>
#include <fstream>

#include "tests.h"
#include "Runs.h"

int main()
{
    //RunAllTests();
    
    //NPR_SimpleG1_WithSubtractions_16cubed();
    
    //NPR_SimpleG1_WithSubtractions_32cubed_Test(SchemeGammaMu, POSITIVE_PARITY);
    //NPR_SimpleG1_WithSubtractions_32cubed_Test(SchemeGammaMu, NEGATIVE_PARITY);
    //NPR_SimpleG1_WithSubtractions_32cubed_Test(SchemeQslash, POSITIVE_PARITY);
    //NPR_SimpleG1_WithSubtractions_32cubed_Test(SchemeQslash, NEGATIVE_PARITY);

    //NPR_SimpleG1_WithSubtractions_24cubed_high(SchemeQslash, POSITIVE_PARITY);
    //NPR_SimpleG1_WithSubtractions_24cubed_high(SchemeQslash, NEGATIVE_PARITY);
    //NPR_SimpleG1_WithSubtractions_32cubed_low_TwoMoms10hits(SchemeGammaMu);

    //printf("ONE MOM FOR G1:\n");
    //NPR_SimpleG1_WithSubtractions_24cubed_low(SchemeGammaMu, NEGATIVE_PARITY);
    //printf("TWO MOMS FOR G1:\n");
    //NPR_SimpleG1_WithSubtractions_24cubed_low_TwoMomsForG1(SchemeGammaMu, NEGATIVE_PARITY);

    //DoStepScaling_24_32(SchemeGammaMu, SchemeGammaMu, NEGATIVE_PARITY);
    DoStepScaling_24_32(SchemeQslash, SchemeQslash, NEGATIVE_PARITY);

    return 0;
}
