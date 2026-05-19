#ifndef DDVCSTOOLS_H
#define DDVCSTOOLS_H

namespace DDVCSTools {

    enum class detType{
        DET_HTCC = 15,
        DET_FTOF = 12,
            panel1a = 1,
            panel1b = 2,
            panel2 = 3,
        ECAL_TYPE = 7,
            PCal_Layer = 1,
            ECin_Layer = 4,
            ECout_Layer = 7
    };

    static double Eb;
    static const double Mprot = 0.9383;
} // DDVCSTools

#endif //DDVCSTOOLS_H
