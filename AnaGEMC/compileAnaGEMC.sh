g++ AnaGEMC.cc -o AnaGEMC.exe `root-config --libs --cflags` -I$CLAS12AnaTools/include -L$CLAS12AnaTools/lib -lclas12AnaTools -I$HIPODIR/hipo4 $HIPODIR/lib/libhipo4.a -lSpectrum
