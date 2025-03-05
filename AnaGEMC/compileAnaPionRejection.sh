g++ AnaPionRejection.cc -o AnaPionRejection.exe `root-config --cflags --libs` -I$CLAS12AnaTools/include -L$CLAS12AnaTools/lib -lclas12AnaTools -I$HIPODIR/hipo4 $HIPODIR/lib/libhipo4.a -lSpectrum
