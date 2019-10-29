void out1() {
gROOT->ProcessLine(".L fit_com_1.cc++");
gSystem->Load("fit_com_1_cc.so");
plot();
}
