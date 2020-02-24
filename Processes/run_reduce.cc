void run_reduce(const char* infile, const char* outfile) {
  gSystem->Load("/home/storage/Users/chenx/atlas/Htautau/delphes/libDelphes.so");
  gSystem->Load("/home/storage/Users/chenx/atlas/Htautau/delphes/reduce_VVV_cc.so");
  reduce(infile,outfile);
}
