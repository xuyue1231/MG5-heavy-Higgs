#!/bin/bash
date
pwd
hostname
rundir=/home/disk-13/xuyue
rm -fr $rundir/tmp_2lep_360
mkdir -p $rundir/tmp_2lep_360
pushd $rundir/tmp_2lep_360
cp /home/storage/Users/chenx/atlas/Htautau/MG5_aMC_v2_5_4/plot_heavy/MG5/LHC_HH_2lep.tgz .
tar xf LHC_HH_2lep.tgz
pushd $rundir/tmp_2lep_360/LHC_HH_2lep
cat Cards/run_card.dat | sed 's/  0   = iseed/  360 = iseed/' > tmp.txt
mv tmp.txt Cards/run_card.dat
rm -f tmp.txt
cat Cards/param_card.dat | sed 's/1 1.000000e+00 # rhoH/1 0.05 # rhoH/'|sed 's/254 5.000000e+02 # MHH/254 300 # MHH/'|sed 's/3 1.300000e+00 # fW/3 160 # fW/'|sed 's/4 1.400000e+00 # fWW/4 960 # fWW/' > tmp.txt
mv tmp.txt Cards/param_card.dat
rm -f tmp.txt
./bin/generate_events -f run1
popd
mv $rundir/tmp_2lep_360/LHC_HH_2lep/Events/run1/unweighted_events.lhe.gz .
gunzip unweighted_events.lhe.gz
mv LHC_HH_2lep/Events/run1/run1_tag_1_banner.txt /home/storage/Users/xuyue/HeavyHiggs/plot_heavy/rootfiles/300GeV/ntuple_2lep/S2_red_360.txt
rm -fr LHC_HH_2lep*
export LD_LIBRARY_PATH=/home/storage/Users/chenx/atlas/Htautau/pythia8219/lib:/usr/local/root/lib:$LD_LIBRARY_PATH
/home/storage/Users/chenx/atlas/Htautau/pythia8219/examples/main_EWZjj /home/storage/Users/chenx/atlas/Htautau/pythia8219/examples/main_EWZjj.cmnd unweighted_events.lhe S2_360.hepmc 360
rm -f S2_360.root
/home/storage/Users/chenx/atlas/Htautau/delphes/DelphesHepMC /home/storage/Users/chenx/atlas/Htautau/delphes/cards/delphes_card_ATLAS_HLVVV.tcl S2_360.root S2_360.hepmc
rm -f S2_360.hepmc
root -l -b -q /home/storage/Users/xuyue/HeavyHiggs/plot_heavy/slurm/run_reduce.cc"(\"S2_360.root\",\"S2_red_360.root\")"
rm -f S2_360.root
mv S2_red_360.root /home/storage/Users/xuyue/HeavyHiggs/plot_heavy/rootfiles/300GeV/ntuple_2lep
popd
rm -fr $rundir/tmp_2lep_360
date
