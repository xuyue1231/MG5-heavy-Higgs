#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  ParticlePropagator

  ChargedHadronMomentumSmearing
  ElectronMomentumSmearing
  MuonMomentumSmearing

  TrackMerger
  IdentificationMap
  ECal

  TreeWriter
}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
  set InputArray Delphes/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius 1.25
  # half-length of the magnetic field coverage, in m
  set HalfLength 1.96

  # magnetic field
  set Bz 1.5
}

########################################
# Momentum resolution for charged tracks
########################################

# https://www.sciencedirect.com/science/article/pii/S0168900200005131

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

  # resolution formula for charged hadrons
  set ResolutionFormula { sqrt(0.003^2 + pt^2*0.001^2) }


}

###################################
# Momentum resolution for electrons
###################################

module MomentumSmearing ElectronMomentumSmearing {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons

  # set ResolutionFormula {resolution formula as a function of eta and energy}

   # resolution formula for charged hadrons
  set ResolutionFormula { sqrt(0.003^2 + pt^2*0.001^2) }
}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray ParticlePropagator/muons
  set OutputArray muons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

   # resolution formula for charged hadrons
  set ResolutionFormula { sqrt(0.003^2 + pt^2*0.001^2) }

}

##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronMomentumSmearing/electrons
  add InputArray MuonMomentumSmearing/muons
  set OutputArray tracks
}

####################################
# Charged hadron PID
####################################

module IdentificationMap IdentificationMap {
  set InputArray TrackMerger/tracks
  set OutputArray tracks

  # {PID in} {PID out} {formula}
  # make sure "PID in" and "PID out" have the same charge (e.g {-13} {211} or {-321} {211})
  # {211} {-13} is equivalent to {-211} {13} (and needs to be written once only...)

  # --- pions ---

  add EfficiencyFormula {211} {211} {      (eta <= -1.317) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt < 0.1) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.1 && pt < 0.7) * (0.99) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.7 && pt < 1.0) * (1.00-0.038) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 1.0) * (1.00-0.013) +
					   (eta > -0.745  && eta <= 1.240) *       (pt < 0.1)* (0.00) +
					   (eta > -0.745  && eta <= 1.240) *       (pt >= 0.1 && pt < 0.7)* (0.995) +
					   (eta > -0.745  && eta <= 1.240) *       (pt >= 0.7 && pt < 1.0)* (0.995-0.038) +
					   (eta > -0.745  && eta <= 1.240) *       (pt >= 1.0)* (0.995-0.013) +
					   (eta > 1.240  && eta <= 1.901) *       (pt < 0.1)* (0.00) +
					   (eta > 1.240  && eta <= 1.901) *       (pt >= 0.1 && pt < 0.7)* (0.99) +
					   (eta > 1.240  && eta <= 1.901) *       (pt >= 0.7 && pt < 1.0)* (0.99-0.038) +
					   (eta > 1.240  && eta <= 1.901) *       (pt >= 1.0)* (0.99-0.013) +
					   (eta > 1.901) * (0.00)}

  add EfficiencyFormula {211} {321} {      (eta <= -1.317) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt < 0.1) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.1 && pt < 0.7) * (0.01) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.7) * (0.00) +
					   (eta > -0.745  && eta <= 1.240) *       (pt < 0.1)* (0.00) +
					   (eta > -0.745  && eta <= 1.240) *       (pt >= 0.1)* (0.005) +
					   (eta > 1.240  && eta <= 1.901) *       (pt < 0.1)* (0.00) +
					   (eta > 1.240  && eta <= 1.901) *       (pt >= 0.1)* (0.01) +
					   (eta > 1.901) * (0.00)}

  add EfficiencyFormula {211} {-13} {      (eta <= -1.317) * (0.00) +
                                           (eta > -1.317  && eta <= 1.901) *       (pt < 0.1) * (0.00) +
                                           (eta > -1.317  && eta <= 1.901) *       (pt >= 0.1 && pt < 0.7) * (0.00) +
                                           (eta > -1.317  && eta <= 1.901) *       (pt >= 0.7 && pt < 1.0) * (0.038) +
                                           (eta > -1.317  && eta <= 1.901) *       (pt >= 1.0) * (0.013) +
					   (eta > 1.901) * (0.00)}

 # --- kaons ---

  add EfficiencyFormula {321} {321} {      (eta <= -1.317) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt < 0.1) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.1 && pt < 0.7) * (0.85) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.7 && pt < 1.0) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 1.0) * (0.00) +
					   (eta > -0.745  && eta <= 1.240) *       (pt < 0.1)* (0.00) +
					   (eta > -0.745  && eta <= 1.240) *       (pt >= 0.1 && pt < 0.7)* (0.99) +
					   (eta > -0.745  && eta <= 1.240) *       (pt >= 0.7 && pt < 1.0)* (0.99-0.038) +
					   (eta > -0.745  && eta <= 1.240) *       (pt >= 1.0)* (0.99-0.013) +
					   (eta > 1.240  && eta <= 1.901) *       (pt < 0.1)* (0.00) +
					   (eta > 1.240  && eta <= 1.901) *       (pt >= 0.1 && pt < 0.7)* (0.96) +
					   (eta > 1.240  && eta <= 1.901) *       (pt >= 0.7 && pt < 1.0)* (0.96-0.038) +
					   (eta > 1.240  && eta <= 1.901) *       (pt >= 1.0)* (0.96-0.013) +
					   (eta > 1.901) * (0.00)}

  add EfficiencyFormula {321} {211} {      (eta <= -1.317) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt < 0.1) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.1 && pt < 0.7) * (0.15) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.7 && pt < 1.0) * (1.00-0.038) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 1.0) * (1.00-0.013) +
					   (eta > -0.745  && eta <= 1.240) *       (pt < 0.1)* (0.00) +
					   (eta > -0.745  && eta <= 1.240) *       (pt >= 0.1)* (0.01) +
					   (eta > 1.240  && eta <= 1.901) *       (pt < 0.1)* (0.00) +
					   (eta > 1.240  && eta <= 1.901) *       (pt >= 0.1)* (0.04) +
					   (eta > 1.901) * (0.00)}

  add EfficiencyFormula {321} {-13} {      (eta <= -1.317) * (0.00) +
                                           (eta > -1.317  && eta <= 1.901) *       (pt < 0.1) * (0.00) +
                                           (eta > -1.317  && eta <= 1.901) *       (pt >= 0.1 && pt < 0.7) * (0.00) +
                                           (eta > -1.317  && eta <= 1.901) *       (pt >= 0.7 && pt < 1.0) * (0.038) +
                                           (eta > -1.317  && eta <= 1.901) *       (pt >= 1.0) * (0.013) +
					   (eta > 1.901) * (0.00)} 

 # --- protons ---

  add EfficiencyFormula {2212} {2212} {    (eta <= -1.317) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt < 0.1) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.1 && pt < 0.7) * (0.85) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.7 && pt < 1.0) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 1.0) * (0.00) +
					   (eta > -0.745  && eta <= 1.901) *       (pt < 0.1)* (0.00) +
					   (eta > -0.745  && eta <= 1.901) *       (pt >= 0.1 && pt < 0.7)* (0.99) +
					   (eta > -0.745  && eta <= 1.901) *       (pt >= 0.7 && pt < 1.0)* (0.99-0.038) +
					   (eta > -0.745  && eta <= 1.901) *       (pt >= 1.0)* (0.99-0.013) +
					   (eta > 1.901) * (0.00)}

  add EfficiencyFormula {2212} {211} {    (eta <= -1.317) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt < 0.1) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.1 && pt < 0.7) * (0.15) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.7 && pt < 1.0) * (1.00-0.038) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 1.0) * (1.00-0.013) +
					   (eta > -0.745  && eta <= 1.901) *       (pt < 0.1)* (0.00) +
					   (eta > -0.745  && eta <= 1.901) *       (pt >= 0.1)* (0.01) +
					   (eta > 1.901) * (0.00)}

  add EfficiencyFormula {2212} {-13} {     (eta <= -1.317) * (0.00) +
                                           (eta > -1.317  && eta <= 1.901) *       (pt < 0.1) * (0.00) +
                                           (eta > -1.317  && eta <= 1.901) *       (pt >= 0.1 && pt < 0.7) * (0.00) +
                                           (eta > -1.317  && eta <= 1.901) *       (pt >= 0.7 && pt < 1.0) * (0.038) +
                                           (eta > -1.317  && eta <= 1.901) *       (pt >= 1.0) * (0.013) +
					   (eta > 1.901) * (0.00)}

 # --- muons ---

  add EfficiencyFormula {-13} {-13} {     (eta <= -1.317) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt < 0.1) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.1 && pt < 0.7) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.7 && pt < 1.0) * (0.76) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 1.0) * (0.89) +
					   (eta > -0.745  && eta <= 1.901) *       (pt < 0.1)* (0.00) +
					   (eta > -0.745  && eta <= 1.901) *       (pt >= 0.1 && pt < 0.7)* (0.60) +
					   (eta > -0.745  && eta <= 1.901) *       (pt >= 0.7 && pt < 1.0)* (0.76) +
					   (eta > -0.745  && eta <= 1.901) *       (pt >= 1.0)* (0.89) +
					   (eta > 1.901) * (0.00)}

  add EfficiencyFormula {-13} {211} {     (eta <= -1.317) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt < 0.1) * (0.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.1 && pt < 0.7) * (1.00) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 0.7 && pt < 1.0) * (0.24) +
                                           (eta > -1.317  && eta <= -0.745) *       (pt >= 1.0) * (0.11) +
					   (eta > -0.745  && eta <= 1.901) *       (pt < 0.1)* (0.00) +
					   (eta > -0.745  && eta <= 1.901) *       (pt >= 0.1 && pt < 0.7)* (0.40) +
					   (eta > -0.745  && eta <= 1.901) *       (pt >= 0.7 && pt < 1.0)* (0.24) +
					   (eta > -0.745  && eta <= 1.901) *       (pt >= 1.0)* (0.11) +
					   (eta > 1.901) * (0.00)}

# --- electrons ---

  add EfficiencyFormula {-11} {-11} {      (eta <= -1.317)                                * (0.00) +
                                           (eta > -1.317  && eta <= 1.901) *      (pt < 0.1)* (0.00) +
                                           (eta > -1.317  && eta <= 1.901) *     (pt >= 0.1)* (0.96) +
                                           (eta > 1.901)                                 * (0.00)}

  add EfficiencyFormula {-11} {211} {     (eta <= -1.317) * (0.00) +
                                           (eta > -1.317  && eta <= 1.901) *       (pt < 0.1) * (0.00) +
                                           (eta > -1.317  && eta <= 1.901) *       (pt >= 0.1) * (0.04) +
					   (eta > 1.901) * (0.00)}

 # efficiency for other charged particles (should be always {0} {0} {formula})

  add EfficiencyFormula {0} {0}     {      (eta <= -1.317)                                * (0.00) +
                                           (eta > -1.317  && eta <= 1.901) *     (pt < 0.1) * (0.00) +
					   (eta > -1.317  && eta <= 1.901) *     (pt > 0.1) * (0.95) +
                                           (eta > 1.901)                                 * (0.00)}

}

#############
#   ECAL
#############

module SimpleCalorimeter ECal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray IdentificationMap/tracks

  set TowerOutputArray ecalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowPhotons

  set IsEcal true 
 
  set EnergyMin 0.02
  set EnergySignificanceMin 0.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # 0.5 degree towers
  set PhiBins {}
  for {set i -360} {$i <= 360} {incr i} {
    add PhiBins [expr {$i * $pi/360.0}]
  }

  # 0.01 unit in eta from eta = -1.317 to eta = 1.901
  for {set i 0} {$i <= 322} {incr i} {
    set eta [expr {-1.317 + $i * 0.01}]
    add EtaPhiBins $eta $PhiBins
  }

  # default energy fractions {abs(PDG code)} {fraction of energy deposited in ECAL}

  add EnergyFraction {0} {0.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {1.0}
  add EnergyFraction {22} {1.0}
  add EnergyFraction {111} {1.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.3}
  add EnergyFraction {3122} {0.3}

  # set ECalResolutionFormula {resolution formula as a function of eta and energy}

  set ResolutionFormula {(eta > -1.317 && eta < 1.901) * sqrt(energy^2*0.012^2 + energy^1.5*0.016^2 + 0.002^2)}
}

##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch Delphes/allParticles Particle GenParticle

  add Branch ECal/eflowTracks EFlowTrack Track
  add Branch ECal/eflowPhotons EFlowPhoton Tower
}
