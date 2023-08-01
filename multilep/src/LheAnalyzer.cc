#include "heavyNeutrino/multilep/interface/LheAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "TLorentzVector.h"

#include <algorithm>
/*
 * Accessing LHE information
 * lheHTIncoming could be used to get the low HT bin for HT-binned samples, e.g. DY
 * Might consider skimming directly here for such samples
 * Also saving the ctau of the heavy neutrino
 * If the storeLheParticles boolean is set, most of the LHE particle information is stored to the tree
 * Also keeping track of LHE taus in the event [to be used in case this is run on a sample where pythia decays all taus leptonically]
 */
LheAnalyzer::LheAnalyzer(const edm::ParameterSet& iConfig, multilep* multilepAnalyzer):
    multilepAnalyzer(multilepAnalyzer)
{};


void LheAnalyzer::beginJob(TTree* outputTree, edm::Service<TFileService>& fs){
    if( multilepAnalyzer->isData() ) return;
    hCounter   = fs->make<TH1D>("hCounter",   "Events counter", 1, 0, 1);
    lheCounter = fs->make<TH1D>("lheCounter", "Lhe weights",    maxNumberOfLheWeights, 0, maxNumberOfLheWeights); //Counter to determine effect of pdf and scale uncertainties on the MC cross section
    psCounter  = fs->make<TH1D>("psCounter",  "Lhe weights",    maxNumberOfPsWeights, 0, maxNumberOfPsWeights);
    tauCounter = fs->make<TH1D>("tauCounter", "Number of taus", 3, 0, 3);

    if (multilepAnalyzer->storeEFTInfo) {
        eftCounter = fs->make<TH1D>("eftCounter", "EFT weights", maxNumberOfEFTWeights, 0, maxNumberOfEFTWeights);
        dynScaleCounter = fs->make<TH1D>("dynScaleCounter", "Dynamic Scale LHE weights", maxNumberOfDynScaleWeights, 0, maxNumberOfDynScaleWeights);
    }

    nTrueInteractions = fs->make<TH1D>("nTrueInteractions", "nTrueInteractions", 100, 0, 100);

    outputTree->Branch("_nTrueInt",      &_nTrueInt,      "_nTrueInt/F");
    outputTree->Branch("_weight",        &_weight,        "_weight/D");
    outputTree->Branch("_lheHTIncoming", &_lheHTIncoming, "_lheHTIncoming/D");
    outputTree->Branch("_ctauHN",        &_ctauHN,        "_ctauHN/D");
    outputTree->Branch("_nLheTau",       &_nTau,          "_nLheTau/i");
    outputTree->Branch("_nLheWeights",   &_nLheWeights,   "_nLheWeights/i");
    outputTree->Branch("_lheWeight",     &_lheWeight,     "_lheWeight[_nLheWeights]/D");
    outputTree->Branch("_nPsWeights",    &_nPsWeights,    "_nPsWeights/i");
    outputTree->Branch("_psWeight",      &_psWeight,      "_psWeight[_nPsWeights]/D");

    if(multilepAnalyzer->storeLheParticles){
        outputTree->Branch("_nLheParticles", &_nLheParticles, "_nLheParticles/i");
        outputTree->Branch("_lheStatus",     &_lheStatus,     "_lheStatus[_nLheParticles]/I");
        outputTree->Branch("_lhePdgId",      &_lhePdgId,      "_lhePdgId[_nLheParticles]/I");
        outputTree->Branch("_lheMother1",    &_lheMother1,    "_lheMother1[_nLheParticles]/I");
        outputTree->Branch("_lheMother2",    &_lheMother2,    "_lheMother2[_nLheParticles]/I");
        outputTree->Branch("_lhePt",         &_lhePt,         "_lhePt[_nLheParticles]/F");
        outputTree->Branch("_lheEta",        &_lheEta,        "_lheEta[_nLheParticles]/F");
        outputTree->Branch("_lhePhi",        &_lhePhi,        "_lhePhi[_nLheParticles]/F");
        outputTree->Branch("_lheE",          &_lheE,          "_lheE[_nLheParticles]/F");
        outputTree->Branch("_lheMass",       &_lheMass,       "_lheMass[_nLheParticles]/F");
    }

    if (multilepAnalyzer->storeEFTInfo) {
        outputTree->Branch("_nDynScaleWeights",   &_nDynScaleWeights,   "_nDynScaleWeights/i");
        outputTree->Branch("_dynScaleWeight",     &_dynScaleWeight,     "_dynScaleWeight[_nLheWeights]/D");
        outputTree->Branch("_nEFTWeights",        &_nEFTWeights,        "_nEFTWeights/i");
        outputTree->Branch("_eftWeight",          &_eftWeight,          "_eftWeight[_nLheWeights]/D");
    }
}

void LheAnalyzer::beginRun(const edm::Run& iRun) {
    // if (! multilepAnalyzer->storeEFTInfo)
    edm::Handle<LHERunInfoProduct> lheRunInfo              = getHandle(iRun, multilepAnalyzer->lheEventInfoToken);

    std::vector<LHERunInfoProduct::Header>::const_iterator it;

    bool localSearchFunction = [](const LHERunInfoProduct::Header& el2) {return el2.tag() == "initrwgt";};
    it = std::find_if(lheRunInfo->headers_begin(), lheRunInfo->headers_end(), localSearchFunction);
    // add conditions here first to find start of relevant weights, we don't care about EFT, just about what IDs are relevant for nominal lheweights -> so 1001 always, but then nondynscal ones

    std::vector<std::string> lheConfiguration = it->lines();

    // find '<weightgroup name="mg_reweighting"' tagged as mg_reweighting. Save all weights there aside from dummy. Can add an explicit check for that

    // loop, find weightgroup:
    // define an nStartEFT (if multilep->storeEFT) is on
    // define a startNominalLHE and StepLHE -> first is just constant offset, second is multiplier in for loop, based on if there are dyn scale objects?
    // Also add an offset pdf, to take into account additional weights for pdfs

    // it->lines() should now start at relevant section
    // in addition, find ids of nominal weights
    
    // also need somewhere here the hyperpoly stuff to initialize
}

void LheAnalyzer::analyze(const edm::Event& iEvent){
    if( multilepAnalyzer->isData() ) return;
    edm::Handle<GenEventInfoProduct> genEventInfo          = getHandle(iEvent, multilepAnalyzer->genEventInfoToken);
    edm::Handle<LHEEventProduct> lheEventInfo              = getHandle(iEvent, multilepAnalyzer->lheEventInfoToken);
    edm::Handle<std::vector<PileupSummaryInfo>> pileUpInfo = getHandle(iEvent, multilepAnalyzer->pileUpToken);

    _nTrueInt = pileUpInfo->begin()->getTrueNumInteractions(); // getTrueNumInteractions is the same for all bunch crossings
    _weight   = genEventInfo->weight();
    hCounter->Fill(0.5, _weight);
    nTrueInteractions->Fill(_nTrueInt, _weight);

    _lheHTIncoming = 0.;
    _ctauHN = 0.;
    _nTau = 0;

    if(!lheEventInfo.isValid()){
        _nLheWeights = 0;
        return;
    }


    // See http://home.thep.lu.se/~leif/LHEF/classLHEF_1_1HEPEUP.html for more detailes
    _nLheParticles = lheEventInfo->hepeup().NUP;
    if(_nLheParticles > nLhe_max){
        throw cms::Exception("maxLheParticles") << "This process has " << _nLheParticles << " lhe particles. Please increase nLhe_max in LheAnalyzer.h.\n";
    }
    for(unsigned i = 0; i < _nLheParticles; ++i){
        _lheStatus[i]          = lheEventInfo->hepeup().ISTUP[i];
        _lhePdgId[i]           = lheEventInfo->hepeup().IDUP[i];
        _lheMother1[i]         = lheEventInfo->hepeup().MOTHUP[i].first-1;
        _lheMother2[i]         = lheEventInfo->hepeup().MOTHUP[i].second-1;
        double px              = lheEventInfo->hepeup().PUP[i][0];
        double py              = lheEventInfo->hepeup().PUP[i][1];
        double pz              = lheEventInfo->hepeup().PUP[i][2];
        _lheE[i]               = lheEventInfo->hepeup().PUP[i][3];
        _lheMass[i]            = lheEventInfo->hepeup().PUP[i][4];
        TLorentzVector vector  = TLorentzVector(px, py, pz, _lheE[i]);
        _lhePt[i]              = vector.Et();
        _lheEta[i]             = vector.Eta();
        _lhePhi[i]             = vector.Phi();

        bool hasIncomingMother = lheEventInfo->hepeup().ISTUP[_lheMother1[i]]==-1 and lheEventInfo->hepeup().ISTUP[_lheMother2[i]]==-1;  // Status -1 means mother is incoming
        bool quarkOrGluon      = (_lhePdgId[i]==21 or (abs(_lhePdgId[i])>0 and abs(_lhePdgId[i]) < 7));

        if(hasIncomingMother and _lheStatus[i]==1 and quarkOrGluon) _lheHTIncoming += _lhePt[i]; // To be used when an inclusive MC sample needs to be combined with HT-splitted samples
        if(_lhePdgId[i]==9900012 or _lhePdgId[i]==9990012) _ctauHN = lheEventInfo->hepeup().VTIMUP[i]; // the Heavy Neutrino lifetime
        if(abs(_lhePdgId[i])==15) ++_nTau;
    }

    tauCounter->Fill(_nTau, _weight);

    /* 
    TODOs:
    - Need to make sure the weighs saved in lheWeights are sensible:
        - muR, muF=1 is central, 1001 variation
    - Ideally 2 branches: dynamic scale and EFT
    - dynamic scale could be 4x9 2D vector?
    - 
    */
    if (multilepAnalyzer->storeEFTInfo) {
        // Fill EFT weights. Change this limit to something dynamic based on the weightgroups.
        _nEFTWeights = maxNumberOfEFTWeights;
        // CHECK: will add a blank in vector probablky
        for (unsigned i = 0; i < _nEFTWeights; i++) {
            if (lheEventInfo->weights()[i].id == "dummy") continue;
            _eftWeight[i] = lheEventInfo->weights()[i].wgt;
            std::cout << lheEventInfo->weights()[i].id << std::endl;
            eftCounter->Fill(i+0.5, _eftWeight[i]*_weight)
        }

        // Fill dynamic scale weights, this also includes the additional explicit alpha_s variations in the EFT MiniAODs from Robert
        // Saved
        unsigned index_true = offsetNominalScaleVar;
        _nDynScaleWeights = maxNumberOfDynScaleWeights;
        for (unsigned i = 0; i < _nDynScaleWeights-2; i++) {
            if (i % stepNominalScaleVar == 0) index_true++;
            _dynScaleWeight[i] = lheEventInfo->weights()[index_true].wgt/lheEventInfo->originalXWGTUP();
            dynScaleCounter->Fill(i+0.5, _dynScaleWeight[i]*_weight)
            std::cout << lheEventInfo->weights()[index_true].id << std::endl;

            index_true++;
        }
        _dynScaleWeight[maxNumberOfDynScaleWeights-2] = lheEventInfo->weights()[82].wgt/lheEventInfo->originalXWGTUP();
        _dynScaleWeight[maxNumberOfDynScaleWeights-1] = lheEventInfo->weights()[83].wgt/lheEventInfo->originalXWGTUP();
        

        // Fill the nominal LHE weights. If there are dynamic scale weights, these are intertwined with nominal Q2 variations, ie need this step
        // Also added an offset to take into account the EFT weights being saved first.
        _nLheWeights = maxNumberOfLheWeights;
        unsigned nominalCounter = 0;
        for(nominalCounter = 0; nominalCounter < _nLheWeights; nominalCounter++){
            _lheWeight[nominalCounter] = lheEventInfo->weights()[nominalCounter * stepNominalScaleVar + offsetNominalScaleVar].wgt/lheEventInfo->originalXWGTUP();
            lheCounter->Fill(nominalCounter + 0.5, _lheWeight[nominalCounter]*_weight);
            std::cout << lheEventInfo->weights()[nominalCounter * stepNominalScaleVar + offsetNominalScaleVar].id << std::endl;

            nominalCounter++;
        }

        // PDF weights
        // Should always be saved after all other weight variations so an offset is sufficient.
        for(nominalCounter; nominalCounter < _nLheWeights; ++nominalCounter){
            _lheWeight[nominalCounter] = lheEventInfo->weights()[nominalCounter+offsetPDFVar].wgt/lheEventInfo->originalXWGTUP();
            lheCounter->Fill(nominalCounter + 0.5, _lheWeight[nominalCounter]*_weight);
            std::cout << lheEventInfo->weights()[nominalCounter+offsetPDFVar].id << std::endl;

            nominalCounter++;
        }
    } else {
        //Store LHE weights to compute pdf and scale uncertainties, as described on https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW
        _nLheWeights = std::min( maxNumberOfLheWeights, static_cast< unsigned >( lheEventInfo->weights().size() ) );
        //most samples: 9 q2 weights + 101 pdfs + 2 alpha s (=pdf variation 101 and 102) = 112 || new ttg and maybe other exceptions: 45 q2 weights -> 148 total
        for(unsigned i = 0; i < _nLheWeights; ++i){
            _lheWeight[i] = lheEventInfo->weights()[i].wgt/lheEventInfo->originalXWGTUP();
            lheCounter->Fill(i + 0.5, _lheWeight[i]*_weight);
        }

        _nDynScaleWeights = 0;
        _nEFTWeights = 0;
    }

    //some tests for PS weight extraction
    const std::vector<double>& psWeights = genEventInfo->weights();
    _nPsWeights = std::min( maxNumberOfPsWeights, static_cast< unsigned >( psWeights.size() ) );
    for(unsigned ps = 0; ps < _nPsWeights; ++ps){
        _psWeight[ps] = psWeights[ps]/_weight;
        psCounter->Fill(ps + 0.5, _psWeight[ps]*_weight);
    }
}

double LheAnalyzer::getWeight() const{
    if( multilepAnalyzer->isData() ) return 1.;
    return _weight;
}

