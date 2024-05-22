
#include"molecularWeight.h"


MolWtVariables molecularWeightExchangeRxn(MolWtVariables mwv) {
    // update the molecularWeight after an exchange reaction

    // for readability, make the following assignments
    double Mi_1=mwv.Mi_1;
    double Mi_2=mwv.Mi_2;
    double Mi_1_old=mwv.Mi_1_old;
    double Mi_2_old=mwv.Mi_2_old;
    const double weightCutoff=mwv.weightCutoff;
    
    if (Mi_1+Mi_2 != Mi_1_old+Mi_2_old) {
        throw std::runtime_error("the code is not yet equipped to handle evolution of small molecule during exchange reaction...");
    }
    if ( (Mi_1+Mi_2)>=weightCutoff ) {
        // check what is tracked now and what was tracked before

        bool trackM1; // is the new fragment with mass M1 above the weight cutoff?
        bool trackM2; // is the new fragment with mass M2 above the weight cutoff?
        bool trackM1_old; // was the old fragment with mass M1_old above the weight cutoff?
        bool trackM2_old; // was the old fragment with mass M2_old above the weight cutoff?

        // ensure that Mi_1 is always greater or equal to Mi_2...
        if (Mi_2>Mi_1) {
            std::swap(Mi_1,Mi_2);
        }
        if (Mi_2_old>Mi_1_old) {
            std::swap(Mi_1_old,Mi_2_old);
        }

        if (Mi_1<weightCutoff) {
            trackM1=false; // the fragment with mass M1 is too small to show up in a GPC trace
        }
        else {
            trackM1=true; // the fragment with mass M1 is still large enough to show up in a GPC trace and be considered polymer
        }
        if (Mi_2<weightCutoff) {
            trackM2=false; // the fragment with mass M2 is too small to show up in a GPC trace
        }
        else {
            trackM2=true; // the fragment with mass M2 is still large enough to show up in a GPC trace and be considered polymer
        }
        trackM1_old=(Mi_1_old<weightCutoff) ? false : true;
        trackM2_old=(Mi_2_old<weightCutoff) ? false : true;

        // Case 1: M1 and M2 tracked before; M1 and M2 still tracked
        if (trackM1 && trackM2 && trackM1_old && trackM2_old) {
            // sumNi does not change
            // sumMiNi does not change
            mwv.sumMi2Ni+=Mi_1*Mi_1+Mi_2*Mi_2-Mi_1_old*Mi_1_old-Mi_2_old*Mi_2_old;
        }
        // Case 2: M1 and M2 tracked before; M1 not M2 tracked now
        else if (trackM1_old && trackM2_old && trackM1 && !trackM2) {
            mwv.sumNi--; // sumNi decreases by one
            mwv.sumMiNi-=Mi_2; // total tracked mass decreases by Mi2
            mwv.sumMi2Ni+=Mi_1*Mi_1-Mi_1_old*Mi_1_old-Mi_2_old*Mi_2_old;
        }
        // Case 3: M1 not M2 tracked before; M1 not M2 tracked now
        else if (trackM1_old && !trackM2_old && trackM1 && !trackM2) {
            // sumNi does not change
            mwv.sumMiNi+=Mi_1-Mi_1_old;
            mwv.sumMi2Ni+=Mi_1*Mi_1-Mi_1_old*Mi_1_old;
        }
        // Case 4: M1 not M2 tracked before; M1 and M2 tracked now
        else if (trackM1_old && !trackM2_old && trackM1 && trackM2) {
            mwv.sumNi++;
            mwv.sumMiNi+=Mi_1+Mi_2-Mi_1_old;
            mwv.sumMi2Ni+=Mi_1*Mi_1+Mi_2*Mi_2-Mi_1_old*Mi_1_old;
        }
        // Case 5: M1 not M2 tracked before; neither M1 nor M2 tracked now
        else if (trackM1_old && !trackM2_old && !trackM1 && !trackM2) {
            mwv.sumNi--;
            mwv.sumMiNi-=Mi_1_old;
            mwv.sumMi2Ni-=Mi_1_old*Mi_1_old;
        }
        // Case 6: Neither M1 nor M2 tracked before; M1 not M2 tracked now
        else if (!trackM1_old && !trackM2_old && trackM1 && !trackM2) {
            mwv.sumNi++;
            mwv.sumMiNi+=Mi_1;
            mwv.sumMi2Ni+=Mi_1*Mi_1;
        }
        // Case 7: Neither M1 nor M2 tracked before; neither M1 nor M2 tracked now
        else if (!trackM1_old && !trackM2_old && !trackM1 && !trackM2) {
            // do nothing
        }
        // if M1 is always greater than M2 (both old and new), then this is all the cases...
        else {
            throw std::runtime_error("something weird happening...");
        }
    }

    mwv.Mn=mwv.sumMiNi/mwv.sumNi; // Mn=sumMiNi/sumNi
    mwv.Mw=mwv.sumMi2Ni/mwv.sumMiNi; // Mw=sumMi2Ni/sumMiNi
    if (mwv.Mw!=0){
        mwv.PDI=mwv.Mw/mwv.Mn; // PDI = Mw/Mn
    }
    return mwv;
}

MolWtVariables molecularWeightBreakingRxn( MolWtVariables mwv) {
    // update the molecular weight and associated variables after a breaking reaction
    
    // for readability, make the following assignments
    double Mi_1=mwv.Mi_1;
    double Mi_2=mwv.Mi_2;
    const double weightCutoff=mwv.weightCutoff;

    if ( (Mi_1+Mi_2)>=weightCutoff) {
        // this is a breaking reaction and the original unbroken chain was tracked

        bool trackM1; // is the new fragment with mass M1 above the weight cutoff?
        bool trackM2; // is the new fragment with mass M2 above the weight cutoff?

        if (Mi_1<weightCutoff) {
            trackM1=false; // the fragment with mass M1 is too small to show up in a GPC trace
        }
        else {
            trackM1=true; // the fragment with mass M1 is still large enough to show up in a GPC trace and be considered polymer
        }
        if (Mi_2<weightCutoff) {
            trackM2=false; // the fragment with mass M2 is too small to show up in a GPC trace
        }
        else {
            trackM2=true; // the fragment with mass M2 is still large enough to show up in a GPC trace and be considered polymer
        }
        if (Mi_1 > 0.00001 && Mi_2 > 0.00001 ) {
            // the if statement checks whether Mi_1 or Mi_2 is 0. If one is 0, this means that the bond breakage did not actually result in a new chain forming
            if (trackM1 && trackM2) {
                mwv.sumNi++; // one new chain has been formed
                // sumMiNi does not change, provided that no small molecule was evolved
                mwv.sumMi2Ni+=Mi_1*Mi_1+Mi_2*Mi_2-(Mi_1+Mi_2)*(Mi_1+Mi_2); // add tracked fragments, decrement by original fragment
            }
            else if (trackM1 & !trackM2) {
                // sumNi does not change
                mwv.sumMiNi-=Mi_2; // decrement by the fragment that is no longer tracked
                mwv.sumMi2Ni+=Mi_1*Mi_1-(Mi_1+Mi_2)*(Mi_1+Mi_2); // decrement by original fragment, add only tracked fragment 
            }
            else if (!trackM1 & trackM2) {
                // sumNi does not change
                mwv.sumMiNi-=Mi_1; // decrement by the fragment that is no longer tracked
                mwv.sumMi2Ni+=Mi_2*Mi_2-(Mi_1+Mi_2)*(Mi_1+Mi_2); // decrement by original fragment, add only tracked fragment 
            }
            else {
                // if neither of the fragments are tracked
                mwv.sumNi--; // total number of tracked chains decreases by 1
                mwv.sumMiNi-=(Mi_1+Mi_2); // total weight decreases by both
                mwv.sumMi2Ni-=(Mi_1+Mi_2)*(Mi_1+Mi_2); // decrement by amt of original chain weight
            }
        }
    }

    mwv.Mn=mwv.sumMiNi/mwv.sumNi; // Mn=sumMiNi/sumNi
    mwv.Mw=mwv.sumMi2Ni/mwv.sumMiNi; // Mw=sumMi2Ni/sumMiNi
    if (mwv.Mw!=0){
        mwv.PDI=mwv.Mw/mwv.Mn; // PDI = Mw/Mn
    }

    return mwv;
}

MolWtVariables molecularWeightLinkingRxn(MolWtVariables mwv) {
    // update the molecular weight and associated variables after a bimolecular linking (polymerization) reaction
    
    // for readability, make the following assignments
    double Mi_1=mwv.Mi_1;
    double Mi_2=mwv.Mi_2;
    const double weightCutoff=mwv.weightCutoff;

    if ((Mi_1+Mi_2)>=weightCutoff) {
        bool trackM1; // is the new fragment with mass M1 above the weight cutoff?
        bool trackM2; // is the new fragment with mass M2 above the weight cutoff?
        if (Mi_1<weightCutoff) {
            trackM1=false; // the fragment with mass M1 is too small to show up in a GPC trace
        }
        else {
            trackM1=true; // the fragment with mass M1 is still large enough to show up in a GPC trace and be considered polymer
        }
        if (Mi_2<weightCutoff) {
            trackM2=false; // the fragment with mass M2 is too small to show up in a GPC trace
        }
        else {
            trackM2=true; // the fragment with mass M2 is still large enough to show up in a GPC trace and be considered polymer
        }
    
        // first scenario is if neither Mi_1 nor Mi_2 were tracked before
        if (trackM1==false && trackM2==false) {
            mwv.sumNi++;
            mwv.sumMiNi+=Mi_1+Mi_2;
            mwv.sumMi2Ni+=(Mi_1+Mi_2)*(Mi_1+Mi_2);
        }
        // second scenario is if M1 was tracked before and M2 was not tracked before
        else if (  trackM1 && !trackM2  ) {
            // sumNi does not change
            mwv.sumMiNi+=Mi_2; // Mi_2 is the weight of the new piece
            mwv.sumMi2Ni+=(Mi_1+Mi_2)*(Mi_1+Mi_2)-Mi_1*Mi_1; // newchainweight^2-oldchainweight^2
        }
        // third scenario is if M2 was tracked before and M1 was not tracked before
        else if ( !trackM1 && trackM2 ) {
            // sumNi does not change
            mwv.sumMiNi+=Mi_1; // Mi_1 is the weight of the additional chain fragment that wasn't tracked previously
            mwv.sumMi2Ni+=(Mi_1+Mi_2)*(Mi_1+Mi_2)-Mi_2*Mi_2; // newchainweight^2-oldchainweight^2
        }
        // fourth scenario is if M1 && M2 were tracked before
        else if ( trackM1 && trackM2 ) {
            // two chains reacted
            mwv.sumNi--;
            // sumMiNi does not change
            mwv.sumMi2Ni+=(Mi_1+Mi_2)*(Mi_1+Mi_2)-Mi_1*Mi_1-Mi_2*Mi_2;
        }
    
        mwv.Mn=mwv.sumMiNi/mwv.sumNi; // Mn=sumMiNi/sumNi
        mwv.Mw=mwv.sumMi2Ni/mwv.sumMiNi; // Mw=sumMi2Ni/sumMiNi
        if (mwv.Mw!=0){
            mwv.PDI=mwv.Mw/mwv.Mn; // PDI = Mw/Mn
        }
    }
    return mwv;
}
MolWtVariables molecularWeight(MolWtVariables mwv) {
    // call the appropriate molecular weight update function and reset all trackers
    if (mwv.sameChain==true) {
        mwv.sameChain=false;
        if (mwv.monomerAdded || mwv.isNewChain) {
            throw std::runtime_error("molwt bools dont make sense");
        }
        return mwv; // nothing changes for an intramolecular reaction in which no small molecule is evolved
    }
    else if (mwv.rxnType=="linking") {
        mwv=molecularWeightLinkingRxn(mwv);
    }
    else if (mwv.rxnType=="breaking") {
        mwv=molecularWeightBreakingRxn(mwv);
    }
    else if (mwv.rxnType=="exchange") {
        // need both the new Mi_1 and Mi_2, but also the old Mi_1 and old Mi_2
        mwv=molecularWeightExchangeRxn(mwv);
    }
    mwv.sameChain=false;
    mwv.monomerAdded=false;
    mwv.isNewChain=false;
    if (mwv.sumNi<0 || mwv.sumMiNi < 0 || mwv.sumMi2Ni < 0) {
        throw std::runtime_error(" this should not happen man");
    }
    return mwv;
}