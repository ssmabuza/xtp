#include "qmnblist.h"
#include "qmbead.h"
#include <votca/csg/nblist.h>
#include "qmtopology.h"

void QMNBList::Generate(BeadList &list1, BeadList &list2, bool do_exclusions)
{
    Cleanup();
    
    _father = dynamic_cast<QMTopology*> (list1.getTopology());
    
    NBList nb;
    nb.setCutoff(_cutoff);
    nb.Generate(list1, list2, do_exclusions);

    for(NBList::iterator iter=nb.begin(); iter!=nb.end();++iter) {
        QMBead *b1=(QMBead*)((*iter)->first);
        QMBead *b2=(QMBead*)((*iter)->second);

        if(b1->getMolecule() == b2->getMolecule()) continue;
        if(b1->GetCrgUnit() == NULL || b2->GetCrgUnit() == NULL) continue;

        CrgUnit *crg1, *crg2;
        crg1 = b1->GetCrgUnit();
        crg2 = b2->GetCrgUnit();

        if(crg1->getId() > crg2->getId())
            swap(crg1, crg2);
        
        if(!FindPair(crg1, crg2))
            AddPair(new QMPair(crg1, crg2, _father));
    }
}

