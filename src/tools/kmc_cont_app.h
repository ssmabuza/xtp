/* 
 * File:   kmc_cont_app.h
 * Author: vehoff
 *
 * Created on March 4, 2010, 4:50 PM
 */

#ifndef _KMC_CONT_APP_H
#define	_KMC_CONT_APP_H

#include "qmapplication.h"
#include <kmc/vertex.h>
#include <kmc/hoppers.h>
#include <kmc/kmc.h>
#include <kmc/graph.h>

class KmcCont : public QMApplication
{
public:
    KmcCont();
    ~KmcCont();

    void HelpText();
    void Initialize();
    bool EvaluateFrame();

private:
    /// electric field
    vec _E;
    /// total simulation time
    double _total_time;
    /// time step
    double _dt;
    /// number of charges simultaneously present in a simultion
    int _ncrg;
    /// number of KMC runs to be performed
    int _nruns;
    ///  output streams for velocity averaging & diffusion
    ofstream _out_cont;
    ofstream _out_diff;

    /// creation of KMC graph

    void make_kmc_graph(graph *a, QMNBList &nblist) {
        cout << "[make_kmc_graph]: Building KMC Graph...";
        /// assign constants
        a->SetField(_E);
        /// set vertices equal to centers of mass
        list < CrgUnit *> listCharges = _qmtop.crglist();
        list < CrgUnit *>::iterator it;
        for (it = listCharges.begin(); it != listCharges.end(); ++it) {
            a->AddVertex((*it)->GetCom(), _E); /// TO DO: remove necessity for E-field at this point
        }
        /// set edges, two edges 1->2 and 2->1 are created per neighboring pair
        for (QMNBList::iterator iter = nblist.begin(); iter != nblist.end(); ++iter) {
            a->AddEdge((*iter)->first->getId(), (*iter)->second->getId(), (*iter)->rate12(), (*iter)->r());
            a->AddEdge((*iter)->second->getId(), (*iter)->first->getId(), (*iter)->rate21(), -(*iter)->r());
        }
        cout << " Done." << endl;
    }
};

#endif	/* _KMC_CONT_APP_H */

