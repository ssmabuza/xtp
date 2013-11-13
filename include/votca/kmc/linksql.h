/*
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _VOTCA_KMC_LINKSQL_H
#define	_VOTCA_KMC_LINKSQL_H

#include <votca/tools/vec.h>
#include <votca/kmc/link.h>

namespace votca { namespace kmc {


/**
 * \brief A link between two nodes
 * 
 * Container for pair properties: rates, couplings, separations 
 * 
 */
class LinkSQL: public Link
{

public:

    LinkSQL( int id, Node* node1, Node* node2, votca::tools::vec r12) : Link(id,node1,node2,r12){
    };

    void setRate(double rate12e, double rate12h, double rate21e, double rate21h){
        _rate12e  = rate12e;
        _rate12h  = rate12h;
        _rate21e  = rate21e;
        _rate21h  = rate21h;
    }
         
    void setJeff2(double Jeff2e, double Jeff2h) {
        _Jeff2e  = Jeff2e;
        _Jeff2h  = Jeff2h;
    }
    
    void setlO(double lOe, double lOh) {
        _lOe = lOe;
        _lOh = lOh;
    }

    
/*    /// forward and backward rates
    const double &rate12() const { return _rate12; }
    const double &rate21() const { return _rate21; }
    
    /// square of effective transfer integral
    const double &jeff2() const { return _jeff2; }
    
    /// outer reorganization energy
    const double &lambda0() const { return _lambda0; }*/

private:

    
    /// forward and backward rates
    double _rate12e;
    double _rate12h;
    double _rate21e;
    double _rate21h;
    
    double _Jeff2e;
    double _Jeff2h;
    
    double _lOe;
    double _lOh;
    
};

}} 

#endif // _VOTCA_KMC_LINKSQL_H
