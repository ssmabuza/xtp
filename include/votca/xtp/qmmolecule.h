/*
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef VOTCA_XTP_QMMOLECULE_H
#define VOTCA_XTP_QMMOLECULE_H

#include <votca/xtp/atomcontainer.h>


namespace votca {
    namespace xtp {

class QMMolecule : public AtomContainer<QMAtom>
{
public:
    QMMolecule(const string& name,int id):_id(id),_name(name){};
    
    void WriteToCpt(CptLoc parent)const;

    void ReadFromCpt(CptLoc parent);

    void 
    
    
    
    
};
        
        
        
        
    }
}

#endif /* VOTCA_XTP_QMMOLECULE_H */
