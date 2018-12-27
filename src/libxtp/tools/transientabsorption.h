/* 
 *            Copyright 2009-2019 The VOTCA Development Team
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

#ifndef _VOTCA_XTP_TRANSIENTABSORPTION_H
#define _VOTCA_XTP_TRANSIENTABSORPTION_H

#include <stdio.h>
#include <votca/ctp/logger.h>
#include <boost/filesystem.hpp>

namespace votca { namespace xtp {
    using namespace std;
    
class TransientAbsorption : public ctp::QMTool
{
public:

    TransientAbsorption () { };
   ~TransientAbsorption () { };

    string Identify() { return "transientabsorption"; }

    void   Initialize(Property *options);
    bool   Evaluate();
    

private:
    
    string      _orbfile;
    string      _output_file;
    
    /* TODO: add variables that are required
     
     */
    
    
    
    ctp::Logger      _log;
    
    
};

void TransientAbsorption::Initialize(Property* options) {
    
            // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options, "xtp" );
    string key = "options." + Identify();
 
    _orbfile      = options->get(key + ".input").as<string> ();
    _output_file  = options->get(key + ".output").as<string> ();
    
    /* TODO: add options for the calculation of TA spectrum
     - which GW-BSE states to include
    */
    
    
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");    
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");

}

bool TransientAbsorption::Evaluate() {
    
    _log.setReportLevel( ctp::logDEBUG );
    _log.setMultithreading( true );
    
    _log.setPreface(ctp::logINFO,    "\n... ...");
    _log.setPreface(ctp::logERROR,   "\n... ...");
    _log.setPreface(ctp::logWARNING, "\n... ...");
    _log.setPreface(ctp::logDEBUG,   "\n... ..."); 

    CTP_LOG(ctp::logDEBUG, _log) << "Converting serialized QM data in " << _orbfile << flush;

    Orbitals _orbitals;
    // load the QM data from serialized orbitals object

    CTP_LOG(ctp::logDEBUG, _log) << " Loading QM data from " << _orbfile << flush;
    _orbitals.ReadFromCpt(_orbfile);

    
    /* TODO: actual calculation
     *  1) calculate free transition dipoles <v'|d|v>, possibly put function 
     *     in BSE object (later)
     *  2) calculate free transition dipoles <c'|d|c>, possibly put function
     *     in BSE object (later)
     *  3) calculate oscillator strength according (TDA or full BSE), possibly 
     *     put function in BSE object (later)
    */
  
    CTP_LOG(ctp::logDEBUG, _log) << "Written charges to " << _output_file << flush;
    
    return true;
}








}}


#endif