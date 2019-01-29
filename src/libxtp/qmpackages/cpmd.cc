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

#include "cpmd.h"

#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/tools/constants.h>
#include <stdio.h>
#include <iomanip>
#include <sys/stat.h>
#include <vector>
#include <functional>

#ifdef DEBUG
#include <votca/xtp/aomatrix.h>
#endif



namespace votca {
    namespace xtp {
        //using namespace std;

        void Cpmd::Initialize(tools::Property &options) {

            // CPMD file names
            std::string fileName = "system";
            _xyz_file_name = fileName + ".xyz";
            _input_file_name = fileName + ".inp";
            _log_file_name = fileName + ".log";
            _orb_file_name = "WFNCOEF" ;

            
            std::string key = "package";
            std::string _name = options.get(key + ".name").as<std::string> ();

            if (_name != "cpmd") {
                std::cerr << "Tried to use " << _name << " package. ";
                throw std::runtime_error("Wrong options file");
            }
            
            
            _executable =       options.get(key + ".executable").as<std::string> ();
            _charge =           options.get(key + ".charge").as<int> ();
            _spin =             options.get(key + ".spin").as<int> ();
            _options =          options.get(key + ".options").as<std::string> ();
            if (options.exists(key + ".memory")) { //no way to specify this to CPMD
                XTP_LOG(logWARNING, *_pLog) << "CPMD: No way to specify maximum "
                        "memory available to CPMD. "
                        "Ignoring the memory option\n"<< std::flush;
            }
            _threads = options.get(key + ".threads").as<int> ();
            _scratch_dir = options.get(key + ".scratch").as<std::string> ();
            _cleanup = options.get(key + ".cleanup").as<std::string> ();
            if(_threads!=1){
                XTP_LOG(logWARNING, *_pLog) << "CPMD: CPMD doesn not fully "
                        "support OpenMP paralelization..\n"
                        "If you want a parallel run of CPMD, use MPI. "
                        "Set the executable name to: mpirun -np <# of threads> cpmd-mpi.x"
                        "Setting number of threads to 1 and continuing." << std::flush;
                _threads=1;
            }



            //restart?
            if (options.exists(key + ".restart")) {
                _rsrt=true;
                _rsrt_kwds = options.get(key + ".restart").as<std::string> ();
            }
            else _rsrt=false;

            //optimize wavefunction?
            if (options.exists(key + ".optimizewf")) {
                _optWF=true;
                _convCutoff = options.get(key + ".optimizewf").as<double> ();
            }
            else _optWF=false;
            
            //use VDW correction?
            if (options.exists(key + ".useGrimmeVDW")) {
                _useGrimmeVDW=true;
            }
            else _useGrimmeVDW=false;
            
            //custom CPMD controls
            if (options.exists(key + ".customCPMDcontrolls")) {
                _custom_CPMD_controlls=options.get(key + ".customCPMDcontrolls").as<std::string> ();
                if(_custom_CPMD_controlls.find("LSD")!=std::string::npos){
                    XTP_LOG(logERROR, *_pLog) << "CPMD: LSD option not supported." << std::endl <<
                            "XTPlib will not be able to read in basisfunction " <<
                            "coefficients after CPMD is finished with an LSD-enabled calculation." << std::endl << std::flush;
                    throw std::runtime_error("LSD is currently not supported");
                }
            }
            else _custom_CPMD_controlls="";


            //functional and pseudopotentials
            if (options.exists(key + ".functional")) {
                if (!(options.exists(key + ".functional"))) throw std::runtime_error("Functional name missing");
                _functional = options.get(key + ".functional.name").as<std::string> ();
                std::list<tools::Property*>::iterator pit;
                std::list<tools::Property *> props = options.Select(key + ".functional.pseudopotentials.*");
                for (pit = props.begin(); pit != props.end(); ++pit)
                {
                    tools::Property* p= *pit;
                    _ppFileNames[p->name()]=p->as<std::string> ();
                }
                props = options.Select(key + ".functional.l.*");
                for (pit = props.begin(); pit != props.end(); ++pit)
                {
                    tools::Property* p= *pit;
                    _ppLData[p->name()]=p->as<std::string> ();
                }
                
            }
            else throw std::runtime_error("No functional and pseudopotentials specified");
            
            
            // check if ecp is specified
            if (!options.exists(key + ".ecp")) {
                _ecp_name = "corelevels.xml";
                XTP_LOG(logDEBUG, *_pLog) << "Effective core potential not specified."
                        << " Setting to "<< _ecp_name
                        << " as those should be compatible with CPMD."
                        << " If they aren't, construct a new xml file from the _ZV vector"
                        << " read from WFNCOEF output file of CPMD."<< std::flush;
            }else{    
                _ecp_name = options.get(key + ".ecp").as<std::string> ();
            }
            
            
            //symmetry
            _symmetry=0;
            if (options.exists(key + ".symmetry")) {
                _symmetry = options.get(key + ".symmetry").as<int> ();
            }
            else
                XTP_LOG(logDEBUG, *_pLog) << "CPMD: no symmetry provided, assuming simple cubic." << std::flush;
            
            //cell
            _cell="20.0   1.0   1.0  0.0  0.0  0.0";
            if (options.exists(key + ".cell")) {
                _cell = options.get(key + ".cell").as<std::string> ();
            }
            else
                XTP_LOG(logWARNING, *_pLog) << "CPMD: no cell provided, assuming cube with side length of 20 Bohr." << std::flush;
            
            //box for wrapping atomic coordinates
            XTP_LOG(logDEBUG, *_pLog) << "Attempting to parse box dimensions"<< 
                " cell definition. Assuming right-angled box & no ABSOLUTE."<<std::flush;
            std::vector<std::string> cell_entries;
            boost::algorithm::split(cell_entries, _cell, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
            double a,sb,sc;
            a = boost::lexical_cast<double>(cell_entries[0]);
            sb= boost::lexical_cast<double>(cell_entries[1]);
            sc= boost::lexical_cast<double>(cell_entries[2]);
            _box=tools::vec(a, a*sb, a*sc);
            XTP_LOG(logDEBUG, *_pLog) << "CPMD: will wrap atoms to a "
                <<_box[0]<<" "<<_box[1]<<" "<<_box[2]
                <<" Bohr box when reading CPMD output."<<std::flush;
            
            
            //plane wave cutoff
            _pwCutoff=80.0;
            if (options.exists(key + ".pwcutoff")) {
                _pwCutoff = options.get(key + ".pwcutoff").as<double> ();
            }
            else
                XTP_LOG(logDEBUG, *_pLog) << "CPMD: no plane wave cutoff provided, assuming "<< _pwCutoff <<" Ry." << std::flush;

            //output electrostatic potential?
            _elpot=false;
            if (options.exists(key + ".elpot")) {
                _elpot=options.get(key + ".elpot").as<bool> ();
            }

            //project wavefunction onto atomic orbitals?
            _projectWF=false;
            if (options.exists(key + ".projectwf")) {
                _projectWF=options.get(key + ".projectwf").as<bool> ();
            }
            
            //basis set file name
            if (options.exists(key + ".basisset")) {
            _basisset_name=options.get(key + ".basisset").as<std::string> ();
            }

            //do population analysis? required to have access to WF coefs in atom-centric basis
            //requires _projectWF
            _popAnalysis=false;
            if (options.exists(key + ".popanalysis")) {
                _popAnalysis=options.get(key + ".popanalysis").as<bool> ();
            }

            //get density and overlap matrices during post processing?
            //requires _popAnalysis and _projectWF
            _getMat=false;
            if (options.exists(key + ".getmatrices")) {
                _getMat=options.get(key + ".getmatrices").as<bool> ();
            }
            
            //Sizes of Fortran int and real CPMD was compiled with.
            //These are needed to correctly read the binary overlap file
            _F_int_size=4;
            _F_real_size=8;
            if (options.exists(key + ".Fortran_int_size")) {
                _F_int_size=options.get(key + ".Fortran_int_size").as<int> ();
            }
            if (options.exists(key + ".Fortran_real_size")) {
                _F_real_size=options.get(key + ".Fortran_real_size").as<int> ();
            }
            

            if(_getMat) _popAnalysis=true;
            if(_popAnalysis) _projectWF=true;
            
            if (options.exists(key + ".outputVxc")) {
                _output_Vxc = options.get(key + ".outputVxc").as<bool> ();
                XTP_LOG(logERROR, *_pLog) << "Error: CPMD interface currently does not support Vxc." << std::flush;
                    throw std::runtime_error("Vxc not supported.\n");
            } else _output_Vxc = false;

            if(_projectWF && _optWF){
                XTP_LOG(logDEBUG, *_pLog) << "CPMD: Splitting run into two steps." << std::endl << std::flush;
                std::cerr << "Warning: Wavefunction optimization and projection onto atom-centric orbitals can not be done together.\nCPMD would crash.\n";
                std::cerr << "Splitting into optimization on first run and projection/population analysis on second.";
            }

        }

        bool Cpmd::WriteInputFile(const Orbitals& orbitals){
            
            std::ofstream _com_file;

            std::string _com_file_name_full = _run_dir + "/" + _input_file_name;
            //if need to run CPMD twice, create new names for input and log files of the first run
            if(_projectWF && _optWF){
                //input
                size_t pointpos = _input_file_name.find('.');
                if (pointpos!=std::string::npos){
                    _wfOpt_input_file_name=_input_file_name.substr(0,pointpos) + "_wfOpt" + _input_file_name.substr(pointpos);
                }else{
                    _wfOpt_input_file_name=_input_file_name + "_wfOpt";
                }
                _com_file_name_full= _run_dir + "/" + _wfOpt_input_file_name;

                //log
                pointpos = _log_file_name.find('.');
                if (pointpos!=std::string::npos){
                    _wfOpt_log_file_name=_log_file_name.substr(0,pointpos) + "_wfOpt" + _log_file_name.substr(pointpos);
                }else{
                _wfOpt_log_file_name=_log_file_name + "_wfOpt";
                }
            }

            _com_file.open(_com_file_name_full.c_str());

            // header
            _com_file << "&INFO\nGenerated by VOTCA\n&END" << std::endl;


            //control
            _com_file << "\n&CPMD\n";
            if(_rsrt) _com_file << "  RESTART " << _rsrt_kwds << std::endl;  //restart
            if(_optWF){                                                 //optimize WF
                _com_file << "  OPTIMIZE WAVEFUNCTION" << std::endl;
                _com_file << "  CONVERGENCE ORBITALS" << std::endl;
                _com_file << "  " << FortranFormat(_convCutoff) << std::endl;
                _com_file << "  PCG MINIMIZE" << std::endl;              //use the more stable optimizer
                _com_file << "  TIMESTEP" << std::endl;
                _com_file << "   20" << std::endl;
                _com_file << "  STORE WAVEFUNCTIONS" <<std::endl;
                _com_file << "   20 SC=20" <<std::endl;


            }
            if(_elpot){                                                 //output electrostatic potential
                _com_file << "  ELECTROSTATIC POTENTIAL" << std::endl;
                _com_file << "  RHOOUT" << std::endl;
            }
            if(_projectWF && !_optWF){
                _com_file << "  PROPERTIES" << std::endl;
            }
            
            _com_file << std::endl <<"  " << _custom_CPMD_controlls << std::endl;        //custom CPMD controls from the .xml file
            if(_useGrimmeVDW){                                          //VDW
                _com_file << "  VDW CORRECTION" << std::endl;
            }
            _com_file << "&END" << std::endl;

            //functional
            _com_file << "\n&DFT\n";
            _com_file << "  FUNCTIONAL " << _functional << std::endl;
            _com_file << "  GC-CUTOFF" << std::endl;
            _com_file << "   1.0d-06" << std::endl;
            _com_file << "&END" << std::endl;
            
            //VDW
            if(_useGrimmeVDW){  
                _com_file << "\n&VDW" << std::endl;
                _com_file << " VDW CORRECTION" << std::endl;
                _com_file << " ALL GRIMME" << std::endl;
                _com_file << " VDW-CELL" << std::endl;
                _com_file << " 1 1 1" << std::endl;
                _com_file << "&END" << std::endl;
            }
            
            //cell
            _com_file << "\n&SYSTEM\n";
            _com_file << "  SYMMETRY" << std::endl;
            _com_file << "   " << _symmetry <<std::endl;
            _com_file << "  CELL" << std::endl;
            _com_file << "   " << _cell <<std::endl;
            _com_file << "  CUTOFF" << std::endl;
            _com_file << "   " << FortranFormat(_pwCutoff) <<std::endl;
            _com_file << "&END" << std::endl;
            
            //properties
            if(_projectWF){
                _com_file << "\n&PROP\n";
                _com_file << "  PROJECT WAVEFUNCTION" << std::endl;
                if(_popAnalysis){
                    _com_file << "  CHARGES" << std::endl;
                    _com_file << "  POPULATION ANALYSIS MULLIKEN" << std::endl;
                }
                _com_file << "&END" << std::endl;
            }
            
            //basis
            if(_projectWF && !_optWF){ //prevent initializing _optWF runs with this basis
                _com_file << "\n&BASIS\n";
                WriteBasisSet(orbitals, _com_file);
                _com_file << "&END" << std::endl;
            }
            
            //atoms
            _com_file << "\n&ATOMS\n";
            //find how many atoms of each element there are
            const QMMolecule& qmatoms = orbitals.QMAtoms();     
            std::vector<std::string> UniqueElements = qmatoms.FindUniqueElements();
            _elements = UniqueElements; //needed in reading log file
            int numAtoms = qmatoms.size();
            
            //find number of atoms of each element
            _nAtomsOfElement.clear();
            for(const std::string& element_name:UniqueElements){
                _nAtomsOfElement.insert(std::make_pair(element_name, 0));
                for(const QMAtom& a:qmatoms){
                    if(a.getElement() == element_name){
                        _nAtomsOfElement[element_name] += 1;
                    }
                }
            }
            
            
            //now loop over elements and store all atoms of that element
            bool atomOrderMapSet = (VOTCA2CPMD_map!=NULL && CPMD2VOTCA_map!=NULL);
            int Vind=0, Cind=0; //atom indexes in VOTCA and CPMD
            if(!atomOrderMapSet){
                VOTCA2CPMD_map=new int[numAtoms];
                CPMD2VOTCA_map=new int[numAtoms];
            }

            for (const std::string& element_name:UniqueElements) {
                if(_ppFileNames.find(element_name)==_ppFileNames.end()) {
                    XTP_LOG(logERROR, *_pLog) << "Error: Element "<<element_name<<
                            " has no pseudopotential specified in CPMD options file." << std::flush;
                    throw std::runtime_error("Encountered element with no pseudopotential.\n");
                }
                else{ //store name of the pseudopotential file and it's read options
                    _com_file << "*" << _ppFileNames[element_name] << std::endl; 
                    if(_ppLData.find(element_name)==_ppLData.end()) {
                        XTP_LOG(logWARNING, *_pLog) << "Element "<<element_name<<
                                " has no angular momentum data (<l></l>) specified in CPMD options file. "<<
                                "Attempting to read it from basis set. This may produce errors." << std::flush;
                        if(_basisset_name.empty()){
                            XTP_LOG(logERROR, *_pLog) << "Error: Basis set file not specified." << std::flush;
                            throw std::runtime_error("Encountered element with no angular momentum data and no basis set is specified.\n");
                        }
                        else{
                            BasisSet _bs;
                            _bs.LoadBasisSet(_basisset_name);
                            int Lmax = 0;
                            //find Lmax by checking all shells of the element
                            const Element& element = _bs.getElement(element_name);
                            for (const Shell& shell:element) {
                                int Ls = shell.getLmax();
                                if(Lmax<Ls) Lmax = Ls;
                            }
                            _com_file << "   "<<Lmax<<" "<<Lmax<<" "<<Lmax<< std::endl; //LMAX LOC SKIP
                        }
                    }
                    else{
                        _com_file << "   "<<_ppLData[element_name]<< std::endl; //LMAX LOC SKIP
                    }
                    _com_file << "   "<< _nAtomsOfElement[element_name] <<std::endl;  //# atoms of element
                    
                    //store atomic positions of atoms of this element
                    for(const QMAtom& a:qmatoms){
                        Vind=0;
                        if(a.getElement() == element_name){ //this element
                            Eigen::Vector3d pos = a.getPos(); //in Bohr
                            _com_file << "   ";
                            _com_file << std::setw(12) << std::setiosflags(std::ios::fixed) << std::setprecision(5) << pos.x() << "   ";
                            _com_file << std::setw(12) << std::setiosflags(std::ios::fixed) << std::setprecision(5) << pos.y() << "   ";
                            _com_file << std::setw(12) << std::setiosflags(std::ios::fixed) << std::setprecision(5) << pos.z() << "   ";
                            _com_file << std::endl;

                            //cache the mapping between VOTCA and CPMD atomic ordering
                            if(!atomOrderMapSet){
                                VOTCA2CPMD_map[Vind]=Cind;
                                CPMD2VOTCA_map[Cind]=Vind;
                                CPMD2TYPE_map[Cind]=(element_name);
                            }
                            Cind++;
                        }
                        Vind++;
                    }
                }
            }
            
            //#warning "TODO: copy pseudopotentials to the _run_dir"
            _com_file << "&END" << std::endl;
            
            
            
            _com_file << std::endl;
            _com_file.close();
            
            
            //now write the output file for the second run, if necessary
            if(_projectWF && _optWF){
                _optWF=false;
                //force reading of restart file on second run
                bool old_rsrt=_rsrt;
                std::string old_rsrt_kwds = _rsrt_kwds;
                _rsrt = true;
                _rsrt_kwds = "WAVEFUNCTION COORDINATES LATEST";
                WriteInputFile(orbitals);
                //reset variables
                _optWF     = true;
                _rsrt      = old_rsrt;
                _rsrt_kwds = old_rsrt_kwds;
            }

            return true;
        }
        
        
        /**
         * Writes the basis set files to disk in a format that CPMD can understand
         */
        void Cpmd::WriteBasisSet(const Orbitals& orbitals, std::ofstream &_com_file) {
            
            const QMMolecule& qmatoms = orbitals.QMAtoms();            
            std::vector<std::string> UniqueElements = qmatoms.FindUniqueElements();
            BasisSet _bs;
            _bs.LoadBasisSet(_basisset_name);
            XTP_LOG(logDEBUG, *_pLog) << "Loaded Basis Set " << _basisset_name << std::flush;
            
            for (const std::string& element_name:UniqueElements) {
                XTP_LOG(logDEBUG, *_pLog) << "CPMD: writing Gaussian basis for element "<< element_name << std::flush;
               
                const Element& element = _bs.getElement(element_name);
                
                std::string _short_el_file_name = element_name + "_" + _basisset_name + ".basis";
                std::string _el_file_name = _run_dir + "/" + _short_el_file_name;


                //write the element to the input file
                _com_file << "*" << _short_el_file_name << " " << element.NumOfShells() << " GAUSSIAN"<<std::endl;
                _com_file << "   ";
                
                //create the .basis file
                std::ofstream _el_file;
                _el_file.open(_el_file_name.c_str());

                //comment
                _el_file << element_name << " with the "<< _basisset_name << " basis." << std::endl;
                
                //Lmax
                int Lmax = 0;
                //find Lmax by checking all shells of the element
                for (const Shell& shell:element) {
                            int Ls = shell.getLmax();
                            if(Lmax<Ls) Lmax = Ls;
                }
                _el_file << Lmax+1 << std::endl; //number of L-values in this .basis file

                //sort shells by L
                for (int L=0; L <= Lmax; L++)
                {
                    std::vector< std::reference_wrapper<const Shell> > Lshells;
                    int ndecays=0;
                    for (const Shell& shell:element) {
                        int Ls = shell.getLmax();
                        if(shell.getType().size()>1){
                            XTP_LOG(logDEBUG, *_pLog) << "CPMD does not support mixed basis function like " << shell.getType() << "." << std::flush;
                            XTP_LOG(logDEBUG, *_pLog) << "Please break the basis set into basis functions with only one L-value each." << std::flush;
                            XTP_LOG(logERROR, *_pLog) << "CPMD: multi-L basis functions not supported." << std::flush;
                            throw std::runtime_error("Unsupported basis function");
                        }

                        //For now assume all shells have only one L-value.
                        //Can decompose the basis set into such shells later, in another place.
                        //Any subsequent analysis will have to use the decomposed basis set too.
                        if (Ls==L) //this shell has the correct L
                        {
                            ndecays += shell.getSize();
                            Lshells.push_back(shell);
                            //write the shell's L-value to the input file
                            _com_file << L << " ";
                        }       
                    }
                    
                    if(!Lshells.empty()){   //only write things if there are shells with this L
                        _el_file << "  Functions for l="<< L << std::endl;
                        _el_file << "  " << Lshells.size() << " " << ndecays << std::endl;
                        _el_file << std::endl;

                        //decays
                        std::ios::fmtflags old_settings = _el_file.flags();
                        _el_file << std::scientific << std::setprecision(6);
                        _el_file << "  ";
                        for (const Shell& s:Lshells)
                        {
                            for (const GaussianPrimitive& g:s) {
                                _el_file << g._decay << "\t";
                            }
                        }
                        _el_file << std::endl;

                        //coefficients (scale*contraction)
                        int gs=0; //number of decays already handled
                        for (const Shell& s:Lshells)
                        {
                            if(s.getSize()!=0) //there are gaussians in this shell
                            {
                                int gi=0; //index of the current decay
                                _el_file << "  ";
                                //output zeros for all decays already handled
                                for(gi=0; gi<gs; gi++)
                                {
                                    _el_file << 0.0 << "\t";
                                }
                                //output coefficients for this shell's gaussians
                                for (const GaussianPrimitive& g:s) {
                                    _el_file << s.getScale() * g._contraction[L] << "\t";
                                    gi++;
                                }
                                gs+=s.getSize();
                                //output zeros till the end of decays
                                for(;gi<ndecays; gi++)
                                {
                                    _el_file << 0.0 << "\t";
                                }
                                _el_file << std::endl;
                            }
                        }
                        _el_file.flags(old_settings);
                    }
                }

                _el_file << std::endl;
                _el_file.close();       //end the .basis file

                _com_file<<std::endl;   //next element definition in the .inp file
            }
        }
        

        /**
         * Runs the CPMD job.
         */
        bool Cpmd::Run() {
            XTP_LOG(logDEBUG, *_pLog) << "CPMD: Run()" << std::flush;
            if(_optWF && _projectWF){ //CPMD needs to run twice, once for _optWF and once for _projectWF
                //_optWF run:
                 XTP_LOG(logDEBUG, *_pLog) << "CPMD: running [" << _executable << " " << _wfOpt_input_file_name << "]" << std::flush;
                 if (std::system(NULL)) {
                    std::string _command;
                    _command = "cd " + _run_dir + "; rm -f LocalError*.log; " + _executable + " " + _wfOpt_input_file_name + " > " + _wfOpt_log_file_name;
                    
                    int check=std::system(_command.c_str());
	            if (check==-1){
    	                XTP_LOG(logERROR, *_pLog) << _input_file_name << " failed to start" << std::flush;
    	                return false;
    	            }
                    
                    check = std::system("mv LocalError*.log LocalError_wfOpt*.log");
                    if (check==0){
    	                XTP_LOG(logWARNING, *_pLog) << "CPMD produced an error log. Moving it to LocalError_wfOpt*.log" << std::flush;
    	            }

                    if (CheckLogFile()) {
                        XTP_LOG(logDEBUG, *_pLog) << "CPMD: finished wavefunction optimization job. Continuing to projection onto AOs." << std::flush;
                    } else {
                        XTP_LOG(logDEBUG, *_pLog) << "CPMD: wavefunction optimization job failed" << std::flush;
                        return false;
                    }
                    
                } else {
                    XTP_LOG(logERROR, *_pLog) << _wfOpt_input_file_name << " failed to start. No shell accessible." << std::flush;
                    return false;
                }
                 
                _optWF = false;
                //continue as usual
            }
            
            //CPMD only needs to run once, or _optWF just finished running
            XTP_LOG(logDEBUG, *_pLog) << "CPMD: running [" << _executable << " " << _input_file_name << "]" << std::flush;

            if (std::system(NULL)) {
                std::string _command;
                _command = "cd " + _run_dir + "; rm -f LocalError*.log; " + _executable + " " + _input_file_name + ">" + _log_file_name;
                int check=std::system(_command.c_str());
                if (check==-1){
                    XTP_LOG(logERROR, *_pLog) << _input_file_name << " failed to start" << std::flush;
                    return false;
                }
                if (CheckLogFile()) {
                    XTP_LOG(logDEBUG, *_pLog) << "CPMD: finished job" << std::flush;
                    return true;
                } else {
                    XTP_LOG(logDEBUG, *_pLog) << "CPMD: job failed" << std::flush;
                }
            } else {
                XTP_LOG(logERROR, *_pLog) << _input_file_name << " failed to start" << std::flush;
                return false;
            }

            return true;

        }
               
        
        

        /**
         * Cleans up after the CPMD job.
         */
        void Cpmd::CleanUp() {
            
            // cleaning up the generated files
            if (_cleanup.size() != 0) {

                XTP_LOG(logDEBUG, *_pLog) << "Removing " << _cleanup << " files" << std::flush;
                tools::Tokenizer tok_cleanup(_cleanup, ", ");
                std::vector <std::string> cleanup_info;
                tok_cleanup.ToVector(cleanup_info);
                for (const std::string& substring:cleanup_info) {

                    if (substring == "inp") {
                        std::string file_name = _run_dir + "/" + _input_file_name;
                        remove(file_name.c_str());
                        if(_projectWF && _optWF)
                        {
                            std::string file_name = _run_dir + "/" + _wfOpt_input_file_name;
                            remove(file_name.c_str()); //also clean up the WF opt input
                        }
//                        if (_output_Vxc) {
//                            std::string file_name = _run_dir + "/" + _input_vxc_file_name;
//                            remove(file_name.c_str());
//                        }
                    }

                    if (substring == "sh") {
                        std::string file_name = _run_dir + "/" + _shell_file_name;
                        remove(file_name.c_str());
                    }

                    if (substring == "log") {
                        std::string file_name = _run_dir + "/" + _log_file_name;
                        remove(file_name.c_str());
                        if(_projectWF && _optWF)
                        {
                            remove(_wfOpt_log_file_name.c_str()); //also clean up the WF opt input
                        }
//                        if (_output_Vxc) {
//                            size_t lastdot = _log_file_name.find_last_of(".");
//                            if (lastdot == std::string::npos) {
//                                std::cerr << std::endl;
//                                std::cerr << "Could not remove Vxc log file" << std::flush;
//                            }
//                            std::string file_name2 = file_name.substr(0, lastdot) + "-2.log";
//                            remove(file_name2.c_str());
//                        }
                    }

                    if (substring == "chk") {
                        std::string file_name = _run_dir + "/" + "LATEST";
                        remove(file_name.c_str());
                        //if user sets custom execution options for CPMD,
                        //they may override this name and its not easy to track.
                        file_name = _run_dir + "/" + "RESTART.1";
                        remove(file_name.c_str());
                    }

                    if (substring == "matrices") {
                        std::string file_name = _run_dir + "/" + "OVERLAP";
                        remove(file_name.c_str());
                        file_name = _run_dir + "/" + "SPINDEN";
                        remove(file_name.c_str());
                        file_name = _run_dir + "/" + "WFNCOEF";
                        remove(file_name.c_str());
                        file_name = _run_dir + "/" + "CHOUT";
                        remove(file_name.c_str());
//                        if (_output_Vxc) {
//                            std::string file_name = _run_dir + "/" + "fort.24";
//                            remove(file_name.c_str());
//                        }
                    }

                    if (substring == "basis" && _write_basis_set) {
                        std::vector<std::string> fileswithfileending;
                        boost::filesystem::recursive_directory_iterator fit(_run_dir);
                        boost::filesystem::recursive_directory_iterator endit;
                        while (fit != endit) {
                            if (boost::filesystem::is_regular_file(* fit) &&
                                    fit->path().extension() == substring)
                            {
                                fileswithfileending.push_back(fit->path().filename().string());
                            }
                            ++fit;
                        }
                        for (const auto filename : fileswithfileending) {
                            std::string file_name = _run_dir + "/" + filename;
                            remove(file_name.c_str());
                        }
                    }
                    
                                        
                    if (substring == "density") {
                        std::string file_name = _run_dir + "/" + "DENSITY";
                        remove(file_name.c_str());
                    }
                    
                    if (substring == "elpot") {
                        std::string file_name = _run_dir + "/" + "ELPOT";
                        remove(file_name.c_str());
                    }

                }
            }
            return;
        }



        bool Cpmd::CheckLogFile() {

            // check if the log file exists
            boost::filesystem::path arg_path;

            std::string _full_name = (arg_path / _run_dir / _log_file_name).c_str();
            if(_optWF && _projectWF){ //CPMD needs to run twice; this is the _optWF run
                _full_name = (arg_path / _run_dir / _wfOpt_log_file_name).c_str();
            }
            std::ifstream _input_file(_full_name.c_str());

            if (_input_file.fail()) {
                XTP_LOG(logERROR, *_pLog) << "CPMD: " << _full_name << " is not found." << std::endl << std::flush;
                return false;
            };

            //Use brute force. Search every line for the termination string.
            //It doesn't appear at the very end, like in gaussian
            std::string::size_type self_energy_pos=std::string::npos;
            std::string _line;
            do {
                getline(_input_file, _line);
                self_energy_pos=_line.find("PROGRAM CPMD ENDED AT");
            } while (self_energy_pos==std::string::npos && !(_input_file.eof()));

            _input_file.close();

            if (self_energy_pos == std::string::npos) {
                XTP_LOG(logERROR, *_pLog) << "CPMD: " << _full_name << " is incomplete."<< std::endl << std::flush;
                return false;
            } else {
                XTP_LOG(logDEBUG,*_pLog) << "CPMD LOG is complete." <<std::endl << std::flush;
                return true;
            }
        }

        /**
         * Parses the CPMD Log file and stores data in the Orbitals object
         */
        bool Cpmd::ParseLogFile(Orbitals& orbitals) {
            std::string _line;
            std::vector<std::string> results;
            std::vector<Eigen::Vector3d> positions;
            
            XTP_LOG(logDEBUG, *_pLog) << "CPMD: parsing " << _log_file_name << std::flush;
            
            std::string _log_file_name_full = _log_file_name;
            if (_run_dir != "") _log_file_name_full = _run_dir + "/" + _log_file_name;
            
            // check if LOG file is complete
            if (!CheckLogFile()) return false;

            // save qmpackage name
            orbitals.setQMpackage("cpmd");
            orbitals.setDFTbasisName(_basisset_name);
            orbitals.setECPName(_ecp_name);
            
            
            std::ifstream _input_file(_log_file_name_full.c_str());
            while (_input_file) {
                
                getline(_input_file, _line);
                boost::trim(_line);
                
                
                /*
                 * number of electrons
                 */
                std::string::size_type electrons_pos = _line.find("alpha electrons");
                if (electrons_pos != std::string::npos) {
                    boost::algorithm::split(results, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    int number_of_electrons = (int) boost::lexical_cast<double>(results.back());
                    orbitals.setNumberOfAlphaElectrons(number_of_electrons);
                    XTP_LOG(logDEBUG, *_pLog) << "Alpha electrons: " << number_of_electrons << std::flush;
                }
                //TODO: beta electrons can be added later. Right now they are disabled by checking for the LSD keyword being passed to CPMD in custom options.
                
                /*
                 * atomic positions
                 */
                if (_line.find("*** ATOMS ***") != std::string::npos) {
                    XTP_LOG(logDEBUG, *_pLog) << "Wrapping atoms into periodic box of size "
                        <<_box[0]<<" "<<_box[1]<<" "<<_box[2]<<" Bohr."<<std::endl<<std::flush;
                    getline(_input_file, _line);
                    do{
                        getline(_input_file, _line);
                        boost::trim(_line);
                        if(_line.find("******") != std::string::npos)
                            break;

                        boost::algorithm::split(results, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                        tools::vec v;
                        v.setX(boost::lexical_cast<double>(results[2]));
                        v.setY(boost::lexical_cast<double>(results[3]));
                        v.setZ(boost::lexical_cast<double>(results[4]));
                        
                        //wrap into box
                        for(int k=0; k<3; k++)
                        {
                            v[k] = fmod(v[k], _box[k]);
                            if(v[k]<0){
                                v[k] = _box[k]+v[k];
                            }
                        }
//                        v=v*tools::conv::bohr2ang;
                        positions.push_back(Eigen::Vector3d(v.getX(),v.getY(),v.getZ())); //store positions in Bohr
                    }while(true);
                }
                
                
                /*
                 * Occupied/unoccupied states
                 */
                if (_line.find("NUMBER OF STATES:") != std::string::npos) {
                    boost::trim(_line);
                    boost::algorithm::split(results, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    int nstates = (boost::lexical_cast<int>(results[3]));
                    getline(_input_file, _line);
                    boost::trim(_line);
                    boost::algorithm::split(results, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    int nelectrons = (int) (boost::lexical_cast<double>(results[3]));
//                    orbitals.setNumberOfLevels(nelectrons/2, nstates-(nelectrons/2));
                    orbitals.setNumberOfOccupiedLevels(nelectrons/2);
                    XTP_LOG(logDEBUG, *_pLog) << "Occupied levels: " << nelectrons/2 << std::flush;
                    XTP_LOG(logDEBUG, *_pLog) << "Unoccupied levels: " << nstates-(nelectrons/2) << std::flush;
                }
                

                
                
            }
            XTP_LOG(logDEBUG, *_pLog) << "Done parsing" << std::flush;
            _input_file.close();
            
            //check that CPMD2TYPE_map is available
            if(CPMD2VOTCA_map==NULL){
                XTP_LOG(logDEBUG, *_pLog) << "CPMD: Can not convert atom order from CPMD to VOTCA." << std::flush;
                XTP_LOG(logDEBUG, *_pLog) << "CPMD: Please rerun with writing CPMD input (<tasks>input, parse</tasks>)." << std::flush;
                throw std::runtime_error("CPMD2TYPE_map unavailable, rerun with <tasks>input, parse</tasks>");
                exit(-1);
            }
            
            //store atoms to Orbitals in VOTCA's order
            bool has_atoms = orbitals.hasQMAtoms();
            for(unsigned int v=0; v<positions.size(); v++){
                int c = ConvAtomIndex_VOTCA2CPMD(v);
                if (has_atoms == false) {
                    //orbitals should not expose its atoms like this. Need a getter/setter, like we had last year.
                    orbitals.QMAtoms().push_back(QMAtom(v,CPMD2TYPE_map[c], positions[c]));
                } else {
                    QMAtom& pAtom = orbitals.QMAtoms().at(v);
                    pAtom.setPos(positions[c]);
                }
            }
            
            if(_projectWF){
                //MO coefficient and overlap matrices
                if(!loadMatrices(orbitals)) return false;

                //Check sanity of ECPs
                if(orbitals.hasQMAtoms()){
                    //lets test in CPMD's order

                    //iterate over elements
                    int i=0; //element index
                    int c=0; //atom index in CPMD
                    for(const std::string& element_name:_elements){
                        for(int a=0; a<_NA[i]; a++){ //check each atom
                            int v = ConvAtomIndex_CPMD2VOTCA(c);
                            QMAtom& pAtom = orbitals.QMAtoms().at(v);
                            if(pAtom.getNuccharge() != _ZV[i]){
                                XTP_LOG(logERROR, *_pLog) << "CPMD: ECP core charge mismatch for element "
                                        << element_name << ". CPMD uses a nuclear charge of " << _ZV[i]
                                        << "and your ECP file (" << _ecp_name << ") uses "
                                        << pAtom.getNuccharge()
                                        <<". Adjust your ECP file!" << std::flush;
                            }
                            c++;
                        }
                        i++;
                    }
                }
            }
            
            //fix order for version 5 of .orb files
            ReorderOutput(orbitals);            
            return true;

        }
        
        bool Cpmd::loadMatrices(Orbitals& orbitals)
        {
            int totAtoms=0;
            
            //check if WFNCOEF exists
            boost::filesystem::path arg_path;
            std::string _full_name = (arg_path / _run_dir / "WFNCOEF").c_str();
            std::ifstream wf_file(_full_name.c_str());
            if(wf_file.fail())
            {
                XTP_LOG(logERROR, *_pLog) << "CPMD: " << _full_name << " is not found." << std::endl << std::flush;
            }
            XTP_LOG(logDEBUG, *_pLog) << "CPMD: parsing " << _full_name << std::flush;
            
            //read WFNCOEF
            //variable names in all CAPS are variables from CPMD source code
            int count=0, endcount=0;
            wf_file.read((char*)&count, _F_int_size); //bytes in this record
            int bl=count;

            int NATTOT=0;
            wf_file.read((char*)&NATTOT, _F_int_size);    //number of basis functions
            bl-=_F_int_size;
            orbitals.setBasisSetSize(NATTOT);
            XTP_LOG(logDEBUG, *_pLog) << "Basis functions: " << NATTOT << std::flush;

            _NSP=0;                  //number of atom types
            wf_file.read((char*)&_NSP, _F_int_size);
            bl-=_F_int_size;

            //double ZV[NSP];             //core charge
            //int NA[NSP], NUMAOR[NSP];
            if(_ZV!=NULL)
            {
                delete[] _ZV;
                delete[] _NA;
                delete[] _NUMAOR;
            }
            _ZV=new double[_NSP];
            _NA=new int[_NSP];
            _NUMAOR=new int[_NSP];
            if(_ZV==NULL || _NA==NULL || _NUMAOR==NULL)
                throw std::runtime_error("Memory allocation failed");
            for(int i=0; i<_NSP; i++)
            {
                wf_file.read((char*)&_ZV[i], _F_real_size);    //core charge of atom type i
                wf_file.read((char*)&_NA[i], _F_int_size);     //number of atoms of type i
                wf_file.read((char*)&_NUMAOR[i], _F_int_size); //number of atomic orbitals of atom type i
                bl-=_F_real_size+_F_int_size*2;
                
                totAtoms+=_NA[i];
            }
            
            //check footer
            wf_file.read((char*)&endcount, _F_int_size);
            if(bl!=0 || endcount!=count){ //number of bytes read was wrong
                XTP_LOG(logERROR, *_pLog) << "CPMD: " << "could not parse record in "<< _full_name << std::flush;
                throw std::runtime_error("IO error");
                return false;
            }
                
            //NEW RECORD
            wf_file.read((char*)&count, _F_int_size); //bytes in this record
            bl=count;

            int NUMORB=count/_F_real_size/NATTOT;          //number of MOs (energy levels))
            XTP_LOG(logDEBUG, *_pLog) << "CPMD: number of energy levels: " << NUMORB << std::flush;
                //resize the coefficient matrix
            
            
            
            
            //map atomic orbitals to (CPMD) atom indeces so we can reorder the MO and Overlap matrices to VOTCA's atom order
            XTP_LOG(logDEBUG, *_pLog) << "CPMD: Reordering orbitals."<< std::flush;
            int* AO_CPMD2VOTCA_map= new int[NATTOT];
            int* VOTCA2numAOs= new int[totAtoms];   //number of AOs for each atom in VOTCA atomic order 
            int* VOTCA2firstAO= new int[totAtoms];  //index of first AO of each atom in VOTCA atomic order 
            if(AO_CPMD2VOTCA_map==NULL || VOTCA2numAOs==NULL || VOTCA2firstAO==NULL){
                throw std::runtime_error("Memory allocation failed");
            }
            int c=0;
            for(int i=0; i<_NSP; i++)
            {
                for(int a=0; a<_NA[i]; a++){
                    VOTCA2numAOs[ConvAtomIndex_CPMD2VOTCA(c)]=_NUMAOR[i];
                    c++; //increment CPMD atom index
                }
            }
                //find the beginning of each atom in the VOTCA order
            int orb=0;
            for(int v=0; v<totAtoms; v++){
                VOTCA2firstAO[v]=orb;
                orb+=VOTCA2numAOs[v];
            }
                //map each CMPD AO to a VOTCA AO
            c=0;    //CPMD atom index
            int co=0;
            int vo=0;
            for(int i=0; i<_NSP; i++)
            {
                for(int a=0; a<_NA[i]; a++){
                    for(int ofst=0; ofst<_NUMAOR[i]; ofst++){
                        int v = ConvAtomIndex_CPMD2VOTCA(c); //VOTCA atom index
                        vo = VOTCA2firstAO[v] + ofst; //VOTCA orbital index
                        AO_CPMD2VOTCA_map[co] = vo;
                        
                        co++;  //increment CPMD orbital index
                    }
                    c++; //increment CPMD atom index
                }
            }

            
            
            
            
            
            XTP_LOG(logDEBUG, *_pLog) << "CPMD: Reading MO coefficients."<< std::flush;
            Eigen::MatrixXd& mo_coefficients = orbitals.MOCoefficients();
            mo_coefficients.resize(NUMORB, NATTOT);
            
            //mo_coefficients need to be in VOTCA's atomic order
            //AO reordering comes later
            double XXMAT;
            for(int i=0; i<NUMORB; i++){
                for(int j=0; j<NATTOT; j++){  
                    wf_file.read((char*)&XXMAT, _F_real_size);
                    mo_coefficients(i,AO_CPMD2VOTCA_map[j])=XXMAT; //(MO, basisfunc in votca atom order)
                    bl-=_F_real_size;
                }
            }
            
            
            //check footer
            wf_file.read((char*)&endcount, _F_int_size);
            if(bl!=0 || endcount!=count){ //number of bytes read was wrong
                XTP_LOG(logERROR, *_pLog) << "CPMD: " << "could not parse record in "<< _full_name << std::endl << std::flush;
                throw std::runtime_error("IO error");
                return false;
            }
            
            wf_file.close();
            XTP_LOG(logDEBUG, *_pLog) << "CPMD: Done parsing" << std::flush;
            
            
            
            
            //check if OVERLAP exists
            _full_name = (arg_path / _run_dir / "OVERLAP").c_str();
            std::ifstream ov_file(_full_name.c_str());
            if(ov_file.fail())
            {
                XTP_LOG(logERROR, *_pLog) << "CPMD: " << _full_name << " is not found." << std::endl << std::flush;
            }
            XTP_LOG(logDEBUG, *_pLog) << "CPMD: parsing " << _full_name << std::flush;
            
            //read OVERLAP
            count=0, endcount=0;
            ov_file.read((char*)&count, 4); //bytes in this record
            bl=count;
            
            if(NATTOT*NATTOT!=count/8)
            {
                XTP_LOG(logERROR, *_pLog) << "CPMD: " << "Number of basis functions in the overlap and coefficient matrices do not match."<< std::endl << std::flush;
                throw std::runtime_error("IO error");
                return false;
            }
            
            
                //resize the overlap matrix
            Eigen::MatrixXd &overlap = orbitals.AOOverlap();
            overlap.resize(NATTOT,NATTOT);
            
            //read
            //Overlap need to be in VOTCA's atomic order
            XTP_LOG(logDEBUG, *_pLog) << "CPMD: Reading Overlap matrix."<< std::flush;
            double XSMAT;
            for(int i=0; i<NATTOT; i++){
                for(int j=0; j<NATTOT; j++){  
                    ov_file.read((char*)&XSMAT, _F_real_size);
                    overlap(AO_CPMD2VOTCA_map[i],AO_CPMD2VOTCA_map[j])=XSMAT;
                    bl-=_F_real_size;
                }
            }
            
            //check footer
            ov_file.read((char*)&endcount, _F_int_size);
            if(bl!=0 || endcount!=count){ //number of bytes read was wrong
                XTP_LOG(logERROR, *_pLog) << "CPMD: " << "could not parse footer in " <<
                        _full_name << ". Either the file is corrupted or " <<
                        "CPMD was compiled with non-standard int and real variable sizes. "<<
                        "Double check that Fortran_int_size and Fortran_real_size "<<
                        "in your .xml file match what CPMD was compiled with."<< std::flush;
                throw std::runtime_error("IO error");
                return false;
            }
            
            ov_file.close();
            XTP_LOG(logDEBUG, *_pLog) << "Done parsing" << std::flush;
            
            delete[] VOTCA2numAOs;
            delete[] VOTCA2firstAO;
            delete[] AO_CPMD2VOTCA_map;
            
            return true;
        }


        
        std::string Cpmd::FortranFormat(const double &number) {
            std::stringstream _ssnumber;
            if (number >= 0) _ssnumber << " ";
            _ssnumber << std::setiosflags(std::ios::fixed) << std::setprecision(8) << std::scientific << number;
            std::string _snumber = _ssnumber.str();
            boost::replace_first(_snumber, "e", "D");
            return _snumber;
        }


        void Cpmd::WriteChargeOption(){
            std::runtime_error("CPMD is a planewave code. It does not support background point charges.");
        }
        
        
    }
}
