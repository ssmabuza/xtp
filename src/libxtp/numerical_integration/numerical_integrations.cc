/* 
 *            Copyright 2009-2016 The VOTCA Development Team
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
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 * _olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
// Overload of uBLAS prod function with MKL/GSL implementations



#include <votca/xtp/numerical_integrations.h>
#include <boost/math/constants/constants.hpp>
#include <votca/xtp/radial_euler_maclaurin_rule.h>
#include <votca/xtp/sphere_lebedev_rule.h>
#include <votca/xtp/aoshell.h>
#include<votca/tools/constants.h>
#include <votca/xtp/elements.h>
#ifdef LIBXC
#include <xc.h>
#endif

#include <votca/xtp/exchange_correlation.h>
#include <fstream>
#include <boost/timer/timer.hpp>
#include <boost/algorithm/string.hpp>
#include <votca/xtp/vxc_functionals.h>
#include <iterator>
#include <string>
#include <math.h>
#include <bits/stl_vector.h>
// #include <xc.h>



namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        double NumericalIntegration::getExactExchange(const string _functional){
#ifdef LIBXC            
        
            double exactexchange=0.0;
            Vxc_Functionals map;
            std::vector<string> strs;
            
            boost::split(strs, _functional, boost::is_any_of(" "));
            if (strs.size()>2 ) {
                throw std::runtime_error("Too many functional names");
            }
            else if (strs.size()<1 ) {
                throw std::runtime_error("Specify at least one funcitonal");
            }
            
            for (unsigned i=0;i<strs.size();i++){
                //cout << strs[i] << endl;
                int func_id = map.getID(strs[i]); 
                if (func_id<0){
                    exactexchange=0.0;
                    break;
                }
                xc_func_type func;
                if (xc_func_init(&func, func_id, XC_UNPOLARIZED) != 0) {
                    fprintf(stderr, "Functional '%d' not found\n", func_id);
                    exit(1);
                }
                if (exactexchange>0 && func.cam_alpha>0){
                    throw std::runtime_error("You have specified two functionals with exact exchange");
                }
                exactexchange+=func.cam_alpha;
            
                
            
            }
            return exactexchange;
            
#else
            return 0.0;
#endif
            
        }
        
        
        ub::matrix<double> NumericalIntegration::IntegrateVXC_Atomblock(ub::matrix<double>& _density_matrix, AOBasis* basis,const string _functional){
            EXC = 0;
            // TODO: switch XC functionals implementation from LIBXC to base own calculation
            ExchangeCorrelation _xc;
            Vxc_Functionals map;
            std::vector<string> strs;           
            boost::split(strs, _functional, boost::is_any_of(" "));
            int xfunc_id = 0;
            
#ifdef LIBXC
            bool _use_votca = false;
            bool _use_separate = false;
            int cfunc_id = 0;

            if (strs.size() == 1) {
                xfunc_id = map.getID(strs[0]);
                if (xfunc_id < 0) _use_votca = true;
            }

            else if (strs.size() == 2) {
                cfunc_id = map.getID(strs[0]);
                xfunc_id = map.getID(strs[1]);
                _use_separate = true;
            }
            else {
                throw std::runtime_error("Please specify one combined or an exchange and a correlation functionals");

            }
            xc_func_type xfunc; // handle for exchange functional
            xc_func_type cfunc; // handle for correlation functional
            if (!_use_votca){
            if (xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED) != 0) {
                fprintf(stderr, "Functional '%d' not found\n", xfunc_id);
                exit(1);
            }
            
            xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED);
            if (xfunc.info->kind!=2 && !_use_separate){
                throw std::runtime_error("Your functional misses either correlation or exchange, please specify another functional, separated by whitespace");
            }
            
            if (_use_separate) {
                if (xc_func_init(&cfunc, cfunc_id, XC_UNPOLARIZED) != 0) {
                    fprintf(stderr, "Functional '%d' not found\n", cfunc_id);
                    exit(1);
                }
                xc_func_init(&cfunc, cfunc_id, XC_UNPOLARIZED);
                xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED);
                if ((xfunc.info->kind+cfunc.info->kind)!=1){
                    throw std::runtime_error("Your functionals are not one exchange and one correlation");
                }
            }
            }
#else
         if (strs.size() == 1) {
                xfunc_id = map.getID(strs[0]);
            }   
         else {
                throw std::runtime_error("Please specify one combined or an exchange and a correlation functionals");
         }
#endif

        
            
            //printf("The exchange functional '%s' is defined in the reference(s):\n%s\n", xfunc.info->name, xfunc.info->refs);
            //printf("The correlation functional '%s' is defined in the reference(s):\n%s\n", cfunc.info->name, cfunc.info->refs);

            // xc_func_end(&xfunc);

           
            // timers for testing
            /*
            boost::timer::cpu_timer cpu_t;
            cpu_t.start();
            double _t_AOvals = 0.0;
            double _t_rho = 0.0;
            double _t_grad_rho = 0.0;
            double _t_vxc =0.0;
            double _t_AOxc_rho=0.0;
            double _t_AOxc_grad=0.0;
            double _t_EXC1=0.0;
            double _t_EXC2=0.0;
            double _t_sum = 0.0;
            double _t_total = 0.0;
             boost::timer::cpu_times tenter = cpu_t.elapsed();
             */
            // generate a list of shells for each atom
            typedef std::vector< AOShell* >::iterator AOShellIterator;
            std::vector< std::vector< AOShellIterator > > _atomshells;
            std::vector< AOShellIterator > _singleatom;

            std::vector < int > _startIdx;
            std::vector < int > _blocksize;

            int _atomindex = 0;
            int _Idx       = 0;
            int _size      = 0;
            
            for (std::vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++) {
                
                
                if ( (*_row)->getIndex() == _atomindex ){
                    
                    _singleatom.push_back(_row);
                    _size += (*_row)->getNumFunc();
                    
                    
                } else {
                    
                    // append _singleatom to _atomshells
                    _atomshells.push_back(_singleatom);
                    _startIdx.push_back( _Idx );
                    _blocksize.push_back(_size);
                    // reset _singleatom
                    _singleatom.clear();
                    _size = (*_row)->getNumFunc();
                    _Idx       = (*_row)->getStartIndex();
                    _singleatom.push_back(_row);
                    _atomindex = (*_row)->getIndex();
                    
                }
                
                
            }
            
            _atomshells.push_back(_singleatom);
            _startIdx.push_back( _Idx );
                    _blocksize.push_back(_size);
           // cout << " Number of atoms " << _atomshells.size() << endl;
            /*
            for ( unsigned iatom = 0 ; iatom < _atomshells.size(); iatom++ ){
                cout << "atom " << iatom << " number of shells " << _atomshells[iatom].size() << " block start " << _startIdx[iatom] << " functions in atom " << _blocksize[iatom] << endl; 
            }
*/
            /*
            std::vector < ub::matrix_range< ub::matrix<double> > > _DMATblocks;
            // get stupid index magic vector of matrix_ranges per atom block
            for ( int rowatom = 0; rowatom < _atomshells.size(); rowatom++){
                for ( int colatom = 0 ; colatom <= rowatom; colatom++ ){
                    _DMATblocks.push_back(ub::subrange(_density_matrix,_startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom], _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom]));
                }
            }
            */
            
            
            // setup a list of min decay constants per atom
            // for every shell
            _atomindex = 0;
            double _decaymin = 1e7;
            std::vector< double > _minimal_decay;
            std::vector < vec > _positions;
            vec _localpos = (*basis->firstShell())->getPos();
            for ( std::vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++   ) {
                
                
                 if ( (*_row)->getIndex() == _atomindex ){
                     
                     // check all decay constants in this shell
                     for (AOShell::GaussianIterator itg = (*_row)->firstGaussian(); itg != (*_row)->lastGaussian(); itg++) {
                         AOGaussianPrimitive* gaussian = *itg;
                         double _decay = gaussian->decay;
                         if (_decay < _decaymin) {
                             _decaymin = _decay;
                         } // decay min check
                     
                     } // Gaussian Primitives 
                     
                 } else {  // if shell belongs to the actual atom
                     // add to mininal_decay vector
                     _minimal_decay.push_back(_decaymin);
                     _positions.push_back( _localpos );
                     // reset counters
                     _decaymin = 1e7;
                     _localpos = (*_row)->getPos();

                     _atomindex++;
                     
                     // check all decay constants in this shell
                     for (AOShell::GaussianIterator itg = (*_row)->firstGaussian(); itg != (*_row)->lastGaussian(); itg++) {
                         AOGaussianPrimitive* gaussian = *itg;
                         double _decay = gaussian->decay;
                         if (_decay < _decaymin) {
                             _decaymin = _decay;
                         } // decay min check
                     
                     } // Gaussian Primitives 
                     
                 
                 }
            } // all shells
                 
            // push final atom
            _minimal_decay.push_back(_decaymin);
            _positions.push_back( _localpos );
             
            /* for ( int i =0; i < _minimal_decay.size(); i++){
                 
                 cout << "Atom " << i << " min decay " << _minimal_decay[i] <<  " at " << _positions[i] << endl; 
                 
             } */
             
             
             // for each gridpoint, check the value of exp(-a*(r-R)^2) < 1e-10
             //                             = alpha*(r-R)^2 >~ 20.7
            
            std::vector< std::vector< std::vector<int> > > _significant_atoms;
            
            // each atomic grid
            for (unsigned i = 0; i < _grid.size(); i++) {
            
                std::vector< std::vector<int> > _significant_atoms_atomgrid;
                
                // each point of the atomic grid
                for (unsigned j = 0; j < _grid[i].size(); j++) {

                    std::vector<int> _significant_atoms_gridpoint;
                    vec grid;
                    grid.setX( _grid[i][j].grid_x);
                    grid.setY( _grid[i][j].grid_y);
                    grid.setZ( _grid[i][j].grid_z);

                    
                    
                    // check all atoms
                    for ( unsigned iatom = 0 ; iatom < _minimal_decay.size(); iatom++){

                        vec dist = grid - _positions[iatom];
                        double distsq = dist.getX()*dist.getX() + dist.getY()*dist.getY()  + dist.getZ()*dist.getZ() ;
                        
                        // if contribution is smaller than -ln(1e-10), add atom to list
                        if ( (_minimal_decay[iatom] * distsq) < 20.7 ){
                            _significant_atoms_gridpoint.push_back(iatom);
                        }
                        
                    } // check all atoms

                    _significant_atoms_atomgrid.push_back(  _significant_atoms_gridpoint );

                    
                } // all points of this atom grid
                
                _significant_atoms.push_back(_significant_atoms_atomgrid);
               
            } // atomic grids
   
            
            /*
            for ( int i = 0; i < _significant_atoms.size(); i++ ){
                
                cout << " Atomgrid: " << i << endl;
                 for ( int j = 0; j < _significant_atoms[i].size(); j++ ){
                     
                     cout << " Atom: " << i << " gridpoint " << j << " significant atoms " << _significant_atoms[i][j].size() << endl;
                     
                 }
                
            }*/
            
            
            
   
             //exit(0);
            
            
            int total_grid =0;
            int significant_grid = 0;
            for ( unsigned i = 0; i < _significant_atoms.size(); i++ ){
                
                total_grid += _grid[i].size(); 
                
                for ( unsigned j = 0; j < _significant_atoms[i].size(); j++ ){
                    
                    int gridpointsize = _significant_atoms[i][j].size();
                    significant_grid += gridpointsize*(gridpointsize+1);
                    
                }
                
                
                
            }
            int natoms = _grid.size();
            
            total_grid = total_grid * ( natoms*(natoms+1) ) / 2;
            
            // cout << "Total number of atom blocks " << total_grid << " after removal of insignificant blocks " << significant_grid/2 << endl;
            //cout << "# Atom blocks by removal of insignificant blocks reduced to " << 50*double(significant_grid)/double(total_grid) << "%" << endl;
            
            ub::matrix<double> XCMAT = ub::zero_matrix<double>(basis->_AOBasisSize, basis->_AOBasisSize);
            
            // parallelization: distribute over threads inside one atom
            int nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif
            
               
               
            // separate storage for each thread
            std::vector< ub::matrix<double> > XCMAT_thread;
            std::vector<double> EXC_thread;
            for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                
                XCMAT_thread.push_back( ub::zero_matrix<double>(basis->_AOBasisSize, basis->_AOBasisSize) );
                EXC_thread.push_back(0.0);
            }
            
            
            
            // for every atom
            for (unsigned i = 0; i < _grid.size(); i++) {
	      // for each point in atom grid
                
                // number of points in this atomgrid
                int atom_points = _grid[i].size();
                // divide among threads
                int atom_points_per_thread = atom_points/nthreads;
                std::vector<int> _thread_start;
                std::vector<int> _thread_stop;
                for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                    _thread_start.push_back( i_thread * atom_points_per_thread );
                    _thread_stop.push_back( (i_thread + 1) * atom_points_per_thread );
                }
                // final stop must be size
                _thread_stop[nthreads-1] = atom_points;
                
               /*
                
               for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                    
                    cout << "Thread " << i_thread << " start " << _thread_start[i_thread] << " stop " << _thread_stop[i_thread] << endl; 
                    
                }
                */ 
                #pragma omp parallel for
                for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                for (int j = _thread_start[i_thread]; j < _thread_stop[i_thread]; j++) {


                 //  boost::timer::cpu_times t0 = cpu_t.elapsed();


                    // get value of orbitals at each gridpoint (vector as 1D boost matrix object -> prod )
                    //ub::matrix<double> AOgrid = ub::zero_matrix<double>(basis->_AOBasisSize, 1);

                   ub::matrix<double> AOgrid = ub::zero_matrix<double>(1, basis->_AOBasisSize); // TRY MORE USEFUL DATA
                   
                   
		    // get value of density gradient at each gridpoint
                    //ub::matrix<double> gradAOgrid = ub::zero_matrix<double>(basis->_AOBasisSize, 3); // for Gradients of AOs

                   ub::matrix<double> gradAOgrid = ub::zero_matrix<double>(3, basis->_AOBasisSize); // for Gradients of AOs
                    
                   ub::matrix<double>  rho_mat = ub::zero_matrix<double>(1,1);
                   //ub::matrix<double> grad_rho = ub::zero_matrix<double>(3,1);
                    
                   ub::matrix<double> grad_rho = ub::zero_matrix<double>(1,3);
		    // evaluate AO Functions for all shells, NOW BLOCKWISE

                    // for each significant atom for this grid point
                    for ( unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size() ; sigrow++){
                    
                        // this atom
                        int rowatom = _significant_atoms[i][j][sigrow];
                        
                    
                    
                        // for each shell in this atom
                        for ( unsigned ishell = 0 ; ishell < _atomshells[rowatom].size() ; ishell++ ){
                         //   boost::timer::cpu_times tstartshells = cpu_t.elapsed();
                            AOShellIterator _row = _atomshells[rowatom][ishell];
                            // for density, fill sub-part of AOatgrid
                            //ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), 0, 1);
                            ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, 0, 1, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc());
                            // (*_row)->EvalAOspace(_AOgridsub, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);

                            // gradient of density
                            //ub::matrix_range< ub::matrix<double> > _gradAO = ub::subrange(gradAOgrid, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), 0, 3);
                            ub::matrix_range< ub::matrix<double> > _gradAO = ub::subrange(gradAOgrid, 0, 3, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc());
                            //(*_row)->EvalAOGradspace(_gradAO, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);
                            (*_row)->EvalAOspace(_AOgridsub, _gradAO , _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);
                         //   boost::timer::cpu_times tendshells = cpu_t.elapsed();
                            
                            // _t_AOvals +=  (tendshells.wall-tstartshells.wall)/1e9;

                        }  // shell in atom
                        
                        /* ub::matrix<double> _temp     = ub::zero_matrix<double>(_blocksize[rowatom],1);
                        ub::matrix<double> _tempgrad = ub::zero_matrix<double>(_blocksize[rowatom],3);
                        
                        ub::matrix_range< ub::matrix<double> > _AOgridrow     = ub::subrange(    AOgrid, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom], 0, 1);

                         * 
                         */
                        
                        ub::matrix<double> _temp     = ub::zero_matrix<double>(1,_blocksize[rowatom]);
                        ub::matrix<double> _tempgrad = ub::zero_matrix<double>(3,_blocksize[rowatom]);
                        
                        ub::matrix_range< ub::matrix<double> > _AOgridrow     = ub::subrange(    AOgrid, 0,1, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);


                        // for each atom
                        // for all significant atoms of triangular matrix
                        for ( unsigned sigcol = 0; sigcol < _significant_atoms[i][j].size() ; sigcol++){
                            int colatom = _significant_atoms[i][j][sigcol];
                            if ( colatom > rowatom ) break;
                            
                            // get the already calculated AO values
                            //ub::matrix_range< ub::matrix<double> >     _AOgridcol = ub::subrange(    AOgrid, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom], 0, 1);
                            //ub::matrix_range< ub::matrix<double> > _gradAOgridcol = ub::subrange(gradAOgrid, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom], 0, 3);

                            ub::matrix_range< ub::matrix<double> >     _AOgridcol = ub::subrange(    AOgrid, 0, 1, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom]);
                            ub::matrix_range< ub::matrix<double> > _gradAOgridcol = ub::subrange(gradAOgrid, 0, 3, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom]);


                            
                            //ub::matrix_range< ub::matrix<double> > DMAT_here = ub::subrange( _density_matrix, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom], _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom]);
                            ub::matrix_range< ub::matrix<double> > DMAT_here = ub::subrange( _density_matrix, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom], _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);
                            
                            // update _temp, careful if diagonal!
                         /*   if ( colatom == rowatom ){
                                _temp     += 0.5 * ub::prod( DMAT_here,_AOgridcol);
                                _tempgrad += 0.5 * ub::prod( DMAT_here,_gradAOgridcol);
                            } else {
                                
                                _temp     += ub::prod( DMAT_here,    _AOgridcol);
                                _tempgrad += ub::prod( DMAT_here,_gradAOgridcol);

                            } */
                            
                            
                            
                             if ( colatom == rowatom ){
                                _temp     += 0.5 * ub::prod( _AOgridcol, DMAT_here);
                                _tempgrad += 0.5 * ub::prod( _gradAOgridcol, DMAT_here);
                            } else {
                                
                                _temp     += ub::prod(  _AOgridcol, DMAT_here);
                                _tempgrad += ub::prod( _gradAOgridcol, DMAT_here);

                            }

                        } //col shells
                        
                        
                        //ub::matrix_range< ub::matrix<double> > _gradAOgridrow = ub::subrange(gradAOgrid, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom], 0, 3);
                        ub::matrix_range< ub::matrix<double> > _gradAOgridrow = ub::subrange(gradAOgrid, 0,3, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);
                        
                        //rho_mat  += ub::prod(ub::trans(    _AOgridrow),_temp);
                        //grad_rho += ub::prod(ub::trans(_gradAOgridrow),_temp) +  ub::prod(ub::trans(_tempgrad),_AOgridrow) ;

                        rho_mat  += ub::prod(_temp, ub::trans( _AOgridrow) );
                        grad_rho += ub::prod(_temp, ub::trans(_gradAOgridrow)) +  ub::prod(_AOgridrow,ub::trans(_tempgrad)) ;

                        
                        
                        
                    } // row shells 


                    double rho      = 2.0 * rho_mat(0,0);
                 //   boost::timer::cpu_times t3 = cpu_t.elapsed();
                    // +=  (t3.wall-t0.wall)/1e9;
                    

		    if ( rho < 1.e-15 ) continue; // skip the rest, if density is very small
                    grad_rho = 2.0 * grad_rho;
                                     

               //     cout << " size1 grad" << grad_rho.size1() << " size2 " << grad_rho.size2() << endl;
                //    cout << "rho ABLK" << rho << " " << grad_rho(0,0) << " "  << grad_rho(0,1) << " "  << grad_rho(0,2)  << endl;
                 //   exit(0);
                    
                    
                    // get XC for this density_at_grid
                    double f_xc;      // E_xc[n] = int{n(r)*eps_xc[n(r)] d3r} = int{ f_xc(r) d3r }
                    double df_drho;   // v_xc_rho(r) = df/drho
                    double df_dsigma; // df/dsigma ( df/dgrad(rho) = df/dsigma * dsigma/dgrad(rho) = df/dsigma * 2*grad(rho))

 #ifdef LIBXC                   
                    if (_use_votca) {
#endif
                        _xc.getXC(xfunc_id, rho, grad_rho(0, 0), grad_rho(0, 1), grad_rho(0, 2), f_xc, df_drho, df_dsigma);
#ifdef LIBXC
                    }                        // evaluate via LIBXC, if compiled, otherwise, go via own implementation

                    else {
                        //double sigma = ub::prod(ub::trans(grad_rho),grad_rho)(0,0);

                        double sigma = ub::prod(grad_rho, ub::trans(grad_rho))(0, 0);

                        double exc[1];
                        double vsigma[1]; // libxc 
                        double vrho[1]; // libxc df/drho
                        switch (xfunc.info->family) {
                            case XC_FAMILY_LDA:
                                xc_lda_exc_vxc(&xfunc, 1, &rho, exc, vrho);
                                break;
                            case XC_FAMILY_GGA:
                            case XC_FAMILY_HYB_GGA:
                                xc_gga_exc_vxc(&xfunc, 1, &rho, &sigma, exc, vrho, vsigma);
                                break;
                        }
                        f_xc = exc[0];
                        df_drho = vrho[0];
                        df_dsigma = vsigma[0];
                        if (_use_separate) {
                            // via libxc correlation part only
                            switch (cfunc.info->family) {
                                case XC_FAMILY_LDA:
                                    xc_lda_exc_vxc(&cfunc, 1, &rho, exc, vrho);
                                    break;
                                case XC_FAMILY_GGA:
                                case XC_FAMILY_HYB_GGA:
                                    xc_gga_exc_vxc(&cfunc, 1, &rho, &sigma, exc, vrho, vsigma);
                                    break;
                            }

                            f_xc += exc[0];
                            df_drho += vrho[0];
                            df_dsigma += vsigma[0];
                        }
                    }
#endif
                    
                    


                 //   boost::timer::cpu_times t4 = cpu_t.elapsed();
                   // _t_vxc += (t4.wall-t3.wall)/1e9;
                    
		    // density part
		    //ub::matrix<double> _addXC = _grid[i][j].grid_weight * df_drho * AOgrid;
                    ub::matrix<double> _addXC = _grid[i][j].grid_weight * df_drho * AOgrid *0.5;
                //    boost::timer::cpu_times t5 = cpu_t.elapsed();
                   // _t_AOxc_rho += (t5.wall-t4.wall)/1e9;
                    
                    
		    // gradient part
                    //_addXC+=  4.0*df_dsigma * _grid[i][j].grid_weight * ub::prod(gradAOgrid,grad_rho);
                    //_addXC+=  4.0*df_dsigma * _grid[i][j].grid_weight * ub::prod(grad_rho,gradAOgrid);
                    _addXC+=  2.0*df_dsigma * _grid[i][j].grid_weight * ub::prod(grad_rho,gradAOgrid);
                //    boost::timer::cpu_times t6 = cpu_t.elapsed();
                 //   _t_AOxc_grad += (t6.wall-t5.wall)/1e9;

		    // finally combine (super-slow...)
                    // XCMAT +=  ub::prod(_addXC,ub::trans(AOgrid));

                    
                    // Exchange correlation energy
                    EXC_thread[i_thread] += _grid[i][j].grid_weight * rho * f_xc;
                 //   boost::timer::cpu_times t6a = cpu_t.elapsed();
                 //   _t_EXC1 += (t6a.wall-t6.wall)/1e9;
                    
                    // combine/sum atom-block wise, only trigonal part, symmetrize later
                    // for each significant atom for this grid point
                    // parallelization only accesses atomblock information (_addXC, AOgrid -> XCmatblock), so no trouble with shared memory access )
                    // #pragma omp parallel for
                    for (unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size(); sigrow++) {
                        
                        // this atom
                        int rowatom = _significant_atoms[i][j][sigrow];
                        // line reference of _addXC
                        //ub::matrix_range< ub::matrix<double> > _rowXC = ub::subrange( _addXC, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom], 0, 1);
                        ub::matrix_range< ub::matrix<double> > _rowXC = ub::subrange( _addXC, 0 , 1, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);    

                        for (unsigned sigcol = 0; sigcol < _significant_atoms[i][j].size(); sigcol++) {
                            int colatom = _significant_atoms[i][j][sigcol];
                            // if (colatom > rowatom) break;

                            // line reference of AOgrid 
                            //ub::matrix_range< ub::matrix<double> > _AOcol = ub::subrange( AOgrid,  _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom], 0, 1);
                            ub::matrix_range< ub::matrix<double> > _AOcol = ub::subrange( AOgrid, 0,1,  _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom]);
                            
                            // update block reference of XCMAT
                            ub::matrix_range<ub::matrix<double> > _XCmatblock = ub::subrange( XCMAT_thread[i_thread],_startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom], _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom] );
                            //_XCmatblock += ub::prod( _rowXC, ub::trans(_AOcol)  );
                            _XCmatblock += ub::prod( ub::trans(_rowXC), _AOcol  );
                            
                            
                            // update the other block
                            
                            
                            
                            
                        } // significant col
                    } // significant row 
                    
                    
                    
                   // boost::timer::cpu_times t7 = cpu_t.elapsed();
                 //   _t_sum += (t7.wall-t6a.wall)/1e9;
                    

                } // j: for each point in atom grid
                }// each thread
            } // i: for each atom grid

            
            // sum thread matrices
            for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                EXC += EXC_thread[i_thread];
                #pragma omp parallel for
                for (unsigned _i = 0; _i < XCMAT.size1(); _i++) {
                    //for (int _j = 0; _j <= _i; _j++) {
                        for (unsigned _j = 0; _j <XCMAT.size2(); _j++) {
                    XCMAT( _i, _j ) += XCMAT_thread[i_thread](_i, _j);
                    }
                }
            }
            
            
            
            XCMAT += ub::trans(XCMAT);
            
            // symmetrize 
            //#pragma omp parallel for
            //for (int _i = 0; _i < XCMAT.size1(); _i++) {
            //    for (int _j = 0; _j < _i; _j++) {
            //        XCMAT(_j, _i) = XCMAT(_i, _j);
            //    }
           // }


            
            
          //  boost::timer::cpu_times t8 = cpu_t.elapsed();
                    
            const ub::vector<double> DMAT_array = _density_matrix.data();
            const ub::vector<double> XCMAT_array = XCMAT.data();
            
            for ( unsigned i = 0; i < DMAT_array.size(); i++ ){
                EXC -= DMAT_array[i] * XCMAT_array[i];
            }

         //   boost::timer::cpu_times t9 = cpu_t.elapsed();
         //   _t_EXC2 += (t9.wall-t8.wall)/1e9;
                    
            
/*
            cout << " ATBLOCK Time AOVals      : " << _t_AOvals << endl;
            cout << " ATBLOCK Time rho         : " << _t_rho << endl;
            cout << " ATBLOCK Time grad rho    : " << _t_grad_rho << endl;
            cout << " ATBLOCK Time Vxc         : " << _t_vxc << endl;
            cout << " ATBLOCK Time AOxc rho    : " << _t_AOxc_rho << endl;
            cout << " ATBLOCK Time AOxc grad   : " << _t_AOxc_grad << endl;
            cout << " ATBLOCK Time AOxc sum    : " << _t_sum << endl;
            cout << " ATBLOCK Time Exc1        : " << _t_EXC1 << endl;
            cout << " ATBLOCK Time Exc2        : " << _t_EXC2 << endl;
 */           
                        // boost::timer::cpu_times texit = cpu_t.elapsed();
                              //  _t_total = (texit.wall-tenter.wall)/1e9;
                    
                                 
                                 
       //     cout << " ATBLOCK TOTAL            : " << _t_total << endl;
            
         
            
            
            
            return XCMAT;

        }
        
        // not used anymore as atomblock is faster
        ub::matrix<double> NumericalIntegration::IntegrateVXC_block(ub::matrix<double>& _density_matrix, AOBasis* basis){
            EXC=0;
            // TODO: switch XC functionals implementation from LIBXC to base own calculation
            #ifdef LIBXC
            xc_func_type xfunc; // handle for exchange functional
            xc_func_type cfunc; // handle for correlation functional

            // define PBE here (should be some optional setting))
            int xfunc_id = XC_GGA_X_PBE;
            int cfunc_id = XC_GGA_C_PBE;
            
            if(xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED) != 0){
               fprintf(stderr, "Functional '%d' not found\n", xfunc_id);
               exit(1);
            }
            
            if(xc_func_init(&cfunc, cfunc_id, XC_UNPOLARIZED) != 0){
               fprintf(stderr, "Functional '%d' not found\n", cfunc_id);
               exit(1);
            }
#else
             ExchangeCorrelation _xc;
#endif
            
            //printf("The exchange functional '%s' is defined in the reference(s):\n%s\n", xfunc.info->name, xfunc.info->refs);
            //printf("The correlation functional '%s' is defined in the reference(s):\n%s\n", cfunc.info->name, cfunc.info->refs);

            // xc_func_end(&xfunc);

           
            // timers for testing
            //boost::timer::cpu_timer cpu_t;
            //cpu_t.start();

            double _t_rho = 0.0;
            double _t_grad_rho = 0.0;
            double _t_vxc =0.0;
            double _t_AOxc_rho=0.0;
            double _t_AOxc_grad=0.0;
            double _t_sum = 0.0;
            double _t_total = 0.0;
            //boost::timer::cpu_times tenter = cpu_t.elapsed();
            
            
            vector < ub::matrix_range< ub::matrix<double> > > _DMATblocks;
            // get stupid index magic vector of matrix_ranges per shell
            for ( vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++ ){
                for ( vector< AOShell* >::iterator _col = basis->firstShell() ; _col != _row; _col++){
                    _DMATblocks.push_back(ub::subrange(_density_matrix,(*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), (*_col)->getStartIndex(), (*_col)->getStartIndex()+(*_col)->getNumFunc()));
                }
                // push same iterator value
                _DMATblocks.push_back(ub::subrange(_density_matrix,(*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc()));
            }
            
            ub::matrix<double> XCMAT = ub::zero_matrix<double>(basis->_AOBasisSize, basis->_AOBasisSize);
            
            // for every atom
            for (unsigned i = 0; i < _grid.size(); i++) {
	      // for each point in atom grid
                for (unsigned j = 0; j < _grid[i].size(); j++) {


                   // boost::timer::cpu_times t0 = cpu_t.elapsed();


                    // get value of orbitals at each gridpoint (vector as 1D boost matrix object -> prod )
                    ub::matrix<double> AOgrid = ub::zero_matrix<double>(basis->_AOBasisSize, 1);

		    // get value of density gradient at each gridpoint
                    ub::matrix<double> gradAOgrid = ub::zero_matrix<double>(basis->_AOBasisSize, 3); // for Gradients of AOs

                    ub::matrix<double> rho_mat  = ub::zero_matrix<double>(1,1); // matrix to use ub::prod
                    ub::matrix<double> grad_rho = ub::zero_matrix<double>(3,1);
                    
		    // evaluate AO Functions for all shells, NOW BLOCKWISE
                      int iterblock = 0;
                    for (vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++) {

                        // for density, fill sub-part of AOatgrid
                        ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), 0, 1);
                        (*_row)->EvalAOspace(_AOgridsub, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);

                        // for density gradient, fill sub-part of gradAOgrid
                        ub::matrix_range< ub::matrix<double> > _gradAO = ub::subrange(gradAOgrid, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), 0, 3);
                        (*_row)->EvalAOGradspace(_gradAO, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);

                        ub::matrix<double> _temp = ub::zero_matrix<double>((*_row)->getNumFunc(),1);
                        ub::matrix<double> _tempgrad = ub::zero_matrix<double>((*_row)->getNumFunc(),3);

                        // now do the same for _col shell, only lower triangle blocks, we can access 
                        for (vector< AOShell* >::iterator _col = basis->firstShell() ; _col != _row; _col++) {
           
                            // for density, get sub-part of AOatgrid
                            ub::matrix_range< ub::matrix<double> > _AOgridcol = ub::subrange(AOgrid, (*_col)->getStartIndex(), (*_col)->getStartIndex()+(*_col)->getNumFunc(), 0, 1);
                            // for gradient, get sub-part of gradAOgrid
                            ub::matrix_range< ub::matrix<double> > _gradAOgridcol = ub::subrange(gradAOgrid,(*_col)->getStartIndex(), (*_col)->getStartIndex()+(*_col)->getNumFunc() , 0, 3);
                            
                            // ub::matrix_range< ub::matrix<double> > _DMAThere = ub::subrange(_density_matrix,(*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), (*_col)->getStartIndex(), (*_col)->getStartIndex()+(*_col)->getNumFunc());
                            
                            
                            // update _temp, careful if diagonal!
                            _temp     += ub::prod(_DMATblocks[iterblock],    _AOgridcol);
                            _tempgrad += ub::prod(_DMATblocks[iterblock],_gradAOgridcol);
                            iterblock++; 
                            
 
                        } //col shells
                          
                        // "diagonal shell"   needs to be added
                        _temp     += 0.5 * ub::prod(_DMATblocks[iterblock],_AOgridsub);
                        _tempgrad += 0.5 * ub::prod(_DMATblocks[iterblock],_gradAO);
                        iterblock++; 
                        
                        rho_mat  += ub::prod(ub::trans(_AOgridsub),_temp);
                        grad_rho += ub::prod(ub::trans(_gradAO),_temp) +  ub::prod(ub::trans(_tempgrad),_AOgridsub) ;
                        
                        
                        
                    } // row shells

                    double rho = 2.0 * rho_mat(0,0);

                    //boost::timer::cpu_times t3 = cpu_t.elapsed();

                    //_t_rho += (t3.wall-t0.wall)/1e9;
                    
		    if ( rho < 1.e-15 ) continue; // skip the rest, if density is very small
                    grad_rho = 2.0 * grad_rho;
                    
                    
                    
                    
                    
                    // get XC for this density_at_grid
                    double f_xc;      // E_xc[n] = int{n(r)*eps_xc[n(r)] d3r} = int{ f_xc(r) d3r }
                    double df_drho;   // v_xc_rho(r) = df/drho
                    double df_dsigma; // df/dsigma ( df/dgrad(rho) = df/dsigma * dsigma/dgrad(rho) = df/dsigma * 2*grad(rho))




                    // evaluate via LIBXC, if compiled, otherwise, go via own implementation
#ifdef LIBXC
		    double sigma = ub::prod(ub::trans(grad_rho),grad_rho)(0,0);

                    double exc[1];
                    double vsigma[1]; // libxc 
                    double vrho[1]; // libxc df/drho
                    switch (xfunc.info->family) {
                        case XC_FAMILY_LDA:
                            xc_lda_exc_vxc(&xfunc, 1, &rho, exc, vrho);
                            break;
                        case XC_FAMILY_GGA:
                        case XC_FAMILY_HYB_GGA:
                            xc_gga_exc_vxc(&xfunc, 1, &rho, &sigma, exc, vrho, vsigma);
                            break;
                    }
                    f_xc      = exc[0];
                    df_drho   = vrho[0];
                    df_dsigma = vsigma[0];

                    // via libxc correlation part only
                    switch (cfunc.info->family) {
                        case XC_FAMILY_LDA:
                            xc_lda_exc_vxc(&cfunc, 1, &rho, exc, vrho);
                            break;
                        case XC_FAMILY_GGA:
                        case XC_FAMILY_HYB_GGA:
                            xc_gga_exc_vxc(&cfunc, 1, &rho, &sigma, exc, vrho, vsigma);
                            break;
                    }

                    f_xc      += exc[0];
                    df_drho   += vrho[0];
                    df_dsigma += vsigma[0];

#else
                    _xc.getXC(-2, rho, grad_rho(0,0), grad_rho(1,0), grad_rho(2,0), f_xc, df_drho, df_dsigma);
#endif

                   // boost::timer::cpu_times t4 = cpu_t.elapsed();
                   // _t_vxc += (t4.wall-t3.wall)/1e9;
                    
		    // density part
		    ub::matrix<double> _addXC = _grid[i][j].grid_weight * df_drho * AOgrid;
                    EXC += _grid[i][j].grid_weight * rho * f_xc;
                    //boost::timer::cpu_times t5 = cpu_t.elapsed();
                   // _t_AOxc_rho += (t5.wall-t4.wall)/1e9;
                    
                    
		    // gradient part
		    int size = gradAOgrid.size1();
                    ub::matrix_range< ub::matrix<double> > _gradAO_x = ub::subrange(gradAOgrid, 0, size , 0, 1);
                    ub::matrix_range< ub::matrix<double> > _gradAO_y = ub::subrange(gradAOgrid, 0, size , 1, 2);
                    ub::matrix_range< ub::matrix<double> > _gradAO_z = ub::subrange(gradAOgrid, 0, size , 2, 3);
		    _addXC += 4.0*df_dsigma * _grid[i][j].grid_weight *(grad_rho(0,0) * _gradAO_x + grad_rho(1,0) * _gradAO_y + grad_rho(2,0) * _gradAO_z );

                    //boost::timer::cpu_times t6 = cpu_t.elapsed();
                    //_t_AOxc_grad += (t6.wall-t5.wall)/1e9;

		    // finally combine
		    XCMAT +=  ub::prod(_addXC,ub::trans(AOgrid));

                    //boost::timer::cpu_times t7 = cpu_t.elapsed();
                    //_t_sum += (t7.wall-t6.wall)/1e9;

                } // j: for each point in atom grid
            } // i: for each atom grid

            
            const ub::vector<double> DMAT_array = _density_matrix.data();
            const ub::vector<double> XCMAT_array = XCMAT.data();
            
            for ( unsigned i = 0; i < DMAT_array.size(); i++ ){
            
                EXC -= DMAT_array[i] * XCMAT_array[i];
                
            }
            
            
            
            

            cout << " BLOCK Time rho         : " << _t_rho << endl;
            cout << " BLOCK Time grad rho    : " << _t_grad_rho << endl;
            cout << " BLOCK Time Vxc         : " << _t_vxc << endl;
            cout << " BLOCK Time AOxc rho    : " << _t_AOxc_rho << endl;
            cout << " BLOCK Time AOxc grad   : " << _t_AOxc_grad << endl;
            cout << " BLOCK Time AOxc sum    : " << _t_sum << endl;


                    
                    
                                 
                                 
            cout << " BLOCK TOTAL            : " << _t_total << endl;
            
            return XCMAT;

        }
                
        // numerically integrate the elements of the AOXC matrix
        ub::matrix<double> NumericalIntegration::IntegrateVXC(ub::matrix<double>& _density_matrix, AOBasis* basis){
                EXC=0;    
            // TODO: switch XC functionals implementation from LIBXC to base own calculation
            #ifdef LIBXC
            xc_func_type xfunc; // handle for exchange functional
            xc_func_type cfunc; // handle for correlation functional

            // define PBE here (should be some optional setting))
            int xfunc_id = XC_GGA_X_PBE;
            int cfunc_id = XC_GGA_C_PBE;
            
            if(xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED) != 0){
               fprintf(stderr, "Functional '%d' not found\n", xfunc_id);
               exit(1);
            }
            
            if(xc_func_init(&cfunc, cfunc_id, XC_UNPOLARIZED) != 0){
               fprintf(stderr, "Functional '%d' not found\n", cfunc_id);
               exit(1);
            }
#else
             ExchangeCorrelation _xc;
#endif
            
            //printf("The exchange functional '%s' is defined in the reference(s):\n%s\n", xfunc.info->name, xfunc.info->refs);
            //printf("The correlation functional '%s' is defined in the reference(s):\n%s\n", cfunc.info->name, cfunc.info->refs);

            // xc_func_end(&xfunc);

           
            // timers for testing
            boost::timer::cpu_timer cpu_t;
            cpu_t.start();

            double _t_AOgrid = 0.0;
            double _t_rho = 0.0;
            double _t_grad_rho = 0.0;
            double _t_vxc =0.0;
            double _t_AOxc_rho=0.0;
            double _t_AOxc_grad=0.0;
            double _t_sum = 0.0;
            double _t_total=0.0;
                                boost::timer::cpu_times tenter = cpu_t.elapsed();
            ub::matrix<double> XCMAT = ub::zero_matrix<double>(basis->_AOBasisSize, basis->_AOBasisSize);
            // for every atom
            for (unsigned i = 0; i < _grid.size(); i++) {
	      // for each point in atom grid
                for (unsigned j = 0; j < _grid[i].size(); j++) {


                   boost::timer::cpu_times t0 = cpu_t.elapsed();


                    // get value of orbitals at each gridpoint (vector as 1D boost matrix object -> prod )
                    ub::matrix<double> AOgrid = ub::zero_matrix<double>(basis->_AOBasisSize, 1);

		    // get value of density gradient at each gridpoint
                    ub::matrix<double> gradAOgrid = ub::zero_matrix<double>(basis->_AOBasisSize, 3); // for Gradients of AOs

                    
		    // evaluate AO Functions for all shells
                    for (vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++) {

                        // for density, fill sub-part of AOatgrid
                        ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), 0, 1);
                        (*_row)->EvalAOspace(_AOgridsub, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);

                        // for density gradient, fill sub-part of gradAOgrid
                        ub::matrix_range< ub::matrix<double> > _gradAO = ub::subrange(gradAOgrid, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), 0, 3);
                        (*_row)->EvalAOGradspace(_gradAO, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);


                    }

                    boost::timer::cpu_times t1 = cpu_t.elapsed();
                    _t_AOgrid += (t1.wall-t0.wall)/1e9;
                    
		    // rho(r) = trans(AOatgrid) * DMAT * AOatgrid ?
		    // rho(r) = sum_{ab}{X_a * D_{ab} * X_b} =sum_a{ X_a * sum_b{D_{ab}*X_b}}
		    ub::matrix<double> _tempmat = ub::prod(_density_matrix,AOgrid); // tempmat can be reused for density gradient
		    double rho = ub::prod(ub::trans(AOgrid),_tempmat)(0,0);

		    if ( rho < 1.e-15 ) continue; // skip the rest, if density is very small

                    boost::timer::cpu_times t2 = cpu_t.elapsed();
                    _t_rho += (t2.wall-t1.wall)/1e9;
                    
                    
                    // density gradient as grad(n) = sum_{ab}[D_{ab} (X_b grad(X_a) + X_a grad(X_b)]
		    // grad(r) = sum_{ab}{X_b grad(X_a)*D_{ab} } + sum_{ab}{X_a grad(X_b)*D_{ab}}
		    //         = sum_{ab}{grad(X_a) * D_{ab} * X_b} + sum_{ab}{grad(X_b) * D_{ab} * X_a}
		    //         = 2.0 * sum_{ab}{grad(X_a) * D_{ab}*X_b}
		    ub::matrix<double> grad_rho = 2.0 * ub::prod(ub::trans(gradAOgrid),_tempmat);

                    //                    cout << rho << " " << grad_rho(0,0) << " " << grad_rho(1,0) << " " << grad_rho(0,1) << endl;
                    //exit(0);
                    
                    boost::timer::cpu_times t3 = cpu_t.elapsed();
                    _t_grad_rho += (t3.wall-t2.wall)/1e9;
                    
   
                    // get XC for this density_at_grid
                    double f_xc;      // E_xc[n] = int{n(r)*eps_xc[n(r)] d3r} = int{ f_xc(r) d3r }
                    double df_drho;   // v_xc_rho(r) = df/drho
                    double df_dsigma; // df/dsigma ( df/dgrad(rho) = df/dsigma * dsigma/dgrad(rho) = df/dsigma * 2*grad(rho))




                    // evaluate via LIBXC, if compiled, otherwise, go via own implementation
#ifdef LIBXC
		    double sigma = ub::prod(ub::trans(grad_rho),grad_rho)(0,0);
                    double exc[1];
                    double vsigma[1]; // libxc 
                    double vrho[1]; // libxc df/drho
                    switch (xfunc.info->family) {
                        case XC_FAMILY_LDA:
                            xc_lda_exc_vxc(&xfunc, 1, &rho, exc, vrho);
                            break;
                        case XC_FAMILY_GGA:
                        case XC_FAMILY_HYB_GGA:
                            xc_gga_exc_vxc(&xfunc, 1, &rho, &sigma, exc, vrho, vsigma);
                            break;
                    }
                    f_xc      = exc[0];
                    df_drho   = vrho[0];
                    df_dsigma = vsigma[0];

                    // via libxc correlation part only
                    switch (cfunc.info->family) {
                        case XC_FAMILY_LDA:
                            xc_lda_exc_vxc(&cfunc, 1, &rho, exc, vrho);
                            break;
                        case XC_FAMILY_GGA:
                        case XC_FAMILY_HYB_GGA:
                            xc_gga_exc_vxc(&cfunc, 1, &rho, &sigma, exc, vrho, vsigma);
                            break;
                    }

                    df_drho   += vrho[0];
                    df_dsigma += vsigma[0];

#else
                    _xc.getXC(-2, rho, grad_rho(0,0), grad_rho(1,0), grad_rho(2,0), f_xc, df_drho, df_dsigma);
#endif

                    boost::timer::cpu_times t4 = cpu_t.elapsed();
                    _t_vxc += (t4.wall-t3.wall)/1e9;
                    
		    // density part
		    ub::matrix<double> _addXC = _grid[i][j].grid_weight * df_drho * AOgrid * 0.5;
                    // Exchange correlation energy
                    EXC += _grid[i][j].grid_weight * rho * f_xc;
                    boost::timer::cpu_times t5 = cpu_t.elapsed();
                    _t_AOxc_rho += (t5.wall-t4.wall)/1e9;
                    
                    
		    // gradient part
		    int size = gradAOgrid.size1();
                    ub::matrix_range< ub::matrix<double> > _gradAO_x = ub::subrange(gradAOgrid, 0, size , 0, 1);
                    ub::matrix_range< ub::matrix<double> > _gradAO_y = ub::subrange(gradAOgrid, 0, size , 1, 2);
                    ub::matrix_range< ub::matrix<double> > _gradAO_z = ub::subrange(gradAOgrid, 0, size , 2, 3);
		    //_addXC += 4.0*df_dsigma * _grid[i][j].grid_weight *(grad_rho(0,0) * _gradAO_x + grad_rho(1,0) * _gradAO_y + grad_rho(2,0) * _gradAO_z );
                    _addXC += 2.0 *df_dsigma * _grid[i][j].grid_weight *(grad_rho(0,0) * _gradAO_x + grad_rho(1,0) * _gradAO_y + grad_rho(2,0) * _gradAO_z );
                    
                    //ub::matrix<double> _XCmat1 = ub::prod(_addXC1,ub::trans(AOgrid));
                    
                    boost::timer::cpu_times t6 = cpu_t.elapsed();
                    _t_AOxc_grad += (t6.wall-t5.wall)/1e9;

		    // finally combine
		    XCMAT +=  ub::prod(_addXC,ub::trans(AOgrid));

                    boost::timer::cpu_times t7 = cpu_t.elapsed();
                    _t_sum += (t7.wall-t6.wall)/1e9;
                    
                    
                } // j: for each point in atom grid
            } // i: for each atom grid

            
            XCMAT += ub::trans(XCMAT);
            
            
            
            const ub::vector<double> DMAT_array = _density_matrix.data();
            const ub::vector<double> XCMAT_array = XCMAT.data();
            for ( unsigned i = 0; i < DMAT_array.size(); i++ ){
            
                EXC -= DMAT_array[i] * XCMAT_array[i];
                
            }
            
            cout << " Time AO functions: " << _t_AOgrid << endl;
            cout << " Time rho         : " << _t_rho << endl;
            cout << " Time grad rho    : " << _t_grad_rho << endl;
            cout << " Time Vxc         : " << _t_vxc << endl;
            cout << " Time AOxc rho    : " << _t_AOxc_rho << endl;
            cout << " Time AOxc grad   : " << _t_AOxc_grad << endl;
            cout << " Time AOxc sum    : " << _t_sum << endl;
                                 boost::timer::cpu_times texit = cpu_t.elapsed();
                                 _t_total += (texit.wall-tenter.wall)/1e9;
                                 
                                 cout << " Time TOTAL   : " << _t_total << endl;
            return XCMAT;

        }

        
 /*       
                
        // numerically integrate the elements of the AOXC matrix
        ub::matrix<double> NumericalIntegration::StupidIntegrateVXC(ub::matrix<double>& _density_matrix, AOBasis* basis){
            
            // TODO: switch XC functionals implementation from LIBXC to base own calculation
            #ifdef LIBXC
            xc_func_type xfunc; // handle for exchange functional
            xc_func_type cfunc; // handle for correleation functions

            // define PBE here (should be some optional setting))
            int xfunc_id = XC_GGA_X_PBE;
            int cfunc_id = XC_GGA_C_PBE;
            
            if(xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED) != 0){
               fprintf(stderr, "Functional '%d' not found\n", xfunc_id);
               exit(1);
            }
            
            if(xc_func_init(&cfunc, cfunc_id, XC_UNPOLARIZED) != 0){
               fprintf(stderr, "Functional '%d' not found\n", cfunc_id);
               exit(1);
            }
#else
             ExchangeCorrelation _xc;
#endif
            
            //printf("The exchange functional '%s' is defined in the reference(s):\n%s\n", xfunc.info->name, xfunc.info->refs);
            //printf("The correlation functional '%s' is defined in the reference(s):\n%s\n", cfunc.info->name, cfunc.info->refs);

            // xc_func_end(&xfunc);

           
            ub::matrix<double> XCMAT = ub::zero_matrix<double>(basis->_AOBasisSize, basis->_AOBasisSize);
            const ub::vector<double> DMAT_array = _density_matrix.data();
            // for every gridpoint
            for (int i = 0; i < _grid.size(); i++) {
                for (int j = 0; j < _grid[i].size(); j++) {
                    // get value of orbitals at each gridpoint
                    ub::matrix<double> AOatgrid = ub::zero_matrix<double>(basis->_AOBasisSize, 1);

                    ub::matrix<double> AODerXatgrid = ub::zero_matrix<double>(basis->_AOBasisSize, 1); // for Gradients of AOs
                    ub::matrix<double> AODerYatgrid = ub::zero_matrix<double>(basis->_AOBasisSize, 1); // for Gradients of AOs
                    ub::matrix<double> AODerZatgrid = ub::zero_matrix<double>(basis->_AOBasisSize, 1); // for Gradients of AOs

                    for (vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++) {

                        // for density
                        ub::matrix_range< ub::matrix<double> > _AOatgridsub = ub::subrange(AOatgrid, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), 0, 0);
                        (*_row)->EvalAOspace(_AOatgridsub, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);

                        // for density gradient  
                        ub::matrix_range< ub::matrix<double> > _AODerXatgridsub = ub::subrange(AODerXatgrid, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), 0, 0);
                        ub::matrix_range< ub::matrix<double> > _AODerYatgridsub = ub::subrange(AODerYatgrid, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), 0, 0);
                        ub::matrix_range< ub::matrix<double> > _AODerZatgridsub = ub::subrange(AODerZatgrid, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), 0, 0);
                        (*_row)->EvalAOGradspace(_AODerXatgridsub, _AODerYatgridsub, _AODerZatgridsub, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);


                    }


                    ub::matrix<double> _AOmatrix_at_grid = ub::prod(AOatgrid, ub::trans(AOatgrid));

                    // density at grid point is sum of element-wise product of density matrix x _AOmatrix
                    ub::vector<double> _AO_array = _AOmatrix_at_grid.data();
                    double density_at_grid = 0.0;
                    for (int _i = 0; _i < DMAT_array.size(); _i++) {
                        density_at_grid += DMAT_array(_i) * _AO_array(_i);
                    }

                    // density gradient as grad(n) = sum_{ab}[D_{ab} (X_b grad(X_a) + X_a grad(X_b)]
                    // x,y-z-components of gradient
                    ub::matrix<double> _AODerXmatrix_at_grid = ub::prod(AODerXatgrid, ub::trans(AOatgrid)) + ub::prod(AOatgrid, ub::trans(AODerXatgrid));
                    ub::matrix<double> _AODerYmatrix_at_grid = ub::prod(AODerYatgrid, ub::trans(AOatgrid)) + ub::prod(AOatgrid, ub::trans(AODerYatgrid));
                    ub::matrix<double> _AODerZmatrix_at_grid = ub::prod(AODerZatgrid, ub::trans(AOatgrid)) + ub::prod(AOatgrid, ub::trans(AODerZatgrid));

                    ub::vector<double> _AODerX_array = _AODerXmatrix_at_grid.data();
                    ub::vector<double> _AODerY_array = _AODerYmatrix_at_grid.data();
                    ub::vector<double> _AODerZ_array = _AODerZmatrix_at_grid.data();

                    double densityDerX_at_grid = 0.0;
                    double densityDerY_at_grid = 0.0;
                    double densityDerZ_at_grid = 0.0;
                    for (int _i = 0; _i < DMAT_array.size(); _i++) {
                        densityDerX_at_grid += DMAT_array(_i) * _AODerX_array(_i);
                        densityDerY_at_grid += DMAT_array(_i) * _AODerY_array(_i);
                        densityDerZ_at_grid += DMAT_array(_i) * _AODerZ_array(_i);
                    }




                    // get XC for this density_at_grid
                    double f_xc; // E_xc[n] = int{n(r)*eps_xc[n(r)] d3r} = int{ f_xc(r) d3r }
                    double df_drho; // v_xc_rho(r) = df/drho
                    double df_dsigma; //df/dsigma ( df/dgrad(rho) = df/dsigma * dsigma/dgrad(rho) = df/dsigma * 2*grad(rho))




                    // evaluate via LIBXC, if compiled, otherwise, go via own implementation
#ifdef LIBXC
                    double sigma_at_grid = densityDerX_at_grid * densityDerX_at_grid + densityDerY_at_grid * densityDerY_at_grid + densityDerZ_at_grid*densityDerZ_at_grid;

                    double vsigma[1]; // libxc 
                    double vrho[1]; // libxc df/drho
                    switch (xfunc.info->family) {
                        case XC_FAMILY_LDA:
                            xc_lda_vxc(&xfunc, 1, &density_at_grid, vrho);
                            break;
                        case XC_FAMILY_GGA:
                        case XC_FAMILY_HYB_GGA:
                            xc_gga_vxc(&xfunc, 1, &density_at_grid, &sigma_at_grid, vrho, vsigma);
                            break;
                    }
                    df_drho = vrho[0];
                    df_dsigma = vsigma[0];

                    // via libxc correlation part only
                    switch (cfunc.info->family) {
                        case XC_FAMILY_LDA:
                            xc_lda_vxc(&cfunc, 1, &density_at_grid, vrho);
                            break;
                        case XC_FAMILY_GGA:
                        case XC_FAMILY_HYB_GGA:
                            xc_gga_vxc(&cfunc, 1, &density_at_grid, &sigma_at_grid, vrho, vsigma);
                            break;
                    }

                    df_drho += vrho[0];
                    df_dsigma += vsigma[0];

#else
                    _xc.getXC("PBE", density_at_grid, densityDerX_at_grid, densityDerY_at_grid, densityDerZ_at_grid, f_xc, df_drho, df_dsigma);
#endif

                    // exit(0);
                    // cout << " out rho : " << density_at_grid << " vxc " << v << endl;
                    XCMAT += _grid[i][j].grid_weight * df_drho * _AOmatrix_at_grid;
                    // gradient corrections
                    XCMAT += _grid[i][j].grid_weight * df_dsigma * 2.0 * (densityDerX_at_grid * _AODerXmatrix_at_grid + densityDerY_at_grid * _AODerYmatrix_at_grid + densityDerZ_at_grid * _AODerZmatrix_at_grid);



                } // j: for each point in atom grid
            } // i: for each atom grid

            return XCMAT;

        } */

        // numerically integrate the electron density
        double NumericalIntegration::IntegrateDensity(ub::matrix<double>& _density_matrix, AOBasis* basis){
            
            double result = 0.0;
            const ub::vector<double> DMAT_array=_density_matrix.data();
             // for every gridpoint
            #pragma omp parallel for 
            for (unsigned i = 0; i < _grid.size(); i++) {
                for (unsigned j = 0; j < _grid[i].size(); j++) {
                    // get value of orbitals at each gridpoint
                    ub::matrix<double> tmat = ub::zero_matrix<double>(basis->_AOBasisSize, 1);

                    for (vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++) {

                        ub::matrix_range< ub::matrix<double> > _submatrix = ub::subrange(tmat,0,1, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc());
                        (*_row)->EvalAOspace(_submatrix, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);
                    }
                    
                    
                    ub::matrix<double> _AOmatrix_at_grid = ub::prod( ub::trans(tmat),tmat);
                    
                    // density at grid point is sum of element-wise product of density matrix x _AOmatrix
                    ub::vector<double> _AO_array  =_AOmatrix_at_grid.data();
                    double density_at_grid = 0.0;
                    for ( unsigned _i =0; _i < DMAT_array.size(); _i++ ){
                        density_at_grid += DMAT_array(_i)*_AO_array(_i);
                    }   
                    

                    
                    _grid[i][j].grid_density=density_at_grid;
                    result += _grid[i][j].grid_weight * density_at_grid;
                }
            } // gridpoints end
            
            density_set=true;
            return result;
        } 
            
        double NumericalIntegration::SetGridToCharges(std::vector< QMAtom* > & _local_atomlist){
            
            _grid.empty();
            int p=1;
            std::vector< GridContainers::integration_grid > _Mad;
            for(std::vector< QMAtom* >::iterator it=_local_atomlist.begin(); it!=_local_atomlist.end(); ++it){
                GridContainers::integration_grid el;
                QMAtom* atom=*it;
                el.grid_density=std::pow(-1.0,p);
                el.grid_weight=1.0;
                el.grid_x=atom->x*tools::conv::ang2bohr;
                el.grid_y=atom->y*tools::conv::ang2bohr;
                el.grid_z=atom->z*tools::conv::ang2bohr;
                //cout<< atom->charge << "@ "<<el.grid_x<<'\t'<<el.grid_y<<'\t'<<el.grid_z<<endl;
                _Mad.push_back(el);
                p++;
            }
            _grid.push_back(_Mad);
            
            return 0.0;
        }
        
        double NumericalIntegration::IntegratePotential(ub::vector<double> rvector){
            
            double result = 0.0;

            if(density_set){
                 for (unsigned i = 0; i < _grid.size(); i++) {
                 for (unsigned j = 0; j < _grid[i].size(); j++) {
                     double dist=sqrt((_grid[i][j].grid_x-rvector(0))*(_grid[i][j].grid_x-rvector(0))+(_grid[i][j].grid_y-rvector(1))*(_grid[i][j].grid_y-rvector(1))+(_grid[i][j].grid_z-rvector(2))*(_grid[i][j].grid_z-rvector(2)));
                     result -= _grid[i][j].grid_weight * _grid[i][j].grid_density/dist;
                     }
                 }
             } 
            else{
                throw std::runtime_error("Density not calculated");
            }
            
            return result;   
        }
        
        
        
        double NumericalIntegration::IntegrateEnergy_w_PBC(ub::vector<double> rvector, double boxLen[3]){
            
            E_rspace = 0.0;
            E_kspace = 0.0;
            double E_r=0, E_k=0;
            double cutoff=min(min(boxLen[0], boxLen[1]), boxLen[2])/2.0;
            double vol=boxLen[0]*boxLen[1]*boxLen[2];
            
            if(density_set){
                //for (unsigned i = 0; i < _grid.size(); i++) {
                    //for (unsigned j = 0; j < _grid[i].size(); j++) {
                #pragma omp parallel for reduction(+:E_r)
                for (unsigned i = 0; i < _Madelung_grid[0].size(); i++) {
                    for (unsigned j = 0; j < _Madelung_grid[0].size(); j++) {

                
                        //charge at this point
                        //double q = -_grid[i][j].grid_weight * _grid[i][j].grid_density; //density is neg of charge
                        double qi = -_Madelung_grid[0][i].grid_weight * _Madelung_grid[0][i].grid_density; //density is neg of charge
                        double qj = -_Madelung_grid[0][j].grid_weight * _Madelung_grid[0][j].grid_density; //density is neg of charge
                        
                        //r-space sum
                        double dif[3];
//                        dif[0] = _grid[i][j].grid_x-rvector(0);
//                        dif[1] = _grid[i][j].grid_y-rvector(1);
//                        dif[2] = _grid[i][j].grid_z-rvector(2);
                        dif[0] = _Madelung_grid[0][j].grid_x-_Madelung_grid[0][i].grid_x;
                        dif[1] = _Madelung_grid[0][j].grid_y-_Madelung_grid[0][i].grid_y;
                        dif[2] = _Madelung_grid[0][j].grid_z-_Madelung_grid[0][i].grid_z;

                        for(int k=0; k<3; k++){
                            if(std::abs(dif[k])>boxLen[k]*0.5) //correct for for bond crossing PBC, if it exists
                                if(dif[k]>0)    //i.x>j.x
                                    dif[k]-=boxLen[k];
                                else            //i.x<j.x
                                    dif[k]+=boxLen[k];
                        }
                        double dist=sqrt((dif[0]*dif[0])+(dif[1]*dif[1])+(dif[2]*dif[2])); //in bohr
                        double ER=0.5*qi*qj*((erfc(alpha*dist)/dist) - tools::conv::Pi/(vol*vol*alpha*alpha));

                        if(dist<1.0e-12){ //point is in the same spot as we are evaluating potential
                            E_r -= (alpha/sqrt(tools::conv::Pi))*qi*qj; //self correction
                        }
                        else if(dist<=cutoff)
                        {
                            E_r += ER;
                        }
						
                    }//j
                }//i 
                E_rspace = E_r;
                cout<< "E: r-space sum: " <<E_rspace<<endl;

                //k-space sum
                #pragma omp parallel for reduction(+:E_k)
                for (unsigned i = 0; i < _Madelung_grid[0].size(); i++) {
                    vec r(_Madelung_grid[0][i].grid_x,_Madelung_grid[0][i].grid_y,_Madelung_grid[0][i].grid_z);
                    double qi = -_Madelung_grid[0][i].grid_weight * _Madelung_grid[0][i].grid_density; //density is neg of charge
                    int nKpoints=numK[0]*numK[1]*numK[2];
                    double* kp;
                    for(int index=1; index<nKpoints; index++){
                            kp = &(Kcoord[index*3]);
                            vec K(kp[0], kp[1], kp[2]);
                            std::complex<double> Kr(0, K*r); // ik dot r

                            
                            double EK=0.5*qi*prefactor[index]*(Rho_k[index]*std::exp(Kr)).real();
                            E_k+=EK;

                    }//index
                }//i
                
                E_kspace = E_k;
                cout<< "E: k-space sum: " <<E_kspace<<endl;
                
                E_erfc=erfc(alpha*cutoff);

            } 
            else{
               throw std::runtime_error("Density not calculated");
            }
            
            return E_kspace+E_rspace;   
        }
        
        
        
        double NumericalIntegration::IntegratePotential_w_PBC(ub::vector<double> rvector, double boxLen[3]){
            
            double result = 0.0;
            double cutoff=min(min(boxLen[0], boxLen[1]), boxLen[2])/2.0;
            double vol=boxLen[0]*boxLen[1]*boxLen[2];
            
            if(density_set){
//                #pragma omp parallel for reduction(+:result)
                for (unsigned i = 0; i < _grid.size(); i++) {
                    for (unsigned j = 0; j < _grid[i].size(); j++) {
//                for (unsigned i = 0; i < _Madelung_grid.size(); i++) {
//                    for (unsigned j = 0; j < _Madelung_grid[i].size(); j++) {

                        //charge at this point
                        double q = -_grid[i][j].grid_weight * _grid[i][j].grid_density; //density is neg of charge
                        //double q = -_Madelung_grid[i][j].grid_weight * _Madelung_grid[i][j].grid_density; //density is neg of charge
                        
                        //r-space sum
                        double dif[3];
                        dif[0] = _grid[i][j].grid_x-rvector(0);
                        dif[1] = _grid[i][j].grid_y-rvector(1);
                        dif[2] = _grid[i][j].grid_z-rvector(2);
//                        dif[0] = _Madelung_grid[i][j].grid_x-rvector(0);
//                        dif[1] = _Madelung_grid[i][j].grid_y-rvector(1);
//                        dif[2] = _Madelung_grid[i][j].grid_z-rvector(2);

                        for(int k=0; k<3; k++){
                            if(std::abs(dif[k])>boxLen[k]*0.5) //correct for for bond crossing PBC, if it exists
                                if(dif[k]>0)    //i.x>j.x
                                    dif[k]-=boxLen[k];
                                else            //i.x<j.x
                                    dif[k]+=boxLen[k];
                        }
                        double dist=sqrt((dif[0]*dif[0])+(dif[1]*dif[1])+(dif[2]*dif[2])); //in bohr
                        double potR=q*((erfc(alpha*dist)/dist) - tools::conv::Pi/(vol*vol*alpha*alpha)); //erfc and shift everage pot to 0
                        
                        if(dist<1.0e-12){ //point is in the same spot as we are evaluating potential
                            result-=2.0*(alpha/sqrt(tools::conv::Pi))*q; //self correction
                        }
                        else if(dist<=cutoff)
                        {
                            //double potR=q/dist;
                            result += potR;
                        }
						
                    }//j
                }//i 
				
                //cout<< "r-space sum: " <<result<<endl;

                //k-space sum
                vec r(rvector(0),rvector(1),rvector(2));
                int nKpoints=numK[0]*numK[1]*numK[2];
                double ksum=0.0;
                //exclude k={0,0,0} (index=0)
//                #pragma omp parallel for reduction(+:ksum)
                for(int index=1; index<nKpoints; index++){
                        double* kp = &(Kcoord[index*3]);
                        vec K(kp[0], kp[1], kp[2]);
                        std::complex<double> Kr(0, K*r); // ik dot r

                        double potK=prefactor[index]*(Rho_k[index]*std::exp(Kr)).real();
                        ksum+=potK;

                }//index
                result+=ksum;


                //cout<< "k-space sum: " <<ksum<<endl;
                //cout<< "final Ewald sum: "<<result<<endl;

                //dipole moment correction
                //This is only usefull if we have a countably infinite
                //sphere of periodic cells embeded in a dielectric.
                //If material is trully periodic ("metallic boundary 
                //conditions"), then this correction is 0
				
                        

            } 
            else{
               throw std::runtime_error("Density not calculated");
            }
            
            return result;   
        }
        
        void NumericalIntegration::FreeKspace(){
            delete[] Rho_k;
            delete[] Kcoord;
            delete[] prefactor;
        }
        
		
        /*
         *	Calculate and return the Ewald coefficient (alpha) from cutoff distance
         *	and requested tolerance.
         */
        void NumericalIntegration::findAlpha(double Rc, double dtol){
                double x=1.0;
                int i=20;
                while (erfc(x*Rc)/Rc > dtol){
                        x*=2.0;
                        i++;
                }
                double lastx=x*2.0;
                double newx;
                while (true){
                    double reldif=(erfc(x*Rc)/Rc - dtol)/dtol;
                    if(std::abs(reldif)<1.0e-3 || i==0)
                        break;
                    
                    if(reldif>0)
                        newx=(x+lastx)*0.5;
                    else
                        newx=x*0.5;
                    lastx=x;
                    x=newx;
                    i--;
                }
                alpha=x;
        }
        
        /**
         * 
         * @param _atoms
         * @param boxLen
         * @param Kspacing in Angstroms
         */
        void NumericalIntegration::PrepKspaceDensity(double boxLen[3], double Kspacing, int natomsonside, double ext_alpha, std::vector< QMAtom* > & _local_atomlist, bool ECP){
            
            cout<<"box is "<< boxLen[0] << " "<< boxLen[1] << " "<< boxLen[2] << endl;
            
            //add nuclear charges to _grid to avoid having to do periodic calculations on them separately;
            //_Madelung_grid.clear();
            std::vector< GridContainers::integration_grid > _Nuc;
            Elements _elements;
            for(std::vector< QMAtom* >::iterator it=_local_atomlist.begin(); it!=_local_atomlist.end(); ++it){
                GridContainers::integration_grid el;
                QMAtom* atom=*it;
                if (ECP) {
                    el.grid_density = -_elements.getNucCrgECP(atom->type);
                } else {
                    el.grid_density = -_elements.getNucCrg(atom->type);
                }
                el.grid_weight=1.0;
                el.grid_x=atom->x*tools::conv::ang2bohr;
                el.grid_y=atom->y*tools::conv::ang2bohr;
                el.grid_z=atom->z*tools::conv::ang2bohr;
                //cout<< atom->charge << "@ "<<el.grid_x<<'\t'<<el.grid_y<<'\t'<<el.grid_z<<endl;
                _Nuc.push_back(el);
            }
            _grid.push_back(_Nuc);
            

//
//            //fill Madelung grid with density grid, all in to _Madelung_grid[0], so that energy calculation foesn't have a quadrupple nested loop
//            _Madelung_grid=_grid;

            
            
            //cout<<"points in _Mad: " <<_Mad.size() << endl;
            //cout<<"_Mads in _Madelung_grid: " <<_Madelung_grid.size() << endl;
            
            
			
			//compute alpha
            double cutoff=min(min(boxLen[0], boxLen[1]), boxLen[2])/2.0;
//            findAlpha(cutoff, 1.0e-10);
            alpha=ext_alpha;
            double fourasq=4.0*alpha*alpha;
            cout<<"found alpha = "<< alpha <<"\t rel err of r-sum ~ " << erfc(alpha*cutoff)/cutoff<< endl;
//            cout<<"found alpha = "<< alpha <<"\t";
			
			
            //this is going to be slow, as points we have density for are not on a periodic grid,
            //so will use simple Ewald summation, not any of the FFT methods.
			
			//find number of k-vectors in each lattice direction
            double minSq=1.0e15;
            for(int i=0; i<3; i++){
                //numK[i]=1+(boxLen[i]/Kspacing);
                //if(numK[i]%2 == 0) numK[i]++; //keep it odd
				
                //find maxK
                int maxK;
                double maxKsq;
                double err;
                double twoPiL=2.0*tools::conv::Pi/boxLen[i];
                for(maxK=2; true; maxK++){
                    maxKsq=maxK*maxK*twoPiL*twoPiL;
                    err=std::exp(-maxKsq/fourasq)/maxKsq;
                    if(err<1.0e-7){
                            break;
                    }
                }
                
                maxK=16;

                numK[i]=2*maxK+1;
                minSq=min(maxKsq, minSq);
            }
            cout<<"numK={"<<numK[0]<<", "<<numK[1]<<", "<<numK[2]<<"}"<<endl;
            cout<<"rel err of k-sum ~ " << std::exp(-minSq/fourasq)/minSq << endl;
				
            
            //allocate
            int nKpoints=numK[0]*numK[1]*numK[2];
            Rho_k=new std::complex<double>[nKpoints];
            Kcoord=new double[nKpoints*3];
            prefactor=new double[nKpoints];
            
            //fill Kcoord and prefactor
            // [0, 1, 2, 3, ..., N/2 -1, N/2, -N/2, -N/2 +1, -N/2 +1, ..., -3, -2, -1] (note no 0)
            // k={0,0,0} is the very first entry (index=0)
            
            
            double invvolume = 1.0/(boxLen[0]*boxLen[1]*boxLen[2]);
            double pre = invvolume*4.0*tools::conv::Pi;
            
            //cout<<"pre=4*Pi/Volume="<<pre<<endl;
            //cout<<"fourasq=4*alpha^2="<<fourasq<<endl;

            #pragma omp parallel for
            for(int l=0; l<numK[0]; l++){
                int L, M, N;
                double ksq;
                if (l>numK[0]/2){ L=l-numK[0];}
                else { L=l;}
                //cout<<L<<endl;
                for(int m=0; m<numK[1]; m++){
                    if (m>numK[1]/2){ M=m-numK[1];}
                    else { M=m;}
                    for(int n=0; n<numK[2]; n++){
                        if (n>numK[2]/2){ N=n-numK[2];}
                        else { N=n;}
						
			int index = l*(numK[1]*numK[2])+m*numK[2]+n;
					  
                        double* kp = &(Kcoord[index*3]);
                        kp[0]=2.0*tools::conv::Pi * L /boxLen[0];
                        kp[1]=2.0*tools::conv::Pi * M /boxLen[1];
                        kp[2]=2.0*tools::conv::Pi * N /boxLen[2];
                        
                        ksq=(kp[0]*kp[0])+(kp[1]*kp[1])+(kp[2]*kp[2]);
                        //cout<<"{L,M,N}={"<<L<<", "<<M<<", "<<N<<"}\n";
//                        cout<<"{l,m,n}={"<<l<<", "<<m<<", "<<n<<"}\t";
//                        cout<<"k={"<<kp[0]<<", "<<kp[1]<<", "<<kp[2]<<"}\tksq="<<ksq<<"\t";
//                        cout<<"prefactor="<<(pre/ksq) * std::exp(-ksq/fourasq)<<endl;
                        prefactor[index]=(pre/ksq) * std::exp(-ksq/fourasq);
                        
                        //index++;
                    }
                }
            }
            
            
            //fill Rho_k
            #pragma omp parallel for
            for(int index=0; index<nKpoints; index++){
                double* kp = &(Kcoord[index*3]);
                vec K(kp[0], kp[1], kp[2]);
                Rho_k[index]=0.0;
                for (unsigned i = 0; i < _grid.size(); i++) {
                    for (unsigned j = 0; j < _grid[i].size(); j++) {
                        vec r(_grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);
                        std::complex<double> nKr(0, -K*r); // -ik dot r
                        Rho_k[index] -= _grid[i][j].grid_weight * _grid[i][j].grid_density * std::exp(nKr); //density is neg of charge
                    }//j
                }//i
                
//                for (unsigned i = 0; i < _Madelung_grid.size(); i++) {
//                    for (unsigned j = 0; j < _Madelung_grid[i].size(); j++) {
//                        vec r(_Madelung_grid[i][j].grid_x, _Madelung_grid[i][j].grid_y, _Madelung_grid[i][j].grid_z);
//                        std::complex<double> nKr(0, -K*r); // -ik dot r
//                        Rho_k[index] -= _Madelung_grid[i][j].grid_weight * _Madelung_grid[i][j].grid_density * std::exp(nKr); //density is neg of charge
////                        cout<<"k={"<<kp[0]<<", "<<kp[1]<<", "<<kp[2]<<"}\t";
////                        cout<<"r={"<<_Madelung_grid[i][j].grid_x<<", "<<_Madelung_grid[i][j].grid_y<<", "<<_Madelung_grid[i][j].grid_z<<"}\t";
////                        cout<<"-ikr={"<<nKr<<"}\t exp(-iKr)="<<std::exp(nKr)<<"\t Rho_k="<<Rho_k[index]<<endl;
//                    }
//                }
                
            }//index
            
        }
        
        
        double NumericalIntegration::IntegrateDensity_Atomblock(ub::matrix<double>& _density_matrix, AOBasis* basis){   
         
            double result=0.0;
            // timers for testing
            //boost::timer::cpu_timer cpu_t;
            //cpu_t.start();
            //double _t_AOvals = 0.0;
            //double _t_rho = 0.0;
            //double _t_grad_rho = 0.0;
            //double _t_vxc =0.0;
            //double _t_AOxc_rho=0.0;
            //double _t_AOxc_grad=0.0;
            //double _t_EXC1=0.0;
            //double _t_EXC2=0.0;
            //double _t_sum = 0.0;
            //double _t_total = 0.0;
            //boost::timer::cpu_times tenter = cpu_t.elapsed();
             
            // generate a list of shells for each atom
            typedef vector< AOShell* >::iterator AOShellIterator;
            vector< vector< AOShellIterator > > _atomshells;
            vector< AOShellIterator > _singleatom;

            vector < int > _startIdx;
            vector < int > _blocksize;

            int _atomindex = 0;
            int _Idx       = 0;
            int _size      = 0;
            
            for (vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++) {
                
                
                if ( (*_row)->getIndex() == _atomindex ){
                    
                    _singleatom.push_back(_row);
                    _size += (*_row)->getNumFunc();
                    
                    
                } else {
                    
                    // append _singleatom to _atomshells
                    _atomshells.push_back(_singleatom);
                    _startIdx.push_back( _Idx );
                    _blocksize.push_back(_size);
                    // reset _singleatom
                    _singleatom.clear();
                    _size = (*_row)->getNumFunc();
                    _Idx       = (*_row)->getStartIndex();
                    _singleatom.push_back(_row);
                    _atomindex = (*_row)->getIndex();
                    
                }
                   
            }
            
            _atomshells.push_back(_singleatom);
            _startIdx.push_back( _Idx );
                    _blocksize.push_back(_size);
            //cout << " Number of atoms " << _atomshells.size() << endl;
            
            //for ( unsigned iatom = 0 ; iatom < _atomshells.size(); iatom++ ){
            //    cout << "atom " << iatom << " number of shells " << _atomshells[iatom].size() << " block start " << _startIdx[iatom] << " functions in atom " << _blocksize[iatom] << endl; 
            //}

           
            // setup a list of min decay constants per atom
            // for every shell
            _atomindex = 0;
            double _decaymin = 1e7;
            vector< double > _minimal_decay;
            vector < vec > _positions;
            vec _localpos = (*basis->firstShell())->getPos();
            for ( vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++   ) {
                               
                 if ( (*_row)->getIndex() == _atomindex ){
                     
                     // check all decay constants in this shell
                     for (AOShell::GaussianIterator itg = (*_row)->firstGaussian(); itg != (*_row)->lastGaussian(); itg++) {
                         AOGaussianPrimitive* gaussian = *itg;
                         double _decay = gaussian->decay;
                         if (_decay < _decaymin) {
                             _decaymin = _decay;
                         } // decay min check
                     
                     } // Gaussian Primitives 
                     
                 } else {  // if shell belongs to the actual atom
                     // add to mininal_decay vector
                     _minimal_decay.push_back(_decaymin);
                     _positions.push_back( _localpos );
                     // reset counters
                     _decaymin = 1e7;
                     _localpos = (*_row)->getPos();

                     _atomindex++;
                     
                     // check all decay constants in this shell
                     for (AOShell::GaussianIterator itg = (*_row)->firstGaussian(); itg != (*_row)->lastGaussian(); itg++) {
                         AOGaussianPrimitive* gaussian = *itg;
                         double _decay = gaussian->decay;
                         if (_decay < _decaymin) {
                             _decaymin = _decay;
                         } // decay min check
                     
                     } // Gaussian Primitives                                       
                 }
            } // all shells
                 
            // push final atom
            _minimal_decay.push_back(_decaymin);
            _positions.push_back( _localpos );
             
            /* for ( int i =0; i < _minimal_decay.size(); i++){
                 
                 cout << "Atom " << i << " min decay " << _minimal_decay[i] <<  " at " << _positions[i] << endl; 
                 
             } */
                          
             // for each gridpoint, check the value of exp(-a*(r-R)^2) < 1e-10
             //                             = alpha*(r-R)^2 >~ 20.7
            
            vector< vector< vector<int> > > _significant_atoms;
            
            // each atomic grid
            for (unsigned i = 0; i < _grid.size(); i++) {
            
                vector< vector<int> > _significant_atoms_atomgrid;
                
                // each point of the atomic grid
                for (unsigned j = 0; j < _grid[i].size(); j++) {

                    vector<int> _significant_atoms_gridpoint;
                    vec grid;
                    grid.setX( _grid[i][j].grid_x);
                    grid.setY( _grid[i][j].grid_y);
                    grid.setZ( _grid[i][j].grid_z);
                    
                    // check all atoms
                    for ( unsigned iatom = 0 ; iatom < _minimal_decay.size(); iatom++){

                        vec dist = grid - _positions[iatom];
                        double distsq = dist.getX()*dist.getX() + dist.getY()*dist.getY()  + dist.getZ()*dist.getZ() ;
                        
                        // if contribution is smaller than -ln(1e-10), add atom to list
                        if ( (_minimal_decay[iatom] * distsq) < 20.7 ){
                            _significant_atoms_gridpoint.push_back(iatom);
                        }
                        
                    } // check all atoms

                    _significant_atoms_atomgrid.push_back(  _significant_atoms_gridpoint );
                   
                } // all points of this atom grid
                
                _significant_atoms.push_back(_significant_atoms_atomgrid);
               
            } // atomic grids
              
            /*
            for ( int i = 0; i < _significant_atoms.size(); i++ ){
                
                cout << " Atomgrid: " << i << endl;
                 for ( int j = 0; j < _significant_atoms[i].size(); j++ ){
                     
                     cout << " Atom: " << i << " gridpoint " << j << " significant atoms " << _significant_atoms[i][j].size() << endl;
                     
                 }
                
            }*/
             //exit(0);           
            
            int total_grid =0;
            int significant_grid = 0;
            for ( unsigned i = 0; i < _significant_atoms.size(); i++ ){
                
                total_grid += _grid[i].size(); 
                
                for ( unsigned j = 0; j < _significant_atoms[i].size(); j++ ){
                    
                    int gridpointsize = _significant_atoms[i][j].size();
                    significant_grid += gridpointsize*(gridpointsize+1);
                    
                } 
            }
            int natoms = _grid.size();
            
            total_grid = total_grid * ( natoms*(natoms+1) ) / 2;
            
            // cout << "Total number of atom blocks " << total_grid << " after removal of insignificant blocks " << significant_grid/2 << endl;
          //  cout << "# Atom blocks by removal of insignificant blocks reduced to " << 50*double(significant_grid)/double(total_grid) << "%" << endl;
  
            // parallelization: distribute over threads inside one atom
            int nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif

            std::vector<double> Density_thread;
            for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){ 
                Density_thread.push_back(0.0);
            }           
            
            // for every atom
            for (unsigned i = 0; i < _grid.size(); i++) {
	      // for each point in atom grid
                
                // number of points in this atomgrid
                int atom_points = _grid[i].size();
                // divide among threads
                int atom_points_per_thread = atom_points/nthreads;
                std::vector<int> _thread_start;
                std::vector<int> _thread_stop;
                for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                    _thread_start.push_back( i_thread * atom_points_per_thread );
                    _thread_stop.push_back( (i_thread + 1) * atom_points_per_thread );
                }
                // final stop must be size
                _thread_stop[nthreads-1] = atom_points;

              //  for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                    
               //     cout << "Thread " << i_thread << " start " << _thread_start[i_thread] << " stop " << _thread_stop[i_thread] << endl;                    
               // }
                
                #pragma omp parallel for
                for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                for (int j = _thread_start[i_thread]; j < _thread_stop[i_thread]; j++) {
                   //boost::timer::cpu_times t0 = cpu_t.elapsed();

                    // get value of orbitals at each gridpoint (vector as 1D boost matrix object -> prod )
                    //ub::matrix<double> AOgrid = ub::zero_matrix<double>(basis->_AOBasisSize, 1);

                   ub::matrix<double> AOgrid = ub::zero_matrix<double>(1, basis->_AOBasisSize); // TRY MORE USEFUL DATA
                   //ub::matrix<double> gradAOgrid = ub::zero_matrix<double>(3, basis->_AOBasisSize); 
                   ub::matrix<double>  rho_mat = ub::zero_matrix<double>(1,1);
                   //ub::matrix<double> grad_rho = ub::zero_matrix<double>(3,1);
                    
		    // evaluate AO Functions for all shells, NOW BLOCKWISE

                    // for each significant atom for this grid point
                    for ( unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size() ; sigrow++){
                    
                        // this atom
                        int rowatom = _significant_atoms[i][j][sigrow];
                                           
                     
                        // for each shell in this atom
                        for ( unsigned ishell = 0 ; ishell < _atomshells[rowatom].size() ; ishell++ ){
                            //boost::timer::cpu_times tstartshells = cpu_t.elapsed();
                            AOShellIterator _row = _atomshells[rowatom][ishell];
                            // for density, fill sub-part of AOatgrid
                            //ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), 0, 1);
                            ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, 0, 1, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc());
                            // (*_row)->EvalAOspace(_AOgridsub, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);

                            //ub::matrix_range< ub::matrix<double> > _gradAO = ub::subrange(gradAOgrid, 0, 3, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc());
                            //(*_row)->EvalAOGradspace(_gradAO, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);
                            (*_row)->EvalAOspace(_AOgridsub, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);
                            //boost::timer::cpu_times tendshells = cpu_t.elapsed();
                             //_t_AOvals +=  (tendshells.wall-tstartshells.wall)/1e9;

                        }  // shell in atom
                    }
                        /* ub::matrix<double> _temp     = ub::zero_matrix<double>(_blocksize[rowatom],1);
                        ub::matrix<double> _tempgrad = ub::zero_matrix<double>(_blocksize[rowatom],3);
                        
                        ub::matrix_range< ub::matrix<double> > _AOgridrow     = ub::subrange(    AOgrid, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom], 0, 1);

                         * 
                         */          
                   for ( unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size() ; sigrow++){
                    
                        // this atom
                        int rowatom = _significant_atoms[i][j][sigrow];
                        ub::matrix<double> _temp     = ub::zero_matrix<double>(1,_blocksize[rowatom]);
                                          
                        ub::matrix_range< ub::matrix<double> > _AOgridrow     = ub::subrange(    AOgrid, 0,1, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);

                        // for each atom
                        
                        for ( unsigned sigcol = 0; sigcol < _significant_atoms[i][j].size() ; sigcol++){
                            int colatom = _significant_atoms[i][j][sigcol];
                            
                            
                            // get the already calculated AO values
                            //ub::matrix_range< ub::matrix<double> >     _AOgridcol = ub::subrange(    AOgrid, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom], 0, 1);
                            //ub::matrix_range< ub::matrix<double> > _gradAOgridcol = ub::subrange(gradAOgrid, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom], 0, 3);

                            ub::matrix_range< ub::matrix<double> >     _AOgridcol = ub::subrange(    AOgrid, 0, 1, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom]);
                        
                            //ub::matrix_range< ub::matrix<double> > DMAT_here = ub::subrange( _density_matrix, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom], _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom]);
                            ub::matrix_range< ub::matrix<double> > DMAT_here = ub::subrange( _density_matrix, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom], _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);
                             
                            _temp     += ub::prod( _AOgridcol, DMAT_here);
                            
                            
                        } //col shells
                        

                        //ub::matrix_range< ub::matrix<double> > _gradAOgridrow = ub::subrange(gradAOgrid, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom], 0, 3);
                        
                        //rho_mat  += ub::prod(ub::trans(    _AOgridrow),_temp);
                        //grad_rho += ub::prod(ub::trans(_gradAOgridrow),_temp) +  ub::prod(ub::trans(_tempgrad),_AOgridrow) ;

                        rho_mat  += ub::prod(_temp, ub::trans( _AOgridrow) );
                                               
                    } // row shells 


                    _grid[i][j].grid_density  =rho_mat(0,0);
                    Density_thread[i_thread] += _grid[i][j].grid_weight * _grid[i][j].grid_density;
                    //boost::timer::cpu_times t3 = cpu_t.elapsed();
                    //_t_rho +=  (t3.wall-t0.wall)/1e9;

                } // j: for each point in atom grid
                }// each thread
            } // i: for each atom grid

             for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                 //cout << result << endl;
                result += Density_thread[i_thread]; 
                }
            density_set=true;
            return result;
         }
        
        // numerically integrate the elements of the AOOverlap matrix as check

        ub::matrix<double> NumericalIntegration::numAOoverlap(AOBasis* basis) {


            ub::matrix<double> OLMAT = ub::zero_matrix<double>(basis->_AOBasisSize, basis->_AOBasisSize);
            // for every gridpoint
            for (unsigned i = 0; i < _grid.size(); i++) {
                for (unsigned j = 0; j < _grid[i].size(); j++) {
                    // get value of orbitals at each gridpoint
                    ub::matrix<double> tmat = ub::zero_matrix<double>(basis->_AOBasisSize, 1);

                    for (vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++) {

                        ub::matrix_range< ub::matrix<double> > _submatrix = ub::subrange(tmat,0,1, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc());
                        (*_row)->EvalAOspace(_submatrix, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);
                    }

                    OLMAT += _grid[i][j].grid_weight * ub::prod( ub::trans(tmat),tmat);
                }
            } // gridpoints end


            return OLMAT;

        } // numAOoverlap
        
        
        
        
        /**
         * Computes electron density belonging to atoms listed in AtomIndeces.
         * Based on IntegrateDensity_Atomblock()
         * 
         * @param _density_matrix
         * @param basis
         * @param AtomIndeces of the atoms in the molecule
         * @return electron density
         */
        double NumericalIntegration::IntegrateDensity_Molecule(ub::matrix<double>& _density_matrix, AOBasis* basis, std::vector<int> AtomIndeces){
            
            double result=0.0;
            
            // generate a list of shells for each atom
            typedef vector< AOShell* >::iterator AOShellIterator;
            vector< vector< AOShellIterator > > _atomshells;
            vector< AOShellIterator > _singleatom;

            vector < int > _startIdx;
            vector < int > _blocksize;

            int _atomindex = 0;
            int _Idx       = 0;
            int _size      = 0;
            
            for (vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++) {
                
                
                if ( (*_row)->getIndex() == _atomindex ){
                    
                    _singleatom.push_back(_row);
                    _size += (*_row)->getNumFunc();
                    
                    
                } else {
                    
                    // append _singleatom to _atomshells
                    _atomshells.push_back(_singleatom);
                    _startIdx.push_back( _Idx );
                    _blocksize.push_back(_size);
                    // reset _singleatom
                    _singleatom.clear();
                    _size = (*_row)->getNumFunc();
                    _Idx       = (*_row)->getStartIndex();
                    _singleatom.push_back(_row);
                    _atomindex = (*_row)->getIndex();
                    
                }
                   
            }
            
            _atomshells.push_back(_singleatom);
            _startIdx.push_back( _Idx );
                    _blocksize.push_back(_size);


           
            // setup a list of min decay constants per atom
            // for every shell
            _atomindex = 0;
            double _decaymin = 1e7;
            vector< double > _minimal_decay;
            vector < vec > _positions;
            vec _localpos = (*basis->firstShell())->getPos();
            for ( vector< AOShell* >::iterator _row = basis->firstShell(); _row != basis->lastShell(); _row++   ) {
                               
                 if ( (*_row)->getIndex() == _atomindex ){
                     
                     // check all decay constants in this shell
                     for (AOShell::GaussianIterator itg = (*_row)->firstGaussian(); itg != (*_row)->lastGaussian(); itg++) {
                         AOGaussianPrimitive* gaussian = *itg;
                         double _decay = gaussian->decay;
                         if (_decay < _decaymin) {
                             _decaymin = _decay;
                         } // decay min check
                     
                     } // Gaussian Primitives 
                     
                 } else {  // if shell belongs to the actual atom
                     // add to mininal_decay vector
                     _minimal_decay.push_back(_decaymin);
                     _positions.push_back( _localpos );
                     // reset counters
                     _decaymin = 1e7;
                     _localpos = (*_row)->getPos();

                     _atomindex++;
                     
                     // check all decay constants in this shell
                     for (AOShell::GaussianIterator itg = (*_row)->firstGaussian(); itg != (*_row)->lastGaussian(); itg++) {
                         AOGaussianPrimitive* gaussian = *itg;
                         double _decay = gaussian->decay;
                         if (_decay < _decaymin) {
                             _decaymin = _decay;
                         } // decay min check
                     
                     } // Gaussian Primitives                                       
                 }
            } // all shells
                 
            // push final atom
            _minimal_decay.push_back(_decaymin);
            _positions.push_back( _localpos );
            
            
            
            // for each gridpoint, check the value of exp(-a*(r-R)^2) < 1e-10
            //                             = alpha*(r-R)^2 >~ 20.7
            
            vector< vector< vector<int> > > _significant_atoms;
            
            // each atomic grid
            for (unsigned i = 0; i < _grid.size(); i++) {
            
                vector< vector<int> > _significant_atoms_atomgrid;
                
                // each point of the atomic grid
                for (unsigned j = 0; j < _grid[i].size(); j++) {

                    vector<int> _significant_atoms_gridpoint;
                    vec grid;
                    grid.setX( _grid[i][j].grid_x);
                    grid.setY( _grid[i][j].grid_y);
                    grid.setZ( _grid[i][j].grid_z);
                    
                    // check all atoms
                    for ( unsigned iatom = 0 ; iatom < _minimal_decay.size(); iatom++){

                        vec dist = grid - _positions[iatom];
                        double distsq = dist.getX()*dist.getX() + dist.getY()*dist.getY()  + dist.getZ()*dist.getZ() ;
                        
                        // if contribution is smaller than -ln(1e-10), add atom to list
                        if ( (_minimal_decay[iatom] * distsq) < 20.7 ){
                            _significant_atoms_gridpoint.push_back(iatom);
                        }
                        
                    } // check all atoms

                    _significant_atoms_atomgrid.push_back(  _significant_atoms_gridpoint );
                   
                } // all points of this atom grid
                
                _significant_atoms.push_back(_significant_atoms_atomgrid);
               
            } // atomic grids
            
            
            
            
            
            // parallelization: distribute over threads inside one atom
            int nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif

            std::vector<double> Density_thread;
            for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){ 
                Density_thread.push_back(0.0);
            }           
            
            // for every atom
            for (unsigned i = 0; i < _grid.size(); i++) {
	      // for each point in atom grid
                
                // number of points in this atomgrid
                int atom_points = _grid[i].size();
                // divide among threads
                int atom_points_per_thread = atom_points/nthreads;
                std::vector<int> _thread_start;
                std::vector<int> _thread_stop;
                for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                    _thread_start.push_back( i_thread * atom_points_per_thread );
                    _thread_stop.push_back( (i_thread + 1) * atom_points_per_thread );
                }
                // final stop must be size
                _thread_stop[nthreads-1] = atom_points;

                
                #pragma omp parallel for
                for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                for (int j = _thread_start[i_thread]; j < _thread_stop[i_thread]; j++) {

                   // get value of orbitals at each gridpoint (vector as 1D boost matrix object -> prod )
                   ub::matrix<double> AOgrid = ub::zero_matrix<double>(1, basis->_AOBasisSize); //basis function values at grid for all significant atoms
                   ub::matrix<double> rho_mat = ub::zero_matrix<double>(1,1);
                    
		    // evaluate AO Functions for all shells, NOW BLOCKWISE

                    // for each significant atom for this grid point
                    for ( unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size() ; sigrow++){
                    
                        // this atom
                        int rowatom = _significant_atoms[i][j][sigrow];                 
                     
                        // for each shell in this atom
                        for ( unsigned ishell = 0 ; ishell < _atomshells[rowatom].size() ; ishell++ ){

                            AOShellIterator _row = _atomshells[rowatom][ishell];
                            AOShell* _shell = *_row;
                            // for density, fill sub-part of AOatgrid
                            ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, 0, 1, _shell->getStartIndex(), _shell->getStartIndex()+_shell->getNumFunc());
                            _shell->EvalAOspace(_AOgridsub, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);

                        }  // shell in atom
                    }
         
                   for ( unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size() ; sigrow++){
                    
                        // this atom
                        int rowatom = _significant_atoms[i][j][sigrow];
                        //only do things if this atom belongs to the molecule of interest and is significant
                        if(std::find(AtomIndeces.begin(), AtomIndeces.end(), rowatom) != AtomIndeces.end()){
                        
                            ub::matrix<double> _temp     = ub::zero_matrix<double>(1,_blocksize[rowatom]);

                            ub::matrix_range< ub::matrix<double> > _AOgridrow     = ub::subrange(    AOgrid, 0,1, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);

                            // for each atom

                            for ( unsigned sigcol = 0; sigcol < _significant_atoms[i][j].size() ; sigcol++){
                                int colatom = _significant_atoms[i][j][sigcol];
                                ub::matrix_range< ub::matrix<double> >     _AOgridcol = ub::subrange(    AOgrid, 0, 1, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom]);
                                ub::matrix_range< ub::matrix<double> > DMAT_here = ub::subrange( _density_matrix, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom], _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);

                                _temp     += ub::prod( _AOgridcol, DMAT_here);


                            } //col shells

                            rho_mat  += ub::prod(_temp, ub::trans( _AOgridrow) );
                        } // atom belongs to molecule
                                               
                    } // row shells 


                    _grid[i][j].grid_density  =rho_mat(0,0);
                    Density_thread[i_thread] += _grid[i][j].grid_weight * _grid[i][j].grid_density;
                } // j: for each point in atom grid
                }// each thread
            } // i: for each atom grid

             for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                 //cout << result << endl;
                result += Density_thread[i_thread]; 
                }
            density_set=true;
            
            return(result);
        }
        
        
        
        
  
        double NumericalIntegration::StupidIntegrate(std::vector<double>& _data){
            
            
            double integral = 0.0;
            int _i_point = 0;
            for ( unsigned i = 0 ; i < _grid.size(); i++){
                for ( unsigned j = 0 ; j < _grid[i].size(); j++){

                    
                    integral += _data[_i_point] * _grid[i][j].grid_weight;
                    
                    _i_point++;

                }
            }
            
            return integral;
            
        }          
        
        void NumericalIntegration::getGridpoints( ub::matrix<double>& _gridpoints ){
            
            _gridpoints = ub::zero_matrix<double>(_totalgridsize,4);
            
            int _i_point = 0;
            for ( unsigned i = 0 ; i < _grid.size(); i++){
                for ( unsigned j = 0 ; j < _grid[i].size(); j++){
                    
                    _gridpoints(_i_point,0) = _grid[i][j].grid_x;
                    _gridpoints(_i_point,1) = _grid[i][j].grid_y;
                    _gridpoints(_i_point,2) = _grid[i][j].grid_z;
                    _gridpoints(_i_point,3) = _grid[i][j].grid_weight;
                    
                    _i_point++;
                }
                }
            
            
            
            
        }
               
        void NumericalIntegration::GridSetup(string type, BasisSet* bs, vector<QMAtom*> _atoms) {
            
            const double pi = boost::math::constants::pi<double>();
            // get GridContainer
            GridContainers _grids;

            // get radial grid per element
            EulerMaclaurinGrid _radialgrid;
            _radialgrid.getRadialGrid(bs, _atoms, type, _grids); // this checks out 1:1 with NWChem results! AWESOME

         //   cout << "Radial grid summary " << endl;
           map<string, GridContainers::radial_grid>::iterator it;
         //   for (it = _grids._radial_grids.begin(); it != _grids._radial_grids.end(); ++it) {
         //       cout << " Element " << it->first << " Number of points " << it->second.radius.size() << endl;
         //  }

            // get angular grid per element
            LebedevGrid _sphericalgrid;
           // cout << "Spherical grid summary " << endl;
            for (it = _grids._radial_grids.begin(); it != _grids._radial_grids.end(); ++it) {
               _sphericalgrid.getSphericalGrid(_atoms, type, _grids);
          //     cout << " Element " << it->first << " Number of points " << _grids._spherical_grids.at(it->first).weight.size() << endl;

            }

            
            // for the partitioning, we need all inter-center distances later, stored in one-directional list
            int ij = 0;
            Rij.push_back(0.0); // 1st center "self-distance"
            Rij_mat = ub::zero_matrix<double>(_atoms.size(),_atoms.size());
            vector< QMAtom* > ::iterator ait;
            vector< QMAtom* > ::iterator bit;
            int i = 1;
            for (ait = _atoms.begin() + 1; ait != _atoms.end(); ++ait) {
                // get center coordinates in Bohr
                double x_a = (*ait)->x * tools::conv::ang2bohr;
                double y_a = (*ait)->y * tools::conv::ang2bohr;
                double z_a = (*ait)->z * tools::conv::ang2bohr;
                int j = 0;
                for (bit = _atoms.begin(); bit != ait; ++bit) {
                    ij++;
                    // get center coordinates in Bohr
                    double x_b = (*bit)->x * tools::conv::ang2bohr;
                    double y_b = (*bit)->y * tools::conv::ang2bohr;
                    double z_b = (*bit)->z * tools::conv::ang2bohr;

                    Rij.push_back(1.0 / sqrt((x_a - x_b)*(x_a - x_b) + (y_a - y_b)*(y_a - y_b) + (z_a - z_b)*(z_a - z_b)));


                    Rij_mat(i,j) = 1.0 / sqrt((x_a - x_b)*(x_a - x_b) + (y_a - y_b)*(y_a - y_b) + (z_a - z_b)*(z_a - z_b));
                                        
                    j++;
                } // atoms
                Rij.push_back(0.0); // self-distance again
                i++;
            } // atoms
            
            //cout << " Determined all inter-center distances " << endl;
            

            // combine the element-based information with the geometry
        

            int i_atom = 0;
            _totalgridsize = 0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                // get center coordinates in Bohr
                std::vector< GridContainers::integration_grid > _atomgrid;
                double x_c = (*ait)->x * tools::conv::ang2bohr;
                double y_c = (*ait)->y * tools::conv::ang2bohr;
                double z_c = (*ait)->z * tools::conv::ang2bohr;
                string name = (*ait)->type;

                // get radial grid information for this atom type
                GridContainers::radial_grid _radial_grid = _grids._radial_grids.at(name);

                
                // get spherical grid information for this atom type
                GridContainers::spherical_grid _spherical_grid = _grids._spherical_grids.at(name);

                // maximum order (= number of points) in spherical integration grid
                int maxorder = _sphericalgrid.Type2MaxOrder(name,type);
                int maxindex = _sphericalgrid.getIndexFromOrder(maxorder);

                // for pruning of integration grid, get interval boundaries for this element
                std::vector<double> PruningIntervals = _radialgrid.getPruningIntervals( name );
              //  cout << " Pruning Intervals: " << PruningIntervals[0] << " " << PruningIntervals[1] << " " << PruningIntervals[2] << " " << PruningIntervals[3] << endl;
                
                int current_order = 0;
                // get spherical grid
                std::vector<double> _theta;
                std::vector<double> _phi;
                std::vector<double> _weight;

                // for each radial value
                for (unsigned _i_rad = 0; _i_rad < _radial_grid.radius.size(); _i_rad++) {
                    double r = _radial_grid.radius[_i_rad];
                    int order;
                    // which Lebedev order for this point?
                    if ( maxindex == 1 ) {
                        // smallest possible grid anyway, nothing to do
                        order = maxorder;
                    } else if ( maxindex == 2 ) {
                        // only three intervals
                        if ( r < PruningIntervals[0] ) {
                            order = _sphericalgrid.getOrderFromIndex(1);//1;
                        } else if ( ( r >= PruningIntervals[0] ) && ( r < PruningIntervals[3] )   ){
                            order = _sphericalgrid.getOrderFromIndex(2);
                        } else {
                            order = _sphericalgrid.getOrderFromIndex(1);
                        } // maxorder == 2
                    } else {
                        // five intervals
                        if ( r < PruningIntervals[0] ) {
                            order = _sphericalgrid.getOrderFromIndex(int(2));
                        } else if ( ( r >= PruningIntervals[0]) && ( r < PruningIntervals[1] ) ) {
                            order = _sphericalgrid.getOrderFromIndex(4);
                        } else if ( ( r >= PruningIntervals[1]) && ( r < PruningIntervals[2] ) ) {
                            order = _sphericalgrid.getOrderFromIndex(max(maxindex-1, 4));
                        } else if ( (r >= PruningIntervals[2]) && ( r < PruningIntervals[3] ) ) {
                            order = maxorder;
                        } else {
                            order = _sphericalgrid.getOrderFromIndex(max(maxindex-1,1));
                        }
                    }                        


                    
                    // get new spherical grid, if order changed
                    if ( order != current_order ){
                        _theta.clear();
                        _phi.clear();
                        _weight.clear();
                        
                        _sphericalgrid.getUnitSphereGrid(order,_theta,_phi,_weight);
                        current_order = order;
                    }
                    
                    // for each (theta,phi)
                    // for (int _i_sph = 0; _i_sph < _spherical_grid.phi.size(); _i_sph++) {

                    for (unsigned _i_sph = 0; _i_sph < _phi.size(); _i_sph++) {
                        /* double p   = _spherical_grid.phi[_i_sph] * pi / 180.0; // back to rad
                        double t   = _spherical_grid.theta[_i_sph] * pi / 180.0; // back to rad
                        double ws  = _spherical_grid.weight[_i_sph];
                         */
                        double p   = _phi[_i_sph] * pi / 180.0; // back to rad
                        double t   = _theta[_i_sph] * pi / 180.0; // back to rad
                        double ws  = _weight[_i_sph];

                        double x_s = sin(p) * cos(t);
                        double y_s = sin(p) * sin(t);
                        double z_s = cos(p);


                        GridContainers::integration_grid _gridpoint;
                        _gridpoint.grid_x = x_c + r * x_s;
                        _gridpoint.grid_y = y_c + r * y_s;
                        _gridpoint.grid_z = z_c + r * z_s;

                        _gridpoint.grid_weight = _radial_grid.weight[_i_rad] * ws;

                        _atomgrid.push_back(_gridpoint);


                    } // spherical gridpoints
                } // radial gridpoint

                // cout << " Constructed full grid of atom " << i_atom << " of size " << _atomgrid.size() <<  endl;
                
                //int fullsize = _atomgrid.size();
                
                
                if ( 0 == 0 ){
                
                // now the partition function magic for this _atomgrid
                // some parameters
                //double eps = 0.002;
                double ass = 0.725;
                
                // get all distances from grid points to centers
                std::vector< std::vector<double> > rq;
                // for each center
                for (bit = _atoms.begin(); bit < _atoms.end(); ++bit) {
                    // get center coordinates
                    double x_b = (*bit)->x * tools::conv::ang2bohr;
                    double y_b = (*bit)->y * tools::conv::ang2bohr;
                    double z_b = (*bit)->z * tools::conv::ang2bohr;

                    std::vector<double> temp;
                    // for each gridpoint
                    for (std::vector<GridContainers::integration_grid >::iterator git = _atomgrid.begin(); git != _atomgrid.end(); ++git) {

                        double x = (*git).grid_x - x_b;
                        double y = (*git).grid_y - y_b;
                        double z = (*git).grid_z - z_b;

                        temp.push_back(sqrt(x * x + y * y + z * z));

                    } // gridpoint of _atomgrid
                    rq.push_back(temp); // rq[center][gridpoint]

                } // centers
                // cout << " Calculated all gridpoint distances to centers for " << i_atom << endl;
                
                // find nearest-neighbor of this atom
                double distNN = 1e10;

                vector< QMAtom* > ::iterator NNit;
                //int i_NN;

                // now check all other centers
                int i_b =0;
                for (bit = _atoms.begin(); bit != _atoms.end(); ++bit) {

                    if (bit != ait) {
                        // get center coordinates
                        double x_b = (*bit)->x * tools::conv::ang2bohr;
                        double y_b = (*bit)->y * tools::conv::ang2bohr;
                        double z_b = (*bit)->z * tools::conv::ang2bohr;

                        double distSQ = (x_c - x_b)*(x_c - x_b) + (y_c - y_b)*(y_c - y_b) + (z_c - z_b)*(z_c - z_b);

                        // update NN distance and iterator
                        if ( distSQ < distNN ) {
                            distNN = distSQ;
                            NNit = bit;
                            //i_NN = i_b;
                        }

                    } // if ( ait != bit) 
                    i_b++;
                }// bit centers
                 // cout << " Nearest neighbor of atom " << i_atom << " is atom " << i_NN << " at distance " << distNN << endl;
                
                //double radwgh = (1.0 - ass ) * sqrt(distNN) * 0.5;
                /* according to SSW scheme, all gridpoints within radwgh 
                 * of its parent center have weight one, and we can skip
                 * calculating the weighting function explicitly.
                 * Since the gridpoints in _atomgrid are sorted with increasing
                 * distance from the center, we can go through the list easily
                 */
                
             /*   int _idx_left = 0;
                for ( int i_grid  = 0 ; i_grid < _atomgrid.size(); i_grid++) {
                    if ( rq[i_atom][i_grid] > (radwgh + eps)  ) {
                        _idx_left = i_grid;
                        break; // out of the for-loop
                    }
                    i_grid++;
                } */
                
                //cout << " First forward non-unity weight is for gridpoint " << _idx_left << endl;
                
                /* Similarly, all gridpoints g for which 
                 * 
                 *      mu_ij = (r_ig - r_jg)/R_ij > a
                 *   
                 *   for i = parent atom and j = next neighbor
                 * 
                 * have zero weight. So we start from the end of the 
                 * gridpoint list and set all weights to zero until 
                 * we find the first non-zero contribution.
                 */

                // update NN distance
             /*   distNN = (ass-eps) * sqrt(distNN) ;
                // reduce "checklist" backward
                int _idx_right;
                for (int i_grid = _atomgrid.size()-1; i_grid >= _idx_left; i_grid--) {
                    cout << i_grid << "  is " <<  rq[i_atom][i_grid] - rq[i_NN][i_grid] << " vs " << distNN << endl;
                    
                    if (  (rq[i_atom][i_grid] - rq[i_NN][i_grid] ) > distNN   ) {
                        // set weight to zero
                        _atomgrid[i_grid].grid_weight = 0.0;
                    } else {
                        _idx_right = i_grid;
                        break;
                    }
                } // right index
                cout << " First backward non-zero weight is for gridpoint " << _idx_right << endl; */
                
                /* only for the remaining gridpoint [_idx_left:_idx_right], we
                 * have to evaluate the weights explicitly
                 */
                //for ( int i_grid = _idx_left; i_grid <= _idx_right ; i_grid++){
                for ( unsigned i_grid = 0; i_grid < _atomgrid.size() ; i_grid++){
                    //cout << " modifying point " << i_grid << endl;
                    // call some shit called grid_ssw0 in NWChem
                    std::vector<double> _p = SSWpartition( _atomgrid.size(), i_grid, _atoms.size(),rq, ass );
                    //cout << " partition for gridpoint " << i_grid << endl;
                    // check weight sum
                    double wsum = 0.0;
                    for (unsigned i =0 ; i < _p.size(); i++ ){
                        wsum += _p[i];
                    }
                    //cout << " sum of partition weights " << wsum << endl;
                    if ( wsum != 0.0 ){
                        
                        // update the weight of this grid point
                        _atomgrid[i_grid].grid_weight = _atomgrid[i_grid].grid_weight * _p[i_atom]/wsum;
                        //cout << " adjusting gridpoint weight "  << endl;
                    } else {
                        
                       cerr << "\nSum of partition weights of grid point " << i_grid << " of atom " << i_atom << " is zero! ";
                       throw std::runtime_error("\nThis should never happen!"); 
                        
                    }
                    

                } // partition weight for each gridpoint

                // now remove points from the grid with negligible weights
                
                for (std::vector<GridContainers::integration_grid >::iterator git = _atomgrid.begin(); git != _atomgrid.end();) {
                    if (git->grid_weight < 1e-13 ) {
                        git = _atomgrid.erase(git);
                    } else {
                        ++git;
                    }
                }
                
                } // 1 == 0
                
               // cout << " Total size of integration grid for atom: " << i_atom << " : " << _atomgrid.size() << " from " << fullsize << endl;

                _totalgridsize += _atomgrid.size() ;
                _grid.push_back(_atomgrid);
                i_atom++;
            } // atoms


            
            ofstream points;
            points.open("molgrid.xyz", ofstream::out);
            points << _totalgridsize << endl;
            points << endl;
            for ( unsigned i = 0 ; i < _grid.size(); i++){
                for ( unsigned j = 0 ; j < _grid[i].size(); j++){
                points << "X " << _grid[i][j].grid_x/tools::conv::ang2bohr << " " << _grid[i][j].grid_y/tools::conv::ang2bohr
                        << " " << _grid[i][j].grid_z/tools::conv::ang2bohr << " "  << _grid[i][j].grid_weight << endl;
                }
            }
            points.close();


        }
    
        std::vector<double> NumericalIntegration::SSWpartition(int ngrid, int igrid, int ncenters, std::vector< std::vector<double> >& rq,  double ass){
            
            // initialize partition vector to 1.0
            std::vector<double> p(ncenters,1.0);
            
            const double tol_scr = 1e-10;
            const double leps    = 1e-6; 
            // go through centers
            for ( int i = 1; i < ncenters; i++ ){
                
                int ij = i*(i+1)/2 -1; // indexing magic
                double rag = rq[i][igrid] ;
                
                // through all other centers (one-directional)
                for (int j = 0; j < i; j++ ){
                    
                    ij++;
                    if ( ( std::abs(p[i]) > tol_scr  ) || ( std::abs(p[j]) > tol_scr  ) ){
                        
                        // cout << "Rij " <<  Rij_mat(i,j) << " or " << Rij[ij] << endl;
                        
                        double mu = ( rag - rq[j][igrid] )*Rij[ij]; 
                        if ( mu > ass ) {
                            p[i] = 0.0;
                        } else if ( mu < -ass ) {
                            p[j] = 0.0;
                        } else {
                            
                            double sk;
                            if (std::abs(mu) < leps ) {
                                sk = -1.88603178008*mu + 0.5;
                            } else {
                                sk = erf1c(mu); 
                            }
                            if ( mu > 0.0 ) sk = 1.0 - sk;
                            p[j] = p[j] * sk;
                            p[i] = p[i] * (1.0-sk);
                            
                            
                        }
                        
                    }
                    
                    
                }
                
                
                
                
            }
            
            return p;
            
            
            
            
            
        }

        double NumericalIntegration::erf1c(double x){
             
            const static double alpha_erf1=1.0/0.30;
            return 0.5*erfcc((x/(1.0-x*x))*alpha_erf1);    
            
        }
              
        double NumericalIntegration::erfcc(double x){
            
            double tau = 1.0/(1.0+0.5*std::abs(x));
            
            return tau*exp(-x*x-1.26551223 + 1.00002368*tau + 0.37409196*tau*tau 
            + 0.09678418*pow(tau,3) - 0.18628806*pow(tau,4) + 0.27886807*pow(tau,5) 
            -1.13520398*pow(tau,6) + 1.48851587*pow(tau,7)  -0.82215223*pow(tau,8) 
            + 0.17087277*pow(tau,9));
            
        }
                                                                                                
    }
}
