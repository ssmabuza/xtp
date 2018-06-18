/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE sigma_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/sigma.h>


using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(sigma_test)

BOOST_AUTO_TEST_CASE(sigma_full){
  
    ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile << " H            .629118     .629118     .629118" << endl;
  xyzfile << " H           -.629118    -.629118     .629118" << endl;
  xyzfile << " H            .629118    -.629118    -.629118" << endl;
  xyzfile << " H           -.629118     .629118    -.629118" << endl;
  xyzfile.close();

  ofstream basisfile("3-21G.xml");
  basisfile <<"<basis name=\"3-21G\">" << endl;
  basisfile << "  <element name=\"H\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"8.245470e-01\">" << endl;
  basisfile << "        <contractions factor=\"9.046910e-01\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.831920e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "  <element name=\"C\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.722560e+02\">" << endl;
  basisfile << "        <contractions factor=\"6.176690e-02\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"2.591090e+01\">" << endl;
  basisfile << "        <contractions factor=\"3.587940e-01\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"5.533350e+00\">" << endl;
  basisfile << "        <contractions factor=\"7.007130e-01\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"3.664980e+00\">" << endl;
  basisfile << "        <contractions factor=\"-3.958970e-01\" type=\"S\"/>" << endl;
  basisfile << "        <contractions factor=\"2.364600e-01\" type=\"P\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"7.705450e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.215840e+00\" type=\"S\"/>" << endl;
  basisfile << "        <contractions factor=\"8.606190e-01\" type=\"P\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"1.958570e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "</basis>" << endl;
  basisfile.close();
  
  Orbitals orbitals;
  orbitals.LoadFromXYZ("molecule.xyz");
  BasisSet basis;
  basis.LoadBasisSet("3-21G.xml");
  
  AOBasis aobasis;
  aobasis.AOBasisFill(basis,orbitals.QMAtoms());
  
  Orbitals orb;
  orb.setBasisSetSize(17);
  orb.setNumberOfLevels(4,13);

Eigen::MatrixXd MOs=Eigen::MatrixXd::Zero(17,17);
MOs<<-0.00761992, -4.69664e-13, 8.35009e-15, -1.15214e-14, -0.0156169, -2.23157e-12, 1.52916e-14, 2.10997e-15, 8.21478e-15, 3.18517e-15, 2.89043e-13, -0.00949189, 1.95787e-12, 1.22168e-14, -2.63092e-15, -0.22227, 1.00844,
0.233602, -3.18103e-12, 4.05093e-14, -4.70943e-14, 0.1578, 4.75897e-11, -1.87447e-13, -1.02418e-14, 6.44484e-14, -2.6602e-14, 6.5906e-12, -0.281033, -6.67755e-12, 2.70339e-14, -9.78783e-14, -1.94373, -0.36629,
-1.63678e-13, -0.22745, -0.054851, 0.30351, 3.78688e-11, -0.201627, -0.158318, -0.233561, -0.0509347, -0.650424, 0.452606, -5.88565e-11, 0.453936, -0.165715, -0.619056, 7.0149e-12, 2.395e-14,
-4.51653e-14, -0.216509, 0.296975, -0.108582, 3.79159e-11, -0.199301, 0.283114, -0.0198557, 0.584622, 0.275311, 0.461431, -5.93732e-11, 0.453057, 0.619523, 0.166374, 7.13235e-12, 2.56811e-14,
-9.0903e-14, -0.21966, -0.235919, -0.207249, 3.75979e-11, -0.199736, -0.122681, 0.255585, -0.534902, 0.362837, 0.461224, -5.91028e-11, 0.453245, -0.453298, 0.453695, 7.01644e-12, 2.60987e-14,
0.480866, 1.8992e-11, -2.56795e-13, 4.14571e-13, 2.2709, 4.78615e-10, -2.39153e-12, -2.53852e-13, -2.15605e-13, -2.80359e-13, 7.00137e-12, 0.145171, -1.96136e-11, -2.24876e-13, -2.57294e-14, 4.04176, 0.193617,
-1.64421e-12, -0.182159, -0.0439288, 0.243073, 1.80753e-10, -0.764779, -0.600505, -0.885907, 0.0862014, 1.10077, -0.765985, 6.65828e-11, -0.579266, 0.211468, 0.789976, -1.41532e-11, -1.29659e-13,
-1.64105e-12, -0.173397, 0.23784, -0.0869607, 1.80537e-10, -0.755957, 1.07386, -0.0753135, -0.989408, -0.465933, -0.78092, 6.72256e-11, -0.578145, -0.790571, -0.212309, -1.42443e-11, -1.31306e-13,
-1.63849e-12, -0.17592, -0.188941, -0.165981, 1.79403e-10, -0.757606, -0.465334, 0.969444, 0.905262, -0.61406, -0.78057, 6.69453e-11, -0.578385, 0.578453, -0.578959, -1.40917e-11, -1.31002e-13,
0.129798, -0.274485, 0.00256652, -0.00509635, -0.0118465, 0.141392, -0.000497905, -0.000510338, -0.000526798, -0.00532572, 0.596595, 0.65313, -0.964582, -0.000361559, -0.000717866, -0.195084, 0.0246232,
0.0541331, -0.255228, 0.00238646, -0.0047388, -0.88576, 1.68364, -0.00592888, -0.00607692, -9.5047e-05, -0.000960887, 0.10764, -0.362701, 1.53456, 0.000575205, 0.00114206, -0.793844, -0.035336,
0.129798, 0.0863299, -0.0479412, 0.25617, -0.0118465, -0.0464689, 0.0750316, 0.110468, -0.0436647, -0.558989, -0.203909, 0.65313, 0.320785, 0.235387, 0.878697, -0.195084, 0.0246232,
0.0541331, 0.0802732, -0.0445777, 0.238198, -0.88576, -0.553335, 0.893449, 1.31541, -0.00787816, -0.100855, -0.0367902, -0.362701, -0.510338, -0.374479, -1.39792, -0.793844, -0.035336,
0.129798, 0.0927742, -0.197727, -0.166347, -0.0118465, -0.0473592, 0.0582544, -0.119815, -0.463559, 0.320126, -0.196433, 0.65313, 0.321765, 0.643254, -0.642737, -0.195084, 0.0246232,
0.0541331, 0.0862654, -0.183855, -0.154677, -0.88576, -0.563936, 0.693672, -1.42672, -0.0836372, 0.0577585, -0.0354411, -0.362701, -0.511897, -1.02335, 1.02253, -0.793844, -0.035336,
0.129798, 0.0953806, 0.243102, -0.0847266, -0.0118465, -0.0475639, -0.132788, 0.00985812, 0.507751, 0.244188, -0.196253, 0.65313, 0.322032, -0.87828, -0.235242, -0.195084, 0.0246232,
0.0541331, 0.088689, 0.226046, -0.0787824, -0.88576, -0.566373, -1.58119, 0.117387, 0.0916104, 0.0440574, -0.0354087, -0.362701, -0.512321, 1.39726, 0.374248, -0.793844, -0.035336;

Eigen::MatrixXd vxc=Eigen::MatrixXd::Zero(17,17);
vxc<<-0.431767, -0.131967, -1.18442e-13, -1.26466e-13, -1.02288e-13, -0.10626, -3.92543e-13, -3.95555e-13, -3.91314e-13, -0.0116413, -0.0478527, -0.0116413, -0.0478527, -0.0116413, -0.0478527, -0.0116413, -0.0478527,
-0.131967, -0.647421, 2.51812e-13, 1.39542e-13, 1.8995e-13, -0.465937, 1.53843e-14, -9.48305e-15, -5.94885e-15, -0.119833, -0.241381, -0.119833, -0.241381, -0.119833, -0.241381, -0.119833, -0.241381,
-1.18442e-13, 2.51812e-13, -0.637843, 1.33983e-13, 9.6584e-14, 5.89028e-14, -0.296161, -4.97511e-13, -5.21849e-13, -0.103175, -0.0760583, -0.103175, -0.0760583, 0.103175, 0.0760583, 0.103175, 0.0760583,
-1.26466e-13, 1.39542e-13, 1.33983e-13, -0.637843, 2.54059e-13, 4.95922e-15, -4.97536e-13, -0.296161, -4.56739e-13, -0.103175, -0.0760583, 0.103175, 0.0760583, 0.103175, 0.0760583, -0.103175, -0.0760583,
-1.02288e-13, 1.8995e-13, 9.6584e-14, 2.54059e-13, -0.637843, 2.5538e-14, -5.21859e-13, -4.56639e-13, -0.296161, -0.103175, -0.0760583, 0.103175, 0.0760583, -0.103175, -0.0760583, 0.103175, 0.0760583,
-0.10626, -0.465937, 5.89028e-14, 4.95922e-15, 2.5538e-14, -0.492236, -6.90263e-14, -8.71169e-14, -9.02027e-14, -0.180782, -0.300264, -0.180782, -0.300264, -0.180782, -0.300264, -0.180782, -0.300264,
-3.92543e-13, 1.53843e-14, -0.296161, -4.97536e-13, -5.21859e-13, -6.90263e-14, -0.375768, -4.87264e-14, -6.59106e-14, -0.147757, -0.122087, -0.147757, -0.122087, 0.147757, 0.122087, 0.147757, 0.122087,
-3.95555e-13, -9.48305e-15, -4.97511e-13, -0.296161, -4.56639e-13, -8.71169e-14, -4.87264e-14, -0.375768, -2.38269e-14, -0.147757, -0.122087, 0.147757, 0.122087, 0.147757, 0.122087, -0.147757, -0.122087,
-3.91314e-13, -5.94885e-15, -5.21849e-13, -4.56739e-13, -0.296161, -9.02027e-14, -6.59106e-14, -2.38269e-14, -0.375768, -0.147757, -0.122087, 0.147757, 0.122087, -0.147757, -0.122087, 0.147757, 0.122087,
-0.0116413, -0.119833, -0.103175, -0.103175, -0.103175, -0.180782, -0.147757, -0.147757, -0.147757, -0.571548, -0.31776, -0.00435077, -0.061678, -0.00435077, -0.061678, -0.00435077, -0.061678,
-0.0478527, -0.241381, -0.0760583, -0.0760583, -0.0760583, -0.300264, -0.122087, -0.122087, -0.122087, -0.31776, -0.353709, -0.061678, -0.149893, -0.061678, -0.149893, -0.061678, -0.149893,
-0.0116413, -0.119833, -0.103175, 0.103175, 0.103175, -0.180782, -0.147757, 0.147757, 0.147757, -0.00435077, -0.061678, -0.571548, -0.31776, -0.00435077, -0.061678, -0.00435077, -0.061678,
-0.0478527, -0.241381, -0.0760583, 0.0760583, 0.0760583, -0.300264, -0.122087, 0.122087, 0.122087, -0.061678, -0.149893, -0.31776, -0.353709, -0.061678, -0.149893, -0.061678, -0.149893,
-0.0116413, -0.119833, 0.103175, 0.103175, -0.103175, -0.180782, 0.147757, 0.147757, -0.147757, -0.00435077, -0.061678, -0.00435077, -0.061678, -0.571548, -0.31776, -0.00435077, -0.061678,
-0.0478527, -0.241381, 0.0760583, 0.0760583, -0.0760583, -0.300264, 0.122087, 0.122087, -0.122087, -0.061678, -0.149893, -0.061678, -0.149893, -0.31776, -0.353709, -0.061678, -0.149893,
-0.0116413, -0.119833, 0.103175, -0.103175, 0.103175, -0.180782, 0.147757, -0.147757, 0.147757, -0.00435077, -0.061678, -0.00435077, -0.061678, -0.00435077, -0.061678, -0.571548, -0.31776,
-0.0478527, -0.241381, 0.0760583, -0.0760583, 0.0760583, -0.300264, 0.122087, -0.122087, 0.122087, -0.061678, -0.149893, -0.061678, -0.149893, -0.061678, -0.149893, -0.31776, -0.353709;

vxc=MOs.transpose()*vxc*MOs;
AOOverlap ov;
ov.Fill(aobasis);
AOCoulomb cou;
cou.Fill(aobasis);

Eigen::VectorXd mo_energy=Eigen::VectorXd::Zero(17);
mo_energy<<-0.612601,-0.341755,-0.341755,-0.341755, 0.137304,  0.16678,  0.16678,  0.16678, 0.671592, 0.671592, 0.671592, 0.974255,  1.01205,  1.01205,  1.01205,  1.64823,  19.4429;
TCMatrix_gwbse Mmn;
Mmn.Initialize(aobasis.AOBasisSize(),0,16,0,16);
Mmn.Fill(aobasis,aobasis,MOs);
Mmn.MultiplyRightWithAuxMatrix(cou.Pseudo_InvSqrt_GWBSE(ov,1e-7));

  RPA rpa;
  rpa.configure(4,0,16);
  
  PPM ppm;
  Eigen::VectorXd screen_r=Eigen::VectorXd::Zero(1);
  screen_r(0)=ppm.getScreening_r();
  Eigen::VectorXd screen_i=Eigen::VectorXd::Zero(1);
  screen_i(0)=ppm.getScreening_i();
  rpa.setScreening(screen_r,screen_i);
  
 
  rpa.calculate_epsilon(mo_energy,Mmn);
  ppm.PPM_construct_parameters(rpa);
  Mmn.MultiplyRightWithAuxMatrix(ppm.getPpm_phi());
  
   votca::ctp::Logger _log;
  Sigma sigma=Sigma(&_log);
  sigma.configure(4,0,16,20,1e-5);
  sigma.setDFTdata(0.0,&vxc,&mo_energy);
  sigma.setGWAEnergies(mo_energy);
  sigma.CalcdiagElements(Mmn,ppm);
  Eigen::VectorXd qp_energy_ref=Eigen::VectorXd::Zero(17);
  qp_energy_ref<< -0.934164,-0.461602,-0.461602,-0.461601,0.00815998,0.323811,0.323811,0.323811,0.901732,0.901732,0.901732,1.2237,1.21009,1.21009,1.21009,1.93666,19.4256;
  bool qp_check=qp_energy_ref.isApprox(sigma.getGWAEnergies(),0.0001);
   if(!qp_check){
      cout<<"qp_energy"<<endl;
      cout<<sigma.getGWAEnergies()<<endl;
      cout<<"qp_energy_ref"<<endl;
      cout<<qp_energy_ref<<endl;
  }
  sigma.CalcOffDiagElements(Mmn,ppm);
  
  Eigen::MatrixXd Hqp=sigma.SetupFullQPHamiltonian();
  
  Eigen::MatrixXd Hqp_ref=Eigen::MatrixXd::Zero(17,17);
  Hqp_ref<<-0.934164,4.16082e-07,-1.3401e-07,1.36475e-07,0.031166,1.20677e-06,-6.81123e-07,-1.22621e-07,-1.83709e-07,1.76372e-08,-1.7807e-07,-0.0220743,-1.4977e-06,-1.75301e-06,-5.29037e-08,0.00737784,-0.00775225,
4.16082e-07,-0.461602,1.12979e-07,-1.47246e-07,-1.3086e-07,0.0443459,0.000553929,0.000427421,8.38616e-05,0.000289144,-0.0101872,-1.28339e-07,0.0141886,-0.000147938,-0.000241557,5.71202e-07,2.1119e-09,
-1.3401e-07,1.12979e-07,-0.461602,1.72197e-07,2.8006e-08,-0.000335948,0.0406153,-0.0178151,0.0101352,0.00106636,0.000113704,1.22667e-07,-0.000128128,-0.0141459,0.00111572,-4.57761e-07,5.12848e-09,
1.36475e-07,-1.47246e-07,1.72197e-07,-0.461601,-4.34283e-08,0.000614026,-0.0178095,-0.0406149,0.00106915,-0.0101316,-0.00027881,4.86348e-08,0.000252415,0.00111443,0.0141441,1.01087e-07,1.3741e-09,
0.031166,-1.3086e-07,2.8006e-08,-4.34283e-08,0.00815998,-1.70198e-07,1.14219e-07,1.10592e-09,-4.81365e-08,2.75431e-09,-2.95191e-08,-0.0166337,5.78666e-08,8.52843e-08,-1.74815e-08,-0.00112475,-0.0204625,
1.20677e-06,0.0443459,-0.000335948,0.000614026,-1.70198e-07,0.323811,1.65813e-07,-1.51122e-08,-2.98465e-05,-0.000191357,0.0138568,2.86823e-07,-0.0372319,6.58278e-05,0.000142268,-2.94575e-07,3.11298e-08,
-6.81123e-07,0.000553929,0.0406153,-0.0178095,1.14219e-07,1.65813e-07,0.323811,-6.98568e-09,-0.0120376,-0.00686446,-0.000120523,-1.7727e-07,0.000108686,0.0351664,0.0122284,1.86591e-07,-1.95807e-08,
-1.22621e-07,0.000427421,-0.0178151,-0.0406149,1.10593e-09,-1.51122e-08,-6.98568e-09,0.323811,0.00686538,-0.0120366,-0.00015138,1.6913e-07,0.000112864,-0.0122286,0.0351659,-2.32341e-08,2.57386e-09,
-1.83709e-07,8.38616e-05,0.0101352,0.00106915,-4.81365e-08,-2.98465e-05,-0.0120376,0.00686538,0.901732,6.12076e-08,-9.96554e-08,2.57089e-07,-1.03264e-05,0.00917151,-0.00170387,-3.30584e-07,-9.14928e-09,
1.76372e-08,0.000289144,0.00106636,-0.0101316,2.75431e-09,-0.000191357,-0.00686446,-0.0120366,6.12076e-08,0.901732,-2.4407e-08,-1.19304e-08,-9.06429e-05,0.00170305,0.00917133,-1.11726e-07,-6.52055e-09,
-1.7807e-07,-0.0101872,0.000113704,-0.00027881,-2.95191e-08,0.0138568,-0.000120523,-0.00015138,-9.96554e-08,-2.4407e-08,0.901732,3.23124e-07,0.00932737,2.69633e-05,8.74181e-05,-4.83481e-07,-1.90439e-08,
-0.0220743,-1.28339e-07,1.22667e-07,4.86348e-08,-0.0166337,2.86823e-07,-1.7727e-07,1.6913e-07,2.57089e-07,-1.19304e-08,3.23124e-07,1.2237,-7.31155e-07,-6.14518e-07,2.79634e-08,-0.042011,0.0229724,
-1.4977e-06,0.0141886,-0.000128128,0.000252415,5.78666e-08,-0.0372319,0.000108686,0.000112864,-1.03264e-05,-9.06429e-05,0.00932737,-7.31155e-07,1.21009,-2.99286e-07,-4.29557e-08,6.13566e-07,-7.73601e-08,
-1.75301e-06,-0.000147938,-0.0141459,0.00111443,8.52843e-08,6.58278e-05,0.0351664,-0.0122286,0.00917151,0.00170305,2.69633e-05,-6.14518e-07,-2.99286e-07,1.21009,2.02234e-07,7.00978e-07,-7.18964e-08,
-5.29037e-08,-0.000241557,0.00111572,0.0141441,-1.74815e-08,0.000142268,0.0122284,0.0351659,-0.00170387,0.00917133,8.74181e-05,2.79634e-08,-4.29557e-08,2.02234e-07,1.21009,3.77938e-08,-4.85316e-09,
0.00737784,5.71202e-07,-4.57761e-07,1.01087e-07,-0.00112475,-2.94575e-07,1.86591e-07,-2.32341e-08,-3.30584e-07,-1.11726e-07,-4.83481e-07,-0.042011,6.13566e-07,7.00978e-07,3.77938e-08,1.93666,0.0330278,
-0.00775225,2.1119e-09,5.12848e-09,1.3741e-09,-0.0204625,3.11298e-08,-1.95807e-08,2.57386e-09,-9.14928e-09,-6.52055e-09,-1.90439e-08,0.0229724,-7.73601e-08,-7.18964e-08,-4.85316e-09,0.0330278,19.4256;

  bool hqp_check=Hqp_ref.isApprox(Hqp,0.0001);
  if(!hqp_check){
      cout<<"hqp"<<endl;
      cout<<Hqp<<endl;
      cout<<"hqp_ref"<<endl;
      cout<<Hqp_ref<<endl;
  }
  
   BOOST_CHECK_EQUAL(hqp_check, true);
 BOOST_CHECK_EQUAL(qp_check, true);


}
        
BOOST_AUTO_TEST_SUITE_END()
