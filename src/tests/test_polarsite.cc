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

#define BOOST_TEST_MODULE polararsite_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/polarsite.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(polararsite_test)

BOOST_AUTO_TEST_CASE(constructors_test) { PolarSite ps(1, "ps1"); }

BOOST_AUTO_TEST_CASE(getters_test) {
  PolarSite ps(1,"ps2");
  BOOST_CHECK_EQUAL(ps.getId(),1);
  BOOST_CHECK_EQUAL(ps.getElement(),"ps2");
}

BOOST_AUTO_TEST_CASE(multipole_test) {
  PolarSite ps(1,"ps2");
  Eigen::VectorXd multipole=Eigen::VectorXd::Zero(9);
  multipole<<1,2,3,4,8,7,2,3.3,-0.5;
  ps.setMultipole(multipole);
  bool check_mpoles=multipole.isApprox(ps.getPermMultipole(),0.0001);
   BOOST_CHECK_EQUAL(check_mpoles,true);
   
   bool check_rank=(ps.getRank()==2);
   BOOST_CHECK_EQUAL(check_rank,true);
  
}

BOOST_AUTO_TEST_CASE(translate_test) {
  PolarSite ps(1,"ps2");
  Eigen::Vector3d shift;
  shift<<0,0,5;
  ps.Translate(shift);
  BOOST_CHECK_EQUAL(shift.isApprox(ps.getPos(),1e-5),true);
}


BOOST_AUTO_TEST_CASE(rotation_test){
  PolarSite ps(1,"ps2",Eigen::Vector3d::UnitY());
  
  Eigen::Matrix3d R=Eigen::Matrix3d::Zero(); //Rotation around z axes
  R << 0, -1, 0,
      1,  0,  0,
      0,  0,  1 ;

  Eigen::VectorXd multipoles=Eigen::VectorXd::Zero(9);
  multipoles<<1,1,0,0,0,1,0,0,0; //q=1, mu_x=1 and Q_21c=1 the rest is 0
  ps.setMultipole(multipoles);
  ps.Rotate(R,Eigen::Vector3d::Zero());
bool equalpos=ps.getPos().isApprox(Eigen::Vector3d(-1,0,0),1e-5);
if(!equalpos){
  std::cout<<"Result "<<std::endl;
  std::cout<<ps.getPos()<<std::endl;
  std::cout<<"Reference"<<std::endl;
  std::cout<<Eigen::Vector3d(-1,0,0)<<std::endl;
}
BOOST_CHECK_EQUAL(equalpos,true); 


Eigen::VectorXd rotmultipoles=Eigen::VectorXd::Zero(9);
rotmultipoles<<1,0,1,0,0,0,1,0,0; //q=1, mu_y=1 and Q_21s=1 is 0
bool equalmultipoles=rotmultipoles.isApprox(ps.getPermMultipole(),1e-5);
if(!equalmultipoles){
  std::cout<<"Result "<<std::endl;
  std::cout<<ps.getPermMultipole()<<std::endl;
  std::cout<<"Reference"<<std::endl;
  std::cout<<rotmultipoles<<std::endl;
}
  BOOST_CHECK_EQUAL(equalmultipoles,true); 
}

BOOST_AUTO_TEST_CASE(interaction_test) {
  PolarSite ps1(1,"ps1");
  PolarSite ps2(2,"ps2",Eigen::Vector3d::UnitX());
  //Charge-Charge Static Interaction Test
  Eigen::VectorXd mp1 = Eigen::VectorXd::Zero(1);
  Eigen::VectorXd mp2 = Eigen::VectorXd::Zero(1);
  mp1<<1;
  mp2<<-1;
  ps1.setPolarisable(false);
  ps2.setPolarisable(false);
  ps1.setMultipole(mp1);
  ps2.setMultipole(mp2);
  
  double Energyref1=-1;
  double Energy1= ps1.InteractStatic(ps2);
  BOOST_CHECK_EQUAL(std::abs(Energy1-Energyref1)<1e-9,true); 
   // Electric Field onto polarsite test
  bool check_field=ps1.getField().isApprox(ps2.getField(),1e-5);
  if(!check_field){
    std::cout<<"Field at ps1"<<std::endl;
    std::cout<<ps1.getField()<<std::endl;
    std::cout<<"Field at ps2"<<std::endl;
    std::cout<<ps2.getField()<<std::endl;
  }
  BOOST_CHECK_EQUAL(check_field,true);
  // Potential onto polarsite test
  bool check_potential=std::abs(ps1.getPotential()+ps2.getPotential())<1e-5;
  if(!check_potential){
     std::cout<<"Potential at ps1"<<std::endl;
    std::cout<<ps1.getPotential()<<std::endl;
     std::cout<<"Potential at ps2"<<std::endl;
    std::cout<<ps2.getPotential()<<std::endl;
  }
  BOOST_CHECK_EQUAL(check_potential,true);
  
  /*    
  PolarSite ps3(3,"ps3");
  PolarSite ps4(4,"ps4",Eigen::Vector3d::UnitZ());
  Eigen::VectorXd multipole=Eigen::VectorXd::Zero(9);
  multipole<<1,2,3,4,8,7,2,3.3,-0.5;
  
  ps3.setPolarisable(true);
  ps4.setPolarisable(true);
  ps3.setMultipole(multipole);
  ps4.setMultipole(multipole);
  ps3.InteractStatic(ps4);
  */
  
  //Dipole-Dipole Interaction Test
  PolarSite ps5(5,"ps5");
  PolarSite ps6(6,"ps6",Eigen::Vector3d::UnitZ());
  Eigen::VectorXd dipole1=Eigen::VectorXd::Zero(4);
  dipole1<<0,0,0,1;  
  
  ps5.setPolarisable(false);
  ps6.setPolarisable(false);
  ps5.setMultipole(dipole1);
  ps6.setMultipole(-dipole1);
  //I put two dipoles anti-parellel to z (pag.52 Stone's book)
  
  double Energyref2 = 2;
  double Energy2 = ps5.InteractStatic(ps6);
  bool check_dipole_dipole = std::abs(Energy2-Energyref2)<1e-9;
  if(!check_dipole_dipole){
    std::cout<<"Dipole1 at"<<std::endl;
    std::cout<<ps5.getPos()<<std::endl;
    std::cout<<"Dipole2 at"<<std::endl;
    std::cout<<ps6.getPos()<<std::endl;
    std::cout<<"Energy"<<std::endl;
    std::cout<<Energy2<<std::endl;
  }
  BOOST_CHECK_EQUAL(std::abs(Energy2-Energyref2)<1e-9,true); 
  
  
  //quadrupole-quadrupole Interaction Test
  PolarSite ps7(7,"ps7");
  PolarSite ps8(8,"ps8",Eigen::Vector3d::UnitZ());
  Eigen::VectorXd multipole1=Eigen::VectorXd::Zero(9);
  multipole1<<0,0,0,0,1,0,0,0,0;
  
  ps7.setPolarisable(false);
  ps8.setPolarisable(false);
  ps7.setMultipole(multipole1);
  ps8.setMultipole(multipole1);
  //I put two quadrupoles on the z axis (see pag. 53 Stone's Book)
  double Energy3 = ps7.InteractStatic(ps8);
  double Energyref3 = 6;
  bool check_quadrupole_quadrupole = std::abs(Energy3-Energyref3)<1e-9;
  if(!check_quadrupole_quadrupole){
      std::cout<<"Quadrupole-Quadrupole interaction energy="<<Energy3<<std::endl;
  }
  BOOST_CHECK_EQUAL(check_quadrupole_quadrupole,true); 
  
  
  //dipole-quadrupole Interaction Test
  PolarSite ps9(9,"ps9");
  PolarSite ps10(8,"ps10",Eigen::Vector3d::UnitZ());
  Eigen::VectorXd dipole2=Eigen::VectorXd::Zero(4);
  dipole2<<0,0,0,1;
  Eigen::VectorXd multipole2=Eigen::VectorXd::Zero(9);
  multipole2<<0,0,0,0,1,0,0,0,0;
  
  ps9.setPolarisable(false);
  ps10.setPolarisable(false);
  ps10.setMultipole(dipole2);
  ps9.setMultipole(multipole2);
  //I put one quadrupole in the origin and one dipole on the z axis parallel to this (pag.55 Stone's book)
  double Energy4 = ps9.InteractStatic(ps10);
  double Energyref4 = -3;
  bool check_dipole_quadrupole = std::abs(Energy4-Energyref4)<1e-9;
  if(!check_dipole_quadrupole){
      std::cout<<"Dipole-Quadrupole interaction energy="<<Energy3<<std::endl;
  }
  BOOST_CHECK_EQUAL(check_dipole_quadrupole,true);    
}


BOOST_AUTO_TEST_CASE(induction_test) {
  PolarSite ps1(1,"ps1");
  PolarSite ps2(2,"ps2",Eigen::Vector3d::UnitX());
  
  Eigen::VectorXd mp1 = Eigen::VectorXd::Zero(4);
  Eigen::VectorXd mp2 = Eigen::VectorXd::Zero(4);
  mp1<<1,0,0,0;
  mp2<<1,0,0,0;
  ps1.setPolarisable(true);
  ps2.setPolarisable(true);
  ps1.setMultipole(mp1);
  ps2.setMultipole(mp2);
  
  std::cout<<"Site 1:"<<std::endl;
  std::cout<<ps1.getDipole()<<std::endl;
  std::cout<<"Site 2:"<<std::endl;
  std::cout<<ps2.getDipole()<<std::endl;
  
  Eigen::Matrix3d poltensor=Eigen::Matrix3d::Zero();
  /*
  poltensor<<2,1,0,
            1,3,1,
            0,1,2.5;
  */
  poltensor<< 1,0,0,
                0,1,0,
                  0,0,1;
  ps1.setPolarisation(poltensor);
  ps2.setPolarisation(poltensor);
  
  double Energy= ps1.InteractStatic(ps2);
  ps1.Induce(1);
  ps2.Induce(1);
  std::cout<<"Site 1:"<<std::endl;
  std::cout<<ps1.getDipole()<<std::endl;
  std::cout<<"Site 2:"<<std::endl;
  std::cout<<ps2.getDipole()<<std::endl;
  
  
  double alpha=0.39;
  double InductionEnergy=ps1.InteractInduction(ps2,alpha);
 
}




BOOST_AUTO_TEST_SUITE_END()
