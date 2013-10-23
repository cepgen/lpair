#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>

using std::string;
using namespace std;

int main(int argc, char* argv[]){

  string line;
  
  int event, npart;

  int indx_pro1, stdhep_pro1, idhep_pro1, jmo1_pro1, jmo2_pro1, jda1_pro1, jda2_pro1;
  float Px_pro1, Py_pro1, Pz_pro1, E_pro1, Mass_pro1, Vx_pro1, Vy_pro1, Vz_pro1, Vt_pro1; 

  int indx_pro2, stdhep_pro2, idhep_pro2, jmo1_pro2, jmo2_pro2, jda1_pro2, jda2_pro2;
  float Px_pro2, Py_pro2, Pz_pro2, E_pro2, Mass_pro2, Vx_pro2, Vy_pro2, Vz_pro2, Vt_pro2; 

  int indx_ele1, stdhep_ele1, idhep_ele1, jmo1_ele1, jmo2_ele1, jda1_ele1, jda2_ele1;
  float Px_ele1, Py_ele1, Pz_ele1, E_ele1, Mass_ele1, Vx_ele1, Vy_ele1, Vz_ele1, Vt_ele1; 

  int indx_ele2, stdhep_ele2, idhep_ele2, jmo1_ele2, jmo2_ele2, jda1_ele2, jda2_ele2;
  float Px_ele2, Py_ele2, Pz_ele2, E_ele2, Mass_ele2, Vx_ele2, Vy_ele2, Vz_ele2, Vt_ele2; 

  int indx_centsys, stdhep_centsys, idhep_centsys, jmo1_centsys, jmo2_centsys, jda1_centsys, jda2_centsys;
  float Px_centsys, Py_centsys, Pz_centsys, E_centsys, Mass_centsys, Vx_centsys, Vy_centsys, Vz_centsys, Vt_centsys; 

  string filename;
  if (argc==2) filename = argv[1];
  else filename = "events.ascii";
  ifstream infile (filename.c_str());
  if (! infile.is_open())
  { cout << "Error opening file"; exit (1); }

  //loop over each event, 
  //reading in all the lines relevant to that event
  //then outputing the relevant info in the output file for that event
  //in each iteration of the loop
  while (! infile.eof() )
  {
    //get the first line in the event
    getline(infile,line);
    sscanf(line.c_str(),"%d %d",&event,&npart);

    //the first incoming proton (I don't need to save the kinematics)
    getline(infile,line);
    getline(infile,line);
    getline(infile,line);

    //the second incoming proton (I don't need to save the kinematics)
    getline(infile,line);
    getline(infile,line);
    getline(infile,line);


    //get the first outgoing proton (always particle 3, the daughter of first incoming proton)
    ///////////////////////////////////////////////////////
    //for inel-el events this is the proton that dissociates
    ///////////////////////////////////////////////////////
    getline(infile,line);
    sscanf(line.c_str(),"%d %d %d %d %d %d %d",
	     &indx_pro1, &idhep_pro1, &stdhep_pro1, &jmo1_pro1, &jmo2_pro1, &jda1_pro1, &jda2_pro1 );
    getline(infile,line);
    sscanf(line.c_str(),"%f %f %f %f %f",
	   &Px_pro1, &Py_pro1, &Pz_pro1, &E_pro1, &Mass_pro1 );
    getline(infile,line);
    sscanf(line.c_str(),"%f %f %f %f",
	   &Vx_pro1, &Vy_pro1, &Vz_pro1, &Vt_pro1 );


    //get the central system (always particle 4)
    getline(infile,line);
    sscanf(line.c_str(),"%d %d %d %d %d %d %d",
	     &indx_centsys, &idhep_centsys, &stdhep_centsys, &jmo1_centsys, &jmo2_centsys, &jda1_centsys, &jda2_centsys );
    getline(infile,line);
    sscanf(line.c_str(),"%f %f %f %f %f",
	   &Px_centsys, &Py_centsys, &Pz_centsys, &E_centsys, &Mass_centsys );
    getline(infile,line);
    sscanf(line.c_str(),"%f %f %f %f",
	   &Vx_centsys, &Vy_centsys, &Vz_centsys, &Vt_centsys );



    //get the second outgoing proton (always particle 5, the daughter of second incoming proton)
    getline(infile,line);
    sscanf(line.c_str(),"%d %d %d %d %d %d %d",
	     &indx_pro2, &idhep_pro2, &stdhep_pro2, &jmo1_pro2, &jmo2_pro2, &jda1_pro2, &jda2_pro2 );
    getline(infile,line);
    sscanf(line.c_str(),"%f %f %f %f %f",
	   &Px_pro2, &Py_pro2, &Pz_pro2, &E_pro2, &Mass_pro2 );
    getline(infile,line);
    sscanf(line.c_str(),"%f %f %f %f",
	   &Vx_pro2, &Vy_pro2, &Vz_pro2, &Vt_pro2 );


    //get the first outgoing electron (always particle 6)
    getline(infile,line);
    sscanf(line.c_str(),"%d %d %d %d %d %d %d",
	     &indx_ele1, &idhep_ele1, &stdhep_ele1, &jmo1_ele1, &jmo2_ele1, &jda1_ele1, &jda2_ele1 );
    getline(infile,line);
    sscanf(line.c_str(),"%f %f %f %f %f",
	   &Px_ele1, &Py_ele1, &Pz_ele1, &E_ele1, &Mass_ele1 );
    getline(infile,line);
    sscanf(line.c_str(),"%f %f %f %f",
	   &Vx_ele1, &Vy_ele1, &Vz_ele1, &Vt_ele1 );


    //get the second outgoing electron (always particle 7)
    getline(infile,line);
    sscanf(line.c_str(),"%d %d %d %d %d %d %d",
	     &indx_ele2, &idhep_ele2, &stdhep_ele2, &jmo1_ele2, &jmo2_ele2, &jda1_ele2, &jda2_ele2 );
    getline(infile,line);
    sscanf(line.c_str(),"%f %f %f %f %f",
	   &Px_ele2, &Py_ele2, &Pz_ele2, &E_ele2, &Mass_ele2 );
    getline(infile,line);
    sscanf(line.c_str(),"%f %f %f %f",
	   &Vx_ele2, &Vy_ele2, &Vz_ele2, &Vt_ele2 );



    //print the first line in the event
    printf("%8d %7d\n",event,npart);


    //print incoming proton 1
    printf("%5d %7d %4d %4d %4d %4d %4d %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n",
	     1, 2212, -1, 0, 0, 0, 0, 0.0, 0.0, 9.799995508E+02, 9.80E+02, 9.385E-01, 0.0, 0.0, 0.0, 0.0 );

    //print incoming proton 2 (the pbar) in -z direction
    printf("%5d %7d %4d %4d %4d %4d %4d %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n",
	     2, -2212, -1, 0, 0, 0, 0, 0.0, 0.0, -9.799995508E+02, 9.80E+02, 9.385E-01, 0.0, 0.0, 0.0, 0.0 );

    //print outgoing proton 1
    printf("%5d %7d %4d %4d %4d %4d %4d %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n",
	     indx_pro1, stdhep_pro1, idhep_pro1,  jmo1_pro1, jmo2_pro1, jda1_pro1, jda2_pro1, Px_pro1, Py_pro1, Pz_pro1, E_pro1, Mass_pro1, Vx_pro1, Vy_pro1, Vz_pro1, Vt_pro1 );

    //print central system
    printf("%5d %7d %4d %4d %4d %4d %4d %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n",
	     indx_centsys, stdhep_centsys, idhep_centsys, jmo1_centsys, jmo2_centsys, jda1_centsys, jda2_centsys, Px_centsys, Py_centsys, Pz_centsys, E_centsys, Mass_centsys, Vx_centsys, Vy_centsys, Vz_centsys, Vt_centsys );

    //print outgoing proton 2 (the pbar)
    printf("%5d %7d %4d %4d %4d %4d %4d %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n",
	     indx_pro2, -1*stdhep_pro2, idhep_pro2,  jmo1_pro2, jmo2_pro2, jda1_pro2, jda2_pro2, Px_pro2, Py_pro2, Pz_pro2, E_pro2, Mass_pro2, Vx_pro2, Vy_pro2, Vz_pro2, Vt_pro2 );

    //print electron 1
    printf("%5d %7d %4d %4d %4d %4d %4d %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n",
	     indx_ele1, stdhep_ele1, idhep_ele1,  jmo1_ele1, jmo2_ele1, jda1_ele1, jda2_ele1, Px_ele1, Py_ele1, Pz_ele1, E_ele1, Mass_ele1, Vx_ele1, Vy_ele1, Vz_ele1, Vt_ele1 );

    //print electron 2 
    printf("%5d %7d %4d %4d %4d %4d %4d %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n",
	     indx_ele2, stdhep_ele2, idhep_ele2,  jmo1_ele2, jmo2_ele2, jda1_ele2, jda2_ele2, Px_ele2, Py_ele2, Pz_ele2, E_ele2, Mass_ele2, Vx_ele2, Vy_ele2, Vz_ele2, Vt_ele2 );

  

  }
  return 0;



}
