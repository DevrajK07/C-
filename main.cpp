//CID : 02094752
//Numerical Methods Coursework Part 2 (C++)
//N-body problem solved numerically using fourth order Runge Kutta scheme
#include "Header.h"
double G = 1.0; //Universal gravitational constant


int main(){
    
    double t0 = 0; // initial time
    
    //Creates output file
    ofstream output;
    output.open("output.txt");
    
    //Ensures file has been opened ok
    if (!output.is_open()) {
        cout << "Error opening file\n";
        return 1;
    }

    //Inputs file data from parameters.txt
    ifstream inFile("parameters.txt");

    //Checks file has been opened OK
    if (!inFile.is_open()) {
        cout << "Error opening test file\n";
        return 1;
    }

    //Initialising linecount to check for number of bodies
    int linecount = 0;
    string line;
    while (getline(inFile,line)){
        linecount++;
    }

    int n = linecount;
    //Setting no. of bodies to number of lines minus one
    n = n - 1;

    //Creating array of body data for each body in system
    body bodies[n];

    //Clears error message while reading
    inFile.clear();
    //Resets pointer to start of line
    inFile.seekg(0, ios::beg);
    //Uses getline() function to get string data from each line
    getline(inFile,line);
    double G,t_final,dt;
    //Splits each lines data using blankspace as point to split each line into multiple strings
    stringstream ss(line);
    string s1,s2,s3;
    //Obtaining data inputs from first line and assigning these values to G, t_final and dt
    getline(ss, s1, ' ');
    //Converting each string segment to a form of double (stod())
    G = stold(s1);
    getline(ss, s2, ' ');
    t_final = stold(s2);
    getline(ss, s3, ' ');
    dt = stold(s3);
    int k = 0;
    //Uses a while loop to cycle through the remaining lines and use these lines to input initial conditions for body data
    while(getline(inFile, line)){
            string myString1, myString2, myString3, myString4, myString5;
            stringstream ss(line);
            getline(ss, myString1, ' ');
            bodies[k].x = stold(myString1);
            getline(ss, myString2, ' ');
            bodies[k].y = stold(myString2);
            getline(ss, myString3, ' ');
            bodies[k].vx = stold(myString3);
            getline(ss, myString4, ' ');
            bodies[k].vy = stold(myString4);
            getline(ss, myString5, ' ');
            bodies[k].m = stold(myString5);
            k++;
    }

    //Use setw() to get nice spacing style for variables in output.txt file (makes it easy to move data into e.g. excel)
    output << "body" << setw(18) << "time" << setw(18) << "x" << setw(18) << "y" << setw(18) << "vx" << setw(18) << "vy" << "\n";

    for(int i=0; i<n; i++)
    {
    output << i << setw(20) << t0 << setw(20) << bodies[i].x << setw(20) << bodies[i].y<< setw(20) << bodies[i].vx<< setw(20) << bodies[i].vy << setw(20) <<"\n";
    }
    inFile.close();
    //loops over the time
    for (double t = t0+dt; t <= t_final+dt/2; t += dt)
    {
        //Calls Runge Kutta function to cycle through bodies and timesteps
        rungeKutta(bodies, n, t, dt,G);
        output.precision(6);
        for (int i = 0; i < n; i++)
        {
            // Prints results to output.txt file
            output << i << setw(20) << t << setw(20) << bodies[i].x << setw(20) << bodies[i].y << setw(20) << bodies[i].vx << setw(20) << bodies[i].vy << "\n";
        }
    }
    output.close();
    return 0;
}
